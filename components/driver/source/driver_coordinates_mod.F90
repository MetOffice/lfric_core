!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief  Module to assign the values of the coordinates of the mesh to a field.
module driver_coordinates_mod

  use base_mesh_config_mod,      only: geometry,                &
                                       geometry_planar,         &
                                       geometry_spherical,      &
                                       topology,                &
                                       topology_fully_periodic, &
                                       topology_non_periodic
  use constants_mod,             only: r_def, i_def, l_def, &
                                       radians_to_degrees,  &
                                       i_halo_index, eps, pi
  use log_mod,                   only: log_event, log_scratch_space, &
                                       log_level_error
  use planet_config_mod,         only: scaled_radius
  use coord_transform_mod,       only: xyz2llr, llr2xyz, identify_panel, &
                                       xyz2alphabetar, alphabetar2xyz,   &
                                       schmidt_transform_xyz,            &
                                       inverse_schmidt_transform_xyz
  use finite_element_config_mod, only: coord_system,            &
                                       coord_system_xyz

  implicit none

  private

  public  :: assign_coordinate_field
! Make procedures public for unit testing
#ifdef UNIT_TEST
  public  :: assign_coordinate_xyz
  public  :: assign_coordinate_lonlatz
  public  :: assign_coordinate_alphabetaz
#else
  private :: assign_coordinate_xyz
  private :: assign_coordinate_lonlatz
  private :: assign_coordinate_alphabetaz
#endif

contains

  !> @brief    Subroutine which assigns the values of the coordinates of the mesh
  !!           to a field.
  !> @details  An array of size 3 for the type field is passed in to be populated.
  !!           The field proxy is used to break encapsulation and access the
  !!           function space and the data attributes of the field so that its
  !!           values can be assigned. Calls two subroutines, 'get_cell_coords'
  !!           from the mesh generator and then 'assign_coordinate' on a column by
  !!           column basis.
  !>
  !> @param[in,out] chi      Model coordinate array of size 3 of fields
  !> @param[in]     panel_id Field giving the ID of mesh panels
  !> @param[in]     mesh     Mesh on which this field is attached
  subroutine assign_coordinate_field(chi, panel_id, mesh)

    use domain_mod,            only: domain_type
    use field_mod,             only: field_type, field_proxy_type
    use reference_element_mod, only: reference_element_type
    use mesh_mod,              only: mesh_type
    use local_mesh_mod,        only: local_mesh_type
    use sci_chi_transform_mod, only: get_inverse_mesh_rotation_matrix, &
                                     get_to_rotate,                    &
                                     get_stretch_factor

    implicit none

    type( field_type ),  intent( inout )        :: chi(3)
    type( field_type ),  intent( inout )        :: panel_id
    type( mesh_type  ),  intent( in ),  pointer :: mesh

    integer(i_def),                     pointer :: map(:,:)          => null()
    integer(i_def),                     pointer :: map_pid(:,:)      => null()
    real(kind=r_def),                   pointer :: dof_coords(:,:)   => null()
    class(reference_element_type),      pointer :: reference_element => null()

    type(field_proxy_type) :: chi_proxy(3)
    type(field_proxy_type) :: panel_id_proxy
    type(domain_type)      :: domain

    real(kind=r_def) :: domain_max_x
    real(kind=r_def) :: domain_min_y

    real(kind=r_def), allocatable :: column_coords(:,:,:)
    real(kind=r_def), allocatable :: dz(:)  ! dz(nlayers) array
    real(kind=r_def), allocatable :: vertex_coords(:,:)

    integer(i_def) :: cell
    integer(i_def) :: undf, ndf, nlayers
    integer(i_def) :: undf_pid, ndf_pid, nlayers_pid
    integer(i_def) :: nverts

    integer(i_def) :: alloc_error
    integer(i_def) :: depth

    integer(i_def), allocatable :: global_dof_id(:)
    integer(i_def) :: panel_ncells
    integer(i_def) :: i
    type(local_mesh_type), pointer :: local_mesh

    logical(l_def)   :: to_rotate
    real(kind=r_def) :: inverse_rot_matrix(3,3)
    real(kind=r_def) :: stretch_factor

    ! Break encapsulation and get the proxy.
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    undf    = chi_proxy(1)%vspace%get_undf()
    ndf     = chi_proxy(1)%vspace%get_ndf()
    nlayers = chi_proxy(1)%vspace%get_nlayers()

    panel_id_proxy = panel_id%get_proxy()
    undf_pid = panel_id_proxy%vspace%get_undf()
    ndf_pid  = panel_id_proxy%vspace%get_ndf()
    nlayers_pid = panel_id_proxy%vspace%get_nlayers()

    map => chi_proxy(1)%vspace%get_whole_dofmap()
    map_pid => panel_id_proxy%vspace%get_whole_dofmap()

    allocate ( dz(nlayers), STAT = alloc_error )
    if ( alloc_error /= 0 ) then
      call log_event( " assign_coordinate_field: Unable to allocate "// &
                      "local array dz(nlayers) ", LOG_LEVEL_ERROR )
    end if

    call mesh%get_dz(dz)

    reference_element => mesh%get_reference_element()
    call reference_element%get_vertex_coordinates( vertex_coords )
    nverts = reference_element%get_number_vertices()

    allocate( column_coords(3,nverts,nlayers ) )
    dof_coords => chi_proxy(1)%vspace%get_nodes( )

    domain = mesh%get_domain()
    if ( domain%is_lonlat()) then
      domain_max_x = radians_to_degrees * domain%maximum_lonlat(axis=1)
      domain_min_y = radians_to_degrees * domain%minimum_lonlat(axis=2)
    else
      domain_max_x = domain%maximum_xy(axis=1)
      domain_min_y = domain%minimum_xy(axis=2)
    end if

    allocate( global_dof_id( panel_id_proxy%vspace%get_ndof_glob() ) )
    do i = 1,panel_id_proxy%vspace%get_ndof_glob()
      global_dof_id(i) = mesh%get_gid_from_lid(i)
    end do
    local_mesh => mesh%get_local_mesh()
    panel_ncells = local_mesh%get_ncells_global_mesh()/local_mesh%get_num_panels_global_mesh()

    ! Get *inverse* rotation matrix and stretching factor here
    stretch_factor = get_stretch_factor()
    inverse_rot_matrix = get_inverse_mesh_rotation_matrix()
    to_rotate = get_to_rotate()

    ! Throw an error if stretching factor is not 1 and not on cubed-sphere
    if ( abs(stretch_factor - 1.0_r_def) > eps .and. .not.                     &
         (geometry == geometry_spherical .and.                                 &
          topology == topology_fully_periodic) ) then
      call log_event(                                                          &
        'driver_coordinates: Cannot determine coordinates if Schmidt ' //      &
        'stretching factor is not 1 and mesh is not cubed-sphere',             &
        log_level_error                                                        &
      )
    end if

    panel_id_proxy%data = 1.0_r_def

    if ( coord_system == coord_system_xyz .or. &
         geometry == geometry_planar ) then

      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call calc_panel_id( nlayers_pid,           &
                            ndf_pid, undf_pid,     &
                            map_pid(:,cell),       &
                            panel_id_proxy%data,   &
                            global_dof_id,         &
                            panel_ncells    )

        call mesh%get_column_coords(cell,column_coords)

        call assign_coordinate_xyz( nlayers,             &
                                    ndf,                 &
                                    nverts,              &
                                    undf,                &
                                    map(:,cell),         &
                                    dz,                  &
                                    chi_proxy(1)%data,   &
                                    chi_proxy(2)%data,   &
                                    chi_proxy(3)%data,   &
                                    column_coords,       &
                                    dof_coords,          &
                                    vertex_coords,       &
                                    domain_max_x,        &
                                    domain_min_y,        &
                                    panel_id_proxy%data, &
                                    ndf_pid,             &
                                    undf_pid,            &
                                    map_pid(:,cell)      )
      end do

    else if ( geometry == geometry_spherical .and. &
              topology /= topology_fully_periodic ) then

      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call calc_panel_id( nlayers_pid,           &
                            ndf_pid, undf_pid,     &
                            map_pid(:,cell),       &
                            panel_id_proxy%data,   &
                            global_dof_id,         &
                            panel_ncells   )

        call mesh%get_column_coords(cell,column_coords)

        call assign_coordinate_lonlatz( nlayers,                 &
                                        ndf,                     &
                                        nverts,                  &
                                        undf,                    &
                                        map(:,cell),             &
                                        chi_proxy(1)%data,       &
                                        chi_proxy(2)%data,       &
                                        chi_proxy(3)%data,       &
                                        column_coords,           &
                                        dof_coords,              &
                                        vertex_coords,           &
                                        to_rotate,               &
                                        inverse_rot_matrix,      &
                                        panel_id_proxy%data,     &
                                        ndf_pid,                 &
                                        undf_pid,                &
                                        map_pid(:,cell)          )
      end do

    else if ( geometry == geometry_spherical .and. &
              topology == topology_fully_periodic ) then

      do cell = 1,chi_proxy(1)%vspace%get_ncell()

        call calc_panel_id( nlayers_pid,           &
                            ndf_pid, undf_pid,     &
                            map_pid(:,cell),       &
                            panel_id_proxy%data,   &
                            global_dof_id,         &
                            panel_ncells    )

        call mesh%get_column_coords(cell,column_coords)

        call assign_coordinate_alphabetaz( nlayers,                 &
                                           ndf,                     &
                                           nverts,                  &
                                           undf,                    &
                                           map(:,cell),             &
                                           chi_proxy(1)%data,       &
                                           chi_proxy(2)%data,       &
                                           chi_proxy(3)%data,       &
                                           column_coords,           &
                                           dof_coords,              &
                                           vertex_coords,           &
                                           to_rotate,               &
                                           inverse_rot_matrix,      &
                                           stretch_factor,          &
                                           panel_id_proxy%data,     &
                                           ndf_pid,                 &
                                           undf_pid,                &
                                           map_pid(:,cell)          )
      end do

    else
      call log_event('This coordinate system has not been implemented yet', LOG_LEVEL_ERROR)
    end if

    ! As we have correctly set the chi fields into their full halos,
    ! mark their halos as clean, out to the full halo depth
    ! This is necessary so that subsequent kernel calls don't try to
    ! halo_swap the Wchi field which is read-only
    depth = mesh%get_halo_depth()
    call chi_proxy(1)%set_clean(depth)
    call chi_proxy(2)%set_clean(depth)
    call chi_proxy(3)%set_clean(depth)

    deallocate ( dz, column_coords, vertex_coords )

  end subroutine assign_coordinate_field

  !> @brief    Assigns the cubed sphere panel ID values to the panel_id field.
  !> @details  A scalar field is passed in and all values are assigned to
  !!           be the panel IDs which are calculated from the coordinates.
  !!           For planar geometry the ID is just 1 everywhere.
  !>
  !> @param[in]   nlayers             Number of layers for the panel_id field
  !> @param[in]   ndf_pid             Number of DoFs per cell for the panel_id field
  !> @param[in]   undf_pid            Universal number of DoFs for the panel_id field
  !> @param[in]   map_pid             DoF map for the panel_id field
  !> @param[out]  panel_id            Field (to be calculated) with the ID of cubed sphere panels
  !> @param[in]   global_dof_id       Array of global id's
  !> @param[in]   panel_ncells        Number of cells per cubed sphere panel
  subroutine calc_panel_id( nlayers,            &
                            ndf_pid,            &
                            undf_pid,           &
                            map_pid,            &
                            panel_id,           &
                            global_dof_id,      &
                            panel_ncells )

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers, ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map_pid(ndf_pid)
    real(kind=r_def),    intent(out) :: panel_id(undf_pid)
    integer(kind=i_def), intent(in)  :: global_dof_id(undf_pid)
    integer(kind=i_def), intent(in)  :: panel_ncells

    ! Internal variables
    integer(kind=i_def) :: vert, k

    if ( geometry == geometry_spherical .and. &
         topology == topology_fully_periodic ) then

      ! The following code assumes that the mesh generator has ordered the
      ! global cell ids panel-by-panel. If this is ever not the case, the
      ! code will produce an incorrect panel_id.
      do k = 0, nlayers-1
        vert = map_pid(1) + k
        panel_id(vert) = real(1+int( (global_dof_id(vert)-1)/panel_ncells , i_def), r_def)
      end do

    end if

  end subroutine calc_panel_id

  !> @brief Determines and assigns the (X,Y,Z) coordinates for a single column.
  !>
  !> @param[in]   nlayers        The number of layers in the mesh
  !> @param[in]   ndf            Number of DoFs per cell for chi field space
  !> @param[in]   nverts         Number of vertices per cell
  !> @param[in]   undf           Number of universal DoFs for chi field space
  !> @param[in]   map            DoF map for chi field
  !> @param[out]  chi_1          1st coordinate field
  !> @param[out]  chi_2          2nd coordinate field
  !> @param[out]  chi_3          3rd coordinate field
  !> @param[in]   column_coords  Coordinates at mesh vertices
  !> @param[in]   chi_hat_node   Reference cell coordinates at the chi space DoFs
  !> @param[in]   chi_hat_vert   Reference cell coordinates at the cell vertices
  !> @param[in]   domain_x       Domain extent in x direction for planar mesh
  !> @param[in]   domain_y       Domain extent in y direction for planar mesh
  !> @param[in]   panel_id       Field giving IDs of mesh panels
  !> @param[in]   ndf_pid        Number of DoFs per cell for panel_id space
  !> @param[in]   undf_pid       Number of universal DoFs for panel_id space
  !> @param[in]   map_pid        DoF map for panel_id space
  subroutine assign_coordinate_xyz( nlayers,       &
                                    ndf,           &
                                    nverts,        &
                                    undf,          &
                                    map,           &
                                    dz,            &
                                    chi_1,         &
                                    chi_2,         &
                                    chi_3,         &
                                    column_coords, &
                                    chi_hat_node,  &
                                    chi_hat_vert,  &
                                    domain_x,      &
                                    domain_y,      &
                                    panel_id,      &
                                    ndf_pid,       &
                                    undf_pid,      &
                                    map_pid        )

    use reference_element_mod, only: SWB, SEB, NEB, NWB, SWT, SET, NET, NWT

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nlayers, ndf, nverts, undf
    integer(kind=i_def), intent(in)  :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map(ndf), map_pid(ndf_pid)
    real(kind=r_def),    intent(in)  :: dz(nlayers)
    real(kind=r_def),    intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def),    intent(in)  :: column_coords(3,nverts,nlayers)
    real(kind=r_def),    intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)
    real(kind=r_def),    intent(in)  :: domain_x, domain_y
    real(kind=r_def),    intent(in)  :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, vert

    real(kind=r_def)    :: interp_weight, x, y, z, radius_correction
    real(kind=r_def)    :: vertex_local_coords(3,nverts)

    radius_correction = 1.0_r_def

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
      vertex_local_coords(:,:) = column_coords(:,:,k+1)
      if ( geometry == geometry_planar .and. &
           topology /= topology_non_periodic ) then
        ! Check if point cell is on right or bottom boundary,
        ! assumes a monotonic coordinate field
        if ( column_coords(1,SEB,k+1) < column_coords(1,SWB,k+1) ) then
        ! On x boundary
          vertex_local_coords(1,SEB) = domain_x
          vertex_local_coords(1,NEB) = domain_x
          vertex_local_coords(1,SET) = domain_x
          vertex_local_coords(1,NET) = domain_x
        end if
        ! Domain does not have N-S boundaries only if topology completely periodic
        if ( column_coords(2,SWB,k+1) > column_coords(2,NWB,k+1) .and. &
             topology == topology_fully_periodic ) then
        ! On y boundary
          vertex_local_coords(2,SWB) = domain_y
          vertex_local_coords(2,SEB) = domain_y
          vertex_local_coords(2,SWT) = domain_y
          vertex_local_coords(2,SET) = domain_y
        end if
      end if

      do df = 1, ndf
        ! Compute interpolation weights
        x = 0.0_r_def
        y = 0.0_r_def
        z = 0.0_r_def
        do vert = 1,nverts
          interp_weight = &
                 (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))

          x = x + interp_weight*vertex_local_coords(1,vert)
          y = y + interp_weight*vertex_local_coords(2,vert)
          z = z + interp_weight*vertex_local_coords(3,vert)
        end do
        ! For spherical domains we need to project x,y,z back onto
        ! spherical shells
        if ( geometry == geometry_spherical ) then
          radius_correction = scaled_radius + &
                              sum(dz(1:k)) + chi_hat_node(3,df)*dz(k+1)
          radius_correction = radius_correction/sqrt(x*x + y*y + z*z)
        end if
        dfk = map(df)+k
        chi_1(dfk) = x*radius_correction
        chi_2(dfk) = y*radius_correction
        chi_3(dfk) = z*radius_correction

      end do
    end do

  end subroutine assign_coordinate_xyz

  !> @brief Determines and assigns the (alpha,beta,height) coordinates for a single column.
  !>
  !> @param[in]   nlayers            The number of layers in the mesh
  !> @param[in]   ndf                Number of DoFs per cell for chi field space
  !> @param[in]   nverts             Number of vertices per cell
  !> @param[in]   undf               Num of universal DoFs for chi field space
  !> @param[in]   map                DoF map for chi field
  !> @param[out]  chi_1              1st coordinate field
  !> @param[out]  chi_2              2nd coordinate field
  !> @param[out]  chi_3              3rd coordinate field
  !> @param[in]   column_coords      Coordinates at mesh vertices
  !> @param[in]   chi_hat_node       Reference cell coords at the chi space DoFs
  !> @param[in]   chi_hat_vert       Reference cell coords at the cell vertices
  !> @param[in]   to_rotate          Whether mesh has been rotated
  !> @param[in]   inverse_rot_matrix Rotation matrix to apply to obtain native
  !!                                 Cartesian coordinates from physical ones
  !> @param[in]   stretch_factor     Stretch factor for Schmidt transform
  !> @param[in]   panel_id           Field giving IDs of mesh panels
  !> @param[in]   ndf_pid            Number of DoFs per cell for panel_id space
  !> @param[in]   undf_pid           Number of universal DoFs for panel_id space
  !> @param[in]   map_pid            DoF map for panel_id space
  subroutine assign_coordinate_alphabetaz( nlayers,            &
                                           ndf,                &
                                           nverts,             &
                                           undf,               &
                                           map,                &
                                           chi_1,              &
                                           chi_2,              &
                                           chi_3,              &
                                           column_coords,      &
                                           chi_hat_node,       &
                                           chi_hat_vert,       &
                                           to_rotate,          &
                                           inverse_rot_matrix, &
                                           stretch_factor,     &
                                           panel_id,           &
                                           ndf_pid,            &
                                           undf_pid,           &
                                           map_pid             )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nlayers, ndf, nverts, undf
    integer(kind=i_def), intent(in)  :: map(ndf)
    real(kind=r_def),    intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def),    intent(in)  :: column_coords(3,nverts,nlayers)
    real(kind=r_def),    intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)
    logical(kind=l_def), intent(in)  :: to_rotate
    real(kind=r_def),    intent(in)  :: inverse_rot_matrix(3,3)
    real(kind=r_def),    intent(in)  :: stretch_factor
    integer(kind=i_def), intent(in)  :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)  :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, vert
    integer(kind=i_def) :: panel
    real(kind=r_def)    :: alpha, beta, radius

    real(kind=r_def) :: interp_weight
    real(kind=r_def) :: v_xyz(3), v_a, v_b, v_r
    real(kind=r_def) :: vertex_local_coords(3,nverts)

    if ( any(map_pid < 0) ) then
      write(log_scratch_space,'(A)' ) &
        'Cell map is pointing off domain'
      call log_event(log_scratch_space, log_level_error)
      return
    else
      panel = int(panel_id(map_pid(1)))
    end if

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
      vertex_local_coords(:,:) = column_coords(:,:,k+1)

      do df = 1, ndf
        dfk = map(df)+k
        ! Compute interpolation weights
        alpha = 0.0_r_def
        beta = 0.0_r_def
        radius = 0.0_r_def
        do vert = 1,nverts
          interp_weight = &
                (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
          v_xyz(:) = vertex_local_coords(:,vert)

          ! Un-rotate coordinates
          if (to_rotate) then
            v_xyz = matmul(inverse_rot_matrix, v_xyz)
          end if

          ! Unstretch coordinates
          if (abs(stretch_factor - 1.0_r_def) > EPS) then
            v_xyz = inverse_schmidt_transform_xyz(v_xyz, stretch_factor)
          end if

          call xyz2alphabetar(v_xyz(1),v_xyz(2),v_xyz(3),panel,v_a,v_b,v_r)
          alpha = alpha + interp_weight*v_a
          beta = beta + interp_weight*v_b
          radius = radius + interp_weight*v_r
        end do

        chi_1(dfk) = alpha
        chi_2(dfk) = beta

        chi_3(dfk) = radius - scaled_radius

      end do

    end do

  end subroutine assign_coordinate_alphabetaz

  !> @brief Determines and assigns the (lon,lat,h) coordinates for a single column.
  !>
  !> @param[in]   nlayers            The number of layers in the mesh
  !> @param[in]   ndf                Number of DoFs per cell for chi field space
  !> @param[in]   nverts             Number of vertices per cell
  !> @param[in]   undf               Num of universal DoFs for chi field space
  !> @param[in]   map                DoF map for chi field
  !> @param[out]  chi_1              1st coordinate field
  !> @param[out]  chi_2              2nd coordinate field
  !> @param[out]  chi_3              3rd coordinate field
  !> @param[in]   column_coords      Coordinates at mesh vertices
  !> @param[in]   chi_hat_node       Reference cell coords at the chi space DoFs
  !> @param[in]   chi_hat_vert       Reference cell coords at the cell vertices
  !> @param[in]   to_rotate          Whether mesh has been rotated
  !> @param[in]   inverse_rot_matrix Rotation matrix to apply to obtain native
  !!                                 Cartesian coordinates from physical ones
  !> @param[in]   panel_id           Field giving IDs of mesh panels
  !> @param[in]   ndf_pid            Number of DoFs per cell for panel_id space
  !> @param[in]   undf_pid           Number of universal DoFs for panel_id space
  !> @param[in]   map_pid            DoF map for panel_id space
  subroutine assign_coordinate_lonlatz( nlayers,            &
                                        ndf,                &
                                        nverts,             &
                                        undf,               &
                                        map,                &
                                        chi_1,              &
                                        chi_2,              &
                                        chi_3,              &
                                        column_coords,      &
                                        chi_hat_node,       &
                                        chi_hat_vert,       &
                                        to_rotate,          &
                                        inverse_rot_matrix, &
                                        panel_id,           &
                                        ndf_pid,            &
                                        undf_pid,           &
                                        map_pid             )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)  :: nlayers, ndf, nverts, undf
    integer(kind=i_def), intent(in)  :: map(ndf)
    real(kind=r_def),    intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def),    intent(in)  :: column_coords(3,nverts,nlayers)
    real(kind=r_def),    intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)
    logical(kind=l_def), intent(in)  :: to_rotate
    real(kind=r_def),    intent(in)  :: inverse_rot_matrix(3,3)
    integer(kind=i_def), intent(in)  :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)  :: map_pid(ndf_pid)
    real(kind=r_def),    intent(in)  :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: k, df, dfk, vert
    real(kind=r_def)    :: longitude, latitude, radius

    real(kind=r_def) :: interp_weight, max_lon, min_lon
    real(kind=r_def) :: v_xyz(3), v_lon, v_lat, v_r
    real(kind=r_def) :: vertex_local_coords(3,nverts)

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
      vertex_local_coords(:,:) = column_coords(:,:,k+1)

      min_lon = 2.0_r_def*PI
      max_lon = -2.0_r_def*PI
      do df = 1, ndf
        dfk = map(df)+k
        ! Compute interpolation weights
        longitude  = 0.0_r_def
        latitude = 0.0_r_def
        radius = 0.0_r_def
        do vert = 1,nverts
          interp_weight = &
                 (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
          v_xyz(:) = vertex_local_coords(:,vert)

          ! Un-rotate coordinates
          if (to_rotate) then
            v_xyz = matmul(inverse_rot_matrix, v_xyz)
          end if

          call xyz2llr(v_xyz(1),v_xyz(2),v_xyz(3),v_lon,v_lat,v_r)
          longitude = longitude + interp_weight*v_lon
          latitude = latitude + interp_weight*v_lat
          radius = radius + interp_weight*v_r
        end do

        chi_1(dfk) = longitude
        chi_2(dfk) = latitude
        chi_3(dfk) = radius - scaled_radius

        min_lon = MIN(longitude, min_lon)
        max_lon = MAX(longitude, max_lon)
      end do

      ! If cell spans the dateline, adjust the negative longitudes
      if ( max_lon - min_lon > PI ) then
        do df = 1, ndf
          if ( chi_1(map(df)+k) < 0.0_r_def ) then
            chi_1(map(df)+k) = chi_1(map(df)+k) + 2.0_r_def*PI
          end if
        end do
      end if
    end do

  end subroutine assign_coordinate_lonlatz

end module driver_coordinates_mod
