!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module for IO subroutines for the LFRic-XIOS interface.
!>
module lfric_xios_io_mod

  use base_mesh_config_mod,          only: geometry,           &
                                           geometry_spherical
  use clock_mod,                     only: clock_type
  use constants_mod,                 only: dp_xios,                       &
                                           i_def, i_halo_index, i_native, &
                                           r_def, r_second, l_def,        &
                                           radians_to_degrees,            &
                                           str_def
  use coord_transform_mod,           only: xyz2llr
  use mesh_constructor_helper_functions_mod, &
                                     only: domain_size_type
  use field_mod,                     only: field_type, field_proxy_type
  use finite_element_config_mod,     only: element_order
  use lfric_xios_clock_mod,          only: lfric_xios_clock_type
  use lfric_xios_context_mod,        only: lfric_xios_context_type
  use lfric_xios_file_mod,           only: xios_file_type
  use function_space_mod,            only: function_space_type, BASIS
  use function_space_collection_mod, only: function_space_collection
  use fs_continuity_mod,             only: W0, W1, W2, W3, Wtheta, W2H, &
                                           name_from_functionspace
  use io_context_mod,                only: io_context_type, &
                                           io_context_initialiser_type
  use linked_list_mod,               only: linked_list_type, &
                                           linked_list_item_type
  use log_mod,                       only: log_event,       &
                                           log_level_error, &
                                           log_level_info,  &
                                           log_scratch_space
  use mesh_mod,                      only: mesh_type
  use mpi_mod,                       only: get_comm_size, &
                                           get_comm_rank, &
                                           all_gather
  use psykal_lite_mod,               only: invoke_nodal_xyz_coordinates_kernel, &
                                           invoke_nodal_coordinates_kernel
  use psykal_builtin_light_mod,      only: invoke_pointwise_convert_xyz2llr
  use xios,                          only: xios_duration,      &
                                           xios_fieldgroup,    &
                                           xios_file,          &
                                           xios_get_attr,      &
                                           xios_get_handle,    &
                                           xios_set_attr,      &
                                           xios_set_axis_attr, &
                                           xios_set_domain_attr

  implicit none
  private
  public :: initialise_xios, populate_filelist_if

  interface
    !> Callback for model to list the files it wants to access.
    !>
    !> @param[in,out] file_list  Collection of xios_file_type objects.
    !> @param[in]     clock      Model's time object.
    !>
    subroutine populate_filelist_if( file_list, clock )
      import clock_type, linked_list_type
      implicit none
      class(linked_list_type), intent(inout) :: file_list
      class(clock_type),       intent(in)    :: clock
    end subroutine populate_filelist_if
  end interface

  !> Sets up context with XIOS specifics.
  !>
  type, extends(io_context_initialiser_type) :: setup_xios_type
    private
    integer           :: mesh_id
    integer           :: twod_mesh_id
    class(field_type), &
      pointer         :: chi(:)            => null()
    class(field_type), &
      pointer         :: panel_id          => null()
    procedure(populate_filelist_if), &
      pointer, nopass :: populate_filelist => null()
    !
    ! An unused allocatable integer that prevents an intenal compiler error
    ! with the GNU Fortran compiler. Adding an allocatable forces the compiler
    ! to accept that the object has a finaliser. It gets confused without it.
    !
    !> @todo This is a workaround for GCC bug id 61767 - when this bug is
    !>       fixed, the integer can be removed.
    !>
    integer(kind=i_def), allocatable :: dummy_for_gcc
  contains
    private
    procedure, public :: initialise => initialise_setup_xios
    procedure, public :: callback   => callback_setup_xios
  end type setup_xios_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Null file list populator. Used when population is not needed.
  !>
  !> @param[in,out] filelist  List to populate with xios_file_type objects.
  !> @param[in]     clock     Model time.
  !>
  subroutine dummy_populate( filelist, clock )
    implicit none
    class(linked_list_type), intent(inout) :: filelist
    class(clock_type),       intent(in)    :: clock
  end subroutine dummy_populate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialises a setup_xios object.
  !>
  !> @param[in] mesh_id            Identifier for the primary mesh.
  !> @param[in] twod_mesh_id       Identifier for the primary 2D mesh.
  !> @param[in] chi                Co-ordinate fields.
  !> @param[in] panel_id           Field with IDs of mesh panels.
  !> @param[in] populate_filelist  Model procedure to fill file list.
  !>
  subroutine initialise_setup_xios( this,         &
                                    mesh_id,      &
                                    twod_mesh_id, &
                                    chi,          &
                                    panel_id,     &
                                    populate_filelist )

    implicit none

    class(setup_xios_type), intent(inout)       :: this
    integer(kind=i_def),    intent(in)          :: mesh_id
    integer(kind=i_def),    intent(in)          :: twod_mesh_id
    class(field_type),      intent(in), target  :: chi(:)
    class(field_type),      intent(in), target  :: panel_id
    procedure(populate_filelist_if), &
                            intent(in), pointer :: populate_filelist

    this%mesh_id           = mesh_id
    this%twod_mesh_id      = twod_mesh_id
    this%chi               => chi
    this%panel_id          => panel_id
    this%populate_filelist => populate_filelist

  end subroutine initialise_setup_xios

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initialises context.
  !>
  !> @param[in] context  Freshly created context object to be initialised.
  !>
  subroutine callback_setup_xios( this, context )

    implicit none

    class(setup_xios_type), intent(inout) :: this
    class(io_context_type), intent(in)    :: context

    type(linked_list_type) :: file_list

    file_list = linked_list_type()
    select type(context)
      class is (lfric_xios_context_type)
        call init_xios_dimensions(this%mesh_id, this%twod_mesh_id, this%chi, this%panel_id)
        call this%populate_filelist( file_list, context%get_clock() )
        call setup_xios_files( file_list, context%get_clock() )

      class default
          write( log_scratch_space, '(A)' ) &
            "XIOS setup callback was passed the context for something else."
          call log_event( log_scratch_space, log_level_error )
    end select

  end subroutine callback_setup_xios

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs XIOS context, dimension and file initialisation.
  !>
  !> @param[out]  context            I/O context to be initialised
  !> @param[in]   identifier         The XIOS context ID
  !> @param[in]   communicator       MPI communicator used by the XIOS context
  !> @param[in]   mesh_id            Mesh ID
  !> @param[in]   twod_mesh_id       2D Mesh ID
  !> @param[in]   chi                Coordinate field
  !> @param[in]   panel_id           Field containing IDs of mesh panels
  !> @param[in]   start_step         The starting model timestep
  !> @param[in]   end_step           The final model timestep
  !> @param[in]   spinup_period      The number of timesteps taken for spin up
  !> @param[in]   seconds_per_step   Seconds per model timestep
  !> @param[in]   timer_flag         Flag for use of subroutine timers
  !> @param[in]   populate_filelist  Optional procedure to build file list
  !>
  subroutine initialise_xios( context,          &
                              identifier,       &
                              communicator,     &
                              mesh_id,          &
                              twod_mesh_id,     &
                              chi,              &
                              panel_id,         &
                              start_step,       &
                              end_step,         &
                              spinup_period,    &
                              seconds_per_step, &
                              timer_flag,       &
                              populate_filelist )

    implicit none

    class(io_context_type),        intent(out), allocatable :: context
    character(len=*),              intent(in)               :: identifier
    integer(kind=i_def),           intent(in)               :: communicator
    integer(kind=i_def),           intent(in)               :: mesh_id
    integer(kind=i_def),           intent(in)               :: twod_mesh_id
    class(field_type),             intent(in)               :: chi(:)
    class(field_type),             intent(in)               :: panel_id
    character(len=*),              intent(in)               :: start_step
    character(len=*),              intent(in)               :: end_step
    real(r_second),                intent(in)               :: spinup_period
    real(r_second),                intent(in)               :: seconds_per_step
    logical(kind=l_def), optional, intent(in)               :: timer_flag
    procedure(populate_filelist_if), &
                optional, pointer, intent(in)               :: populate_filelist

    procedure(populate_filelist_if), pointer :: dummy_pointer

    type(setup_xios_type) :: callback

    integer :: rc

    if (present (populate_filelist)) then
      call callback%initialise( mesh_id, twod_mesh_id, chi, panel_id, populate_filelist )
    else
      dummy_pointer => dummy_populate
      call callback%initialise( mesh_id, twod_mesh_id, chi, panel_id, dummy_pointer )
    end if

    allocate( lfric_xios_context_type::context, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate XIOS context.", log_level_error )
    end if
    select type(context)
      class is (lfric_xios_context_type)
        call context%initialise( identifier,       &
                                 communicator,     &
                                 callback,         &
                                 start_step,       &
                                 end_step,         &
                                 spinup_period,    &
                                 seconds_per_step, &
                                 timer_flag )
      ! No need to default the select as we allocated just above.
    end select

  end subroutine initialise_xios

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief    Performs XIOS domain and axis initialisation.
  !> @details  Calculates the coordinates and bounds for the different kinds
  !!           of XIOS dimensionality (domains, axes, etc) and initialised the
  !!           corresponding XIOS objects.
  !>
  !> @param[in]  mesh_id       Mesh id
  !> @param[in]  twod_mesh_id  2D Mesh id
  !> @param[in]  chi           Coordinate field
  !> @param[in]  panel_id      Field with IDs of mesh panels
  !>
  subroutine init_xios_dimensions(mesh_id, twod_mesh_id, chi, panel_id)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: mesh_id
    integer(kind=i_def), intent(in) :: twod_mesh_id
    type(field_type),    intent(in) :: chi(:)
    type(field_type),    intent(in) :: panel_id

    ! Local variables
    integer(kind=i_def) :: i

    ! Node domain (W0)
    integer(kind=i_def)           :: coord_dim_full
    integer(kind=i_def)           :: coord_dim_owned
    real(kind=r_def), allocatable :: nodes_lon_full(:)
    real(kind=r_def), allocatable :: nodes_lat_full(:)
    real(dp_xios),    allocatable :: nodes_lon(:)
    real(dp_xios),    allocatable :: nodes_lat(:)
    real(dp_xios),    allocatable :: bnd_nodes_lon(:,:)
    real(dp_xios),    allocatable :: bnd_nodes_lat(:,:)

    ! Face domain (W3)
    real(dp_xios),allocatable :: bnd_faces_lon(:,:)
    real(dp_xios),allocatable :: bnd_faces_lat(:,:)

    ! Edge domain on half levels (W2H)
    real(dp_xios),allocatable :: bnd_edges_lon(:,:)
    real(dp_xios),allocatable :: bnd_edges_lat(:,:)

    ! Levels variables
    integer(kind=i_def) :: nfull_levels

    ! Checkpoint domain parameters
    character(len=str_def)            :: domain_name, domain_fs_name
    integer(kind=i_native), parameter :: domain_function_spaces(5) &
                                          = (/W0, W1, W2, W3, Wtheta/)
    integer(kind=i_native) :: fs_index

    ! Variables needed to compute output domain coordinates in lat-long
    type( field_type ) :: sample_chi(3)
    type( field_type ) :: coord_output(3)
    ! Field proxies (to calculate domain coordinate info)
    type(field_proxy_type), target  :: proxy_coord_output(3)


    ! Variables for local and global mesh information
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def)      :: num_face_local
    integer(kind=i_def)      :: nodes_per_edge
    integer(kind=i_def)      :: nodes_per_face
    integer(kind=i_def)      :: num_edge_local

    type(function_space_type), pointer :: output_field_fs   => null()
    type(function_space_type), pointer :: w2h_fs   => null()

    ! Variables for the gather to determine global domain sizes
    ! from the local partitioned ones
    integer(kind=i_def), allocatable :: local_undf(:)
    integer(kind=i_def) :: use_i_index(2) = (/ W3, Wtheta /)

    ! Factor to convert coords from radians to degrees if needed
    ! set as 1.0 for planar mesh
    real(kind=r_def) :: r2d

    if ( geometry == geometry_spherical ) then
     r2d = radians_to_degrees
    else
     r2d = 1.0_r_def
    end if

    ! Set up array to hold number of dofs for local domains
    allocate(local_undf(1))

    ! Set up fields to hold the output coordinates
    output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W0 )
    do i = 1,3
      call coord_output(i)%initialise( vector_space = output_field_fs )
      proxy_coord_output(i) = coord_output(i)%get_proxy()
    end do

    ! Get mesh information
    mesh => coord_output(1)%get_mesh()
    num_face_local = mesh%get_last_edge_cell()
    nodes_per_face = mesh%get_nverts_per_cell_2d()
    nodes_per_edge = mesh%get_nverts_per_edge()

    ! Calculate the local size of a W2H fs in order to determine
    ! how many edge dofs for the current partition
    w2h_fs => function_space_collection%get_fs( mesh_id, element_order, W2H )
    num_edge_local = w2h_fs%get_last_dof_owned()/size(w2h_fs%get_levels())

    ! Get the local value for last owned dof
    local_undf(1) = proxy_coord_output(1)%vspace%get_last_dof_owned()

    ! Get the unique fractional levels to set up vertical output domain
    nfull_levels = size( proxy_coord_output(1)%vspace%get_levels() )

    ! Get sizes of full nodal coordinate field as well as local partition size
    coord_dim_full = size(proxy_coord_output(1)%data) / nfull_levels
    coord_dim_owned = local_undf(1) / nfull_levels

    ! Obtain sample_chi, which will be used for setting up XIOS coordinates
    if ( geometry == geometry_spherical ) then
      ! Sample chi on W0 function space to prevent "unzipping" of cubed-sphere mesh
      do i = 1,3
        call sample_chi(i)%initialise( vector_space = output_field_fs )
      end do
      ! Convert to (X,Y,Z) coordinates
      call invoke_nodal_xyz_coordinates_kernel(sample_chi, chi, panel_id)
    else
      ! For planar geometries just re-use existing chi which are already (X,Y,Z)
      do i = 1,3
        call chi(i)%copy_field(sample_chi(i))
      end do
    end if

    ! Allocate coordinate arrays
    allocate(nodes_lon_full(coord_dim_full))
    allocate(nodes_lat_full(coord_dim_full))

    allocate(nodes_lon( coord_dim_owned ))
    allocate(nodes_lat( coord_dim_owned ))

    allocate(bnd_nodes_lon(1,size(nodes_lon)))
    allocate(bnd_nodes_lat(1,size(nodes_lat)))

    allocate(bnd_faces_lon(nodes_per_face,num_face_local))
    allocate(bnd_faces_lat(nodes_per_face,num_face_local))

    allocate(bnd_edges_lon(nodes_per_edge,num_edge_local))
    allocate(bnd_edges_lat(nodes_per_edge,num_edge_local))

    ! Calculate the node coords arrays and also the face and edge bounds
    call calc_xios_domain_coords(mesh, coord_output, sample_chi, &
                                 nfull_levels, num_face_local,   &
                                 nodes_lon_full, nodes_lat_full, &
                                 bnd_faces_lon, bnd_faces_lat,   &
                                 bnd_edges_lon, bnd_edges_lat)

    ! Construct node bounds arrays
    bnd_nodes_lon=(reshape(nodes_lon, (/1, size(nodes_lon)/) ) )
    bnd_nodes_lat=(reshape(nodes_lat, (/1, size(nodes_lat)/) ) )

    ! Initialise XIOS UGRID domains
    call init_xios_ugrid_domain( "node", mesh_id, W0,  sample_chi, bnd_nodes_lon, bnd_nodes_lat )
    call init_xios_ugrid_domain( "face", mesh_id, W3,  sample_chi, bnd_faces_lon, bnd_faces_lat )
    call init_xios_ugrid_domain( "edge", mesh_id, W2H, sample_chi, bnd_edges_lon, bnd_edges_lat )

    ! Initialise XIOS axes
    call init_xios_axis( "vert_axis_full_levels", mesh_id, W0 )
    call init_xios_axis( "vert_axis_half_levels", mesh_id, W3 )

    ! Create all the regular checkpoint domains based on current function spaces
    ! Loop over function spaces we need to create domains for:
    do fs_index = lbound(domain_function_spaces, 1), &
                  ubound(domain_function_spaces, 1)

      domain_fs_name = name_from_functionspace(domain_function_spaces(fs_index))
      domain_name = "checkpoint_" // trim(domain_fs_name)

      ! Enable use of the XIOS i_index for relevant function spaces
      if (any( use_i_index == domain_function_spaces(fs_index) )) then
        call checkpoint_domain_init(domain_function_spaces(fs_index), &
                                    trim(domain_name), mesh_id, sample_chi, .true.)
      else
        call checkpoint_domain_init(domain_function_spaces(fs_index), &
                                    trim(domain_name), mesh_id, sample_chi, .false.)
      end if

    end do

    ! Clean up things that are not needed after dimension setup
    if ( allocated(bnd_nodes_lon) ) deallocate(bnd_nodes_lon)
    if ( allocated(bnd_nodes_lat) ) deallocate(bnd_nodes_lat)
    if ( allocated(bnd_edges_lon) ) deallocate(bnd_edges_lon)
    if ( allocated(bnd_edges_lat) ) deallocate(bnd_edges_lat)
    if ( allocated(bnd_faces_lon) ) deallocate(bnd_faces_lon)
    if ( allocated(bnd_faces_lat) ) deallocate(bnd_faces_lat)

    return

  end subroutine init_xios_dimensions

  !> @brief  Sets up XIOS file context information from list of file objects
  !>
  !> @param[in]  files_list  List of file objects
  !> @param[in]  clock       Clock object
  !>
  subroutine setup_xios_files(files_list, clock)

    implicit none

    type(linked_list_type), intent(in) :: files_list
    type(clock_type),       intent(in) :: clock

    type(xios_file)        :: file_hdl
    type(xios_duration)    :: file_freq
    type(xios_fieldgroup)  :: field_group_hdl
    character(len=str_def) :: field_group_id
    character(len=str_def) :: file_mode

    type(linked_list_item_type), pointer :: loop  => null()
    type(xios_file_type),        pointer :: file  => null()

    ! Start at the head of the time_axis linked list
    loop => files_list%get_head()
    do
      ! If list is empty or we're at the end of list, return a null pointer
      if ( .not. associated(loop) ) then
        nullify(file)
        exit
      end if

      ! tmp_ptr is a dummy pointer used to 'cast' to the xios_file_type so that
      ! we can get at the information in the payload
      select type( tmp_ptr => loop%payload )
        type is (xios_file_type)
          file => tmp_ptr

          ! Get file handle from XIOS and set attributes
          call xios_get_handle( file%get_xios_id(), file_hdl )
          call xios_set_attr( file_hdl, name=file%get_path() )

          ! Set XIOS duration object second value equal to file output frequency
          if ( .not. file%get_output_freq() == -999 ) then
            file_freq%second = file%get_output_freq() * clock%get_seconds_per_step()
            call xios_set_attr( file_hdl, output_freq=file_freq )
          end if

          call xios_set_attr( file_hdl, enabled=.true. )

          ! If there is an associated field group, enable it
          field_group_id = file%get_field_group()

          if ( .not. field_group_id == "unset" ) then
            call xios_get_handle( field_group_id, field_group_hdl )
            call xios_set_attr( field_group_hdl, enabled=.true. )
          end if

          ! If file is not in "read" mode switch time-counter to exclusive
          call xios_get_attr( file_hdl, mode=file_mode )
          if ( .not. file_mode == "read" ) then
            call xios_set_attr( file_hdl, time_counter="exclusive", &
                                          time_counter_name="time" )
          end if

      end select
      loop => loop%next
    end do

    nullify(loop)
    nullify(file)

  end subroutine setup_xios_files

  !> @brief   Compute the node domain coords for this partition
  !> @details Samples the chi field at nodal points, calculates cartesian coordinates.
  !>          For spherical geometry, converts to lat-lon in degrees for specified layer
  !>
  !> @param[in]     mesh            The id of the partitioned mesh
  !> @param[in]     nodal_coords          Input field
  !> @param[in]     chi                   Input coordinate field
  !> @param[in]     nlayers               The number of layers data is output on
  !> @param[in]     ncells                The number of cells on the partition
  !> @param[out]    lon_coords            Array of longitude coordinates for the nodes
  !> @param[out]    lat_coords            Array of latitude coordinates for the nodes
  !> @param[inout]  face_bnds_lon_coords  Array of longitude coords making up the faces
  !> @param[inout]  face_bnds_lat_coords  Array of latitude coords making up the faces
  !> @param[inout]  edge_bnds_lon_coords  Array of coords making up the edges
  !> @param[inout]  edge_bnds_lat_coords  Array of coords making up the edges
  !>
  subroutine calc_xios_domain_coords(mesh, nodal_coords, chi, &
                                     nlayers, ncells,               &
                                     lon_coords, lat_coords,        &
                                     face_bnds_lon_coords,          &
                                     face_bnds_lat_coords,          &
                                     edge_bnds_lon_coords,          &
                                     edge_bnds_lat_coords)

    implicit none

    type(mesh_type), pointer, intent(in)    :: mesh
    type(field_type),         intent(in)    :: nodal_coords(3)
    type(field_type),         intent(in)    :: chi(:)
    integer(kind=i_def),      intent(in)    :: nlayers
    integer(kind=i_def),      intent(in)    :: ncells
    real(kind=r_def),         intent(out)   :: lon_coords(:), lat_coords(:)
    real(kind=dp_xios),       intent(inout) :: face_bnds_lon_coords(:,:)
    real(kind=dp_xios),       intent(inout) :: face_bnds_lat_coords(:,:)
    real(kind=dp_xios),       intent(inout) :: edge_bnds_lon_coords(:,:)
    real(kind=dp_xios),       intent(inout) :: edge_bnds_lat_coords(:,:)

    type(field_proxy_type) :: x_p(3), chi_p(3)

    integer(kind=i_def) :: cell, edge_count
    integer(kind=i_def) :: ndf_chi, ndf_x
    integer(kind=i_def) :: dim_chi
    integer(kind=i_def) :: df_x, df_chi, i
    integer(kind=i_def) :: edge1, edge2
    real(kind=r_def)    :: xyz(3)
    real(kind=r_def)    :: llr(3)

    integer(kind=i_def), pointer :: map_chi(:)   => null()
    integer(kind=i_def), pointer :: map_x(:)     => null()
    real(kind=r_def),    pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    type(domain_size_type) :: domain_size
    real(r_def), parameter :: vertex_location_tol = 1e-8

    ! Factor to convert coords from radians to degrees if needed
    ! set as 1.0 for planar mesh
    real(kind=r_def) :: r2d

    edge_count = 0

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    ndf_x  = x_p(1)%vspace%get_ndf( )
    nodes_x => x_p(1)%vspace%get_nodes()
    ndf_chi  = chi_p(1)%vspace%get_ndf( )

    dim_chi = chi_p(1)%vspace%get_dim_space( )

    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))

    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
      call chi_p(1)%halo_exchange(depth=1)
    end if

    if (chi_p(2)%is_dirty(depth=1)) then
      call chi_p(2)%halo_exchange(depth=1)
    end if

    if (chi_p(3)%is_dirty(depth=1)) then
      call chi_p(3)%halo_exchange(depth=1)
    end if

    ! Get domains size to determine boundary vertex locations
    domain_size = mesh%get_domain_size()

    ! Loop over cells
    do cell = 1, ncells

      map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
      map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )

      ! Loop over bottom half of the cell dofmap for the given layer
      do df_x = 1,(ndf_x/2)
        xyz(:) = 0.0_r_def
        do df_chi = 1, (ndf_chi/2)
          xyz(1) = xyz(1) + chi_p(1)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
          xyz(2) = xyz(2) + chi_p(2)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
          xyz(3) = xyz(3) + chi_p(3)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
        end do

        ! Convert to lat-lon in degrees if required
        if ( geometry == geometry_spherical ) then
          r2d = radians_to_degrees
          call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))

          lon_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(1)*r2d
          lat_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(2)*r2d

          face_bnds_lon_coords(df_x,cell) = llr(1)*r2d
          face_bnds_lat_coords(df_x,cell) = llr(2)*r2d
        else
          r2d = 1.0_r_def

          lon_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(1)*r2d
          lat_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(2)*r2d

          ! For biperiodic meshes, ensure that nodes on mesh boundary have
          ! consistent coordinates
          if ( mesh%is_topology_periodic() ) then
            if (abs(lon_coords( ((map_x(df_x)-1)/nlayers)+1) - domain_size%minimum%x) &
                  < vertex_location_tol) then
              lon_coords( ((map_x(df_x)-1)/nlayers)+1) = domain_size%maximum%x
            end if
            if (abs(lat_coords( ((map_x(df_x)-1)/nlayers)+1) - domain_size%minimum%y) &
                  < vertex_location_tol) then
              lat_coords( ((map_x(df_x)-1)/nlayers)+1) = domain_size%maximum%y
            end if
          end if

          face_bnds_lon_coords(df_x,cell) = lon_coords( ((map_x(df_x)-1)/nlayers)+1)
          face_bnds_lat_coords(df_x,cell) = lat_coords( ((map_x(df_x)-1)/nlayers)+1)

        end if
      end do ! Loop over bottom layer dofs

      ! For this cell compute the edge-bounds coordinates from the face-bounds coordinates
      do df_x = 1,(ndf_x/2)

        ! Retrieve the lat / lon coords of the points bounding the edge
        if (df_x == ndf_x/2) then
          edge1 = df_x
          edge2 = 1
        else
          edge1 = df_x
          edge2 = df_x + 1
        end if

        ! Is the edge owned by this cell?
        if (mesh%get_edge_cell_owner(df_x, cell) == cell) then
          edge_count = edge_count + 1

          edge_bnds_lon_coords(1,edge_count) = face_bnds_lon_coords(edge1,cell)
          edge_bnds_lon_coords(2,edge_count) = face_bnds_lon_coords(edge2,cell)
          edge_bnds_lat_coords(1,edge_count) = face_bnds_lat_coords(edge1,cell)
          edge_bnds_lat_coords(2,edge_count) = face_bnds_lat_coords(edge2,cell)
        end if ! Edge is owned by this cell

      end do ! loop over edges

    end do ! loop over cells

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)

    nullify( map_chi, map_x, nodes_x )

  end subroutine calc_xios_domain_coords

  !> @brief    Performs XIOS checkpoint domain initialisation
  !>
  !> @param[in]  fs_id        Function space id
  !> @param[in]  domain_name  XIOS domain name
  !> @param[in]  mesh_id      Mesh id
  !> @param[in]  chi          Coordinate field
  !> @param[in]  use_index    Flag to specify use of domain index
  !>                          to preserve order over decomposition
  !> @param[in]  k_order      Function space order (optional,
  !>                          default = 0)
  !>
  subroutine checkpoint_domain_init(fs_id, domain_name, mesh_id, chi, &
                                         use_index, k_order)

    implicit none

    ! Arguments
    integer(kind=i_def),           intent(in) :: fs_id
    character(len=*),              intent(in) :: domain_name
    integer(kind=i_def),           intent(in) :: mesh_id
    type(field_type),              intent(in) :: chi(3)
    logical(kind=l_def),           intent(in) :: use_index
    integer(kind=i_def), optional, intent(in) :: k_order

    ! Local variables
    integer(kind=i_def) :: i
    integer(kind=i_def) :: k_ord

    ! Checkpoint domain
    integer(kind=i_def)                     :: ibegin_checkpoint
    real(dp_xios), allocatable              :: checkpoint_lon(:)
    real(dp_xios), allocatable              :: checkpoint_lat(:)
    real(dp_xios), allocatable              :: bnd_checkpoint_lon(:,:)
    real(dp_xios), allocatable              :: bnd_checkpoint_lat(:,:)
    integer(kind=i_halo_index), allocatable :: domain_index(:)


    ! Variables needed to compute output domain coordinates in lat-long
    type( field_type ) :: coord_output(3)
    type(field_proxy_type), target  :: proxy_coord_output(3)
    type(function_space_type), pointer :: output_field_fs   => null()

    ! Variables for the gather to determine global domain sizes
    ! from the local partitioned ones
    integer(kind=i_def)              :: global_undf_checkpoint
    integer(kind=i_def), allocatable :: local_undf(:)
    integer(kind=i_def), allocatable :: all_undfs_checkpoint_domain(:)

    ! Factor to convert coords from radians to degrees if needed
    ! set as 1.0 for planar mesh
    real(kind=r_def) :: r2d

    if ( geometry == geometry_spherical ) then
     r2d = radians_to_degrees
    else
     r2d = 1.0_r_def
    end if

    ! Set k order value to 0 if unassigned
    if (present(k_order))then
      k_ord=k_order
    else
      k_ord=0
    end if

    ! Set up arrays to hold number of dofs for local and global domains
    allocate(local_undf(1))
    allocate(all_undfs_checkpoint_domain(get_comm_size()))

    all_undfs_checkpoint_domain = 0

    ! Create appropriate function space in order to be able to get the
    ! physical coordinates
    output_field_fs => function_space_collection%get_fs( mesh_id, &
                                                         k_ord, &
                                                         fs_id)

    ! Calculate the nodal coords for a field on the function space

    ! Set up fields to hold the output coordinates
    do i = 1,3
      call coord_output(i)%initialise( vector_space = output_field_fs )
    end do

    ! Convert field to physical nodal output & sample chi on nodal points
    call invoke_nodal_coordinates_kernel(coord_output, chi)

    ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
    if ( geometry == geometry_spherical ) then
      call invoke_pointwise_convert_xyz2llr(coord_output)
    end if

    ! Get proxies for coordinates so we can access them
    do i = 1,3
      proxy_coord_output(i) = coord_output(i)%get_proxy()
    end do

    ! Get the local value for undf
    local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

    !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!

    call all_gather ( local_undf, all_undfs_checkpoint_domain, 1 )

    ! Now get the global sum of undf across all ranks to set the global domain sizes
    ! for checkpoint domain
    global_undf_checkpoint = sum(all_undfs_checkpoint_domain)

    ! Calculate ibegin for each rank as we have the array of undfs in order
    ! we can just sum to get it.
    if (get_comm_rank() == 0) then
      ibegin_checkpoint = 0
    else
      ibegin_checkpoint = sum(all_undfs_checkpoint_domain(1:get_comm_rank()))
    end if

    ! Allocate coordinate arrays to be the size required for checkpoint domain.
    ! Essentially up to last owned dof of the current partition.
    allocate( checkpoint_lon( size( proxy_coord_output(1)%data(1: local_undf(1)))) )
    allocate( checkpoint_lat( size( proxy_coord_output(2)%data(1: local_undf(1)))) )

    ! Populate the arrays with data
    checkpoint_lon =  proxy_coord_output(1)%data(1: local_undf(1)) * r2d
    checkpoint_lat =  proxy_coord_output(2)%data(1: local_undf(1)) * r2d

    allocate(bnd_checkpoint_lon(1,size(checkpoint_lon)))
    allocate(bnd_checkpoint_lat(1,size(checkpoint_lat)))

    ! Construct bounds arrays
    bnd_checkpoint_lon=(reshape(checkpoint_lon, (/1, size(checkpoint_lon)/) ) )
    bnd_checkpoint_lat=(reshape(checkpoint_lat, (/1, size(checkpoint_lat)/) ) )

    ! Give coordinate information to the XIOS domain
    call xios_set_domain_attr(trim(domain_name), ni_glo=global_undf_checkpoint,    &
                              ibegin=ibegin_checkpoint, ni=local_undf(1),          &
                              type='unstructured')
    call xios_set_domain_attr(trim(domain_name), lonvalue_1d=checkpoint_lon,       &
                              latvalue_1d=checkpoint_lat)
    call xios_set_domain_attr(trim(domain_name), bounds_lon_1d=bnd_checkpoint_lon, &
                              bounds_lat_1d=bnd_checkpoint_lat)

    ! If we have requested to use domain index then get it and use it
    if (use_index) then

      ! Allocate domain_index - it is of size ndof_glob
      allocate(domain_index(output_field_fs%get_ndof_glob()))

      ! Populate domain_index for this rank
      call output_field_fs%get_global_dof_id(domain_index)

      ! temporary fix for higher-order domain decomposition
      if (k_ord > 0) then
        domain_index = domain_index/2
      end if

      ! Pass local portion of domain_index (up to undf)
      call xios_set_domain_attr(domain_name, i_index=int(domain_index(1:local_undf(1))))

    end if

    if ( allocated(checkpoint_lon) )     deallocate(checkpoint_lon)
    if ( allocated(checkpoint_lat) )     deallocate(checkpoint_lat)
    if ( allocated(domain_index) )    deallocate(domain_index)
    if ( allocated(bnd_checkpoint_lon) ) deallocate(bnd_checkpoint_lon)
    if ( allocated(bnd_checkpoint_lat) ) deallocate(bnd_checkpoint_lat)
    if ( allocated(local_undf) )      deallocate(local_undf)
    if ( allocated(all_undfs_checkpoint_domain) ) deallocate(all_undfs_checkpoint_domain)

    nullify( output_field_fs )
    return
  end subroutine checkpoint_domain_init

  !> @brief   Initialises unstructured XIOS domain from function space
  !> @details Calculates coordinates from function space and local mesh, obtains
  !>          local portion of the domain index and passes information to XIOS
  !>          domain object
  !>
  !> @param[in]  domain_id   The name of the XIOS domain
  !> @param[in]  mesh_id     The id of the partitioned mesh
  !> @param[in]  fs_id       The id of the function space corresponding to the domain
  !> @param[in]  chi         Input coordinate field
  !> @param[in]  lon_bounds  Array of longitude coords making up the domain bounds
  !> @param[in]  lat_bounds  Array of latitude coords making up the domain bounds
  !>
  subroutine init_xios_ugrid_domain( domain_id, mesh_id, fs_id, chi, lon_bounds, lat_bounds )

    implicit none

    character(len=*),       intent(in) :: domain_id
    integer(kind=i_def),    intent(in) :: mesh_id
    integer(kind=i_native), intent(in) :: fs_id
    type(field_type),       intent(in) :: chi(:)

    type(function_space_type), pointer :: domain_fs   => null()

    type( field_type )              :: coord_output(3)
    type(field_proxy_type), target  :: proxy_coord_output(3)

    ! Variables for the gather to determine global domain sizes
    ! from the local partitioned ones
    integer(kind=i_def)              :: global_undf, n_levels, i, ibegin, local_domain_size
    integer(kind=i_def), allocatable :: local_undf(:), all_undfs(:)
    real(dp_xios),       allocatable :: dp_levels(:)
    real(dp_xios),       allocatable :: lat_data(:)
    real(dp_xios),       allocatable :: lon_data(:)
    real(dp_xios),       allocatable :: lat_bounds(:,:)
    real(dp_xios),       allocatable :: lon_bounds(:,:)
    integer(kind=i_def), allocatable :: domain_index(:)

    ! Factor to convert coords from radians to degrees if needed
    ! set as 1.0 for planar mesh
    real(kind=r_def) :: r2d

    if ( geometry == geometry_spherical ) then
      r2d = radians_to_degrees
    else
      r2d = 1.0_r_def
    end if

    ! Set up arrays for AllGather
    allocate(local_undf(1))
    allocate(all_undfs(get_comm_size()))

    all_undfs = 0

    ! Here we use information from the input function space to calculate the
    ! physical coordinates for the horizontal domain
    domain_fs => function_space_collection%get_fs( mesh_id, element_order, fs_id )

    ! Get the function space levels information
    dp_levels = real( domain_fs%get_levels(), kind=dp_xios )
    n_levels = size(dp_levels)

    ! Set up fields to hold the output coordinates
    do i = 1,3
      call coord_output(i)%initialise( vector_space = domain_fs )
    end do

    !=========================================================================
    ! These calls need to be replaced with a working infrastructure alternative:
    ! Convert field to physical nodal output & sample chi on nodal points
    call invoke_nodal_coordinates_kernel(coord_output, chi)

    ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
    if ( geometry == geometry_spherical ) then
      call invoke_pointwise_convert_xyz2llr(coord_output)
    end if
    !=========================================================================

    ! Get proxies for coordinates so we can access them
    do i = 1,3
      proxy_coord_output(i) = coord_output(i)%get_proxy()
    end do

    ! Get the local value for undf
    local_undf(1)  = domain_fs%get_last_dof_owned()
    local_domain_size = local_undf(1)/n_levels

    allocate( lat_data(local_domain_size) )
    allocate( lon_data(local_domain_size) )

    lat_data = proxy_coord_output(2)%data(1: local_undf(1):n_levels) * r2d
    lon_data = proxy_coord_output(1)%data(1: local_undf(1):n_levels) * r2d

    !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call all_gather ( local_undf, all_undfs, 1 )

    ! Adjust size of data taking into account how many levels we have (same for each
    ! partition as we only partition horizontally)
    all_undfs = all_undfs/n_levels

    ! Now get the global sum of undf across all ranks to set the global domain sizes
    ! for xios node domain
    global_undf = sum(all_undfs)

    ! Calculate ibegin for each rank (as we have the array of undfs in order
    ! we can just sum to get it)
    if (get_comm_rank() == 0) then
      ibegin = 0
    else
      ibegin = sum(all_undfs(1:get_comm_rank()))
    end if

    ! Populate domain_index for this rank for relevant function space
    allocate( domain_index(local_domain_size) )
    if (fs_id == W0) then
      call domain_fs%get_global_vert_dof_id_2d(domain_index)
    else if (fs_id == W2H) then
      call domain_fs%get_global_edge_dof_id_2d(domain_index)
    else if (fs_id == W3) then
      call domain_fs%get_global_cell_dof_id_2d(domain_index)
    end if

    ! Pass domain attibutes to XIOS
    call xios_set_domain_attr( trim(domain_id), ni_glo=global_undf, &
                               ibegin=ibegin, &
                               ni=local_domain_size, &
                               type='unstructured' )
    call xios_set_domain_attr( trim(domain_id), lonvalue_1d=lon_data, &
                               latvalue_1d=lat_data )
    call xios_set_domain_attr( trim(domain_id), bounds_lon_1d=lon_bounds, &
                               bounds_lat_1d=lat_bounds )

    ! Pass local portion of domain_index to XIOS
    call xios_set_domain_attr( domain_id, i_index=int( domain_index( 1 : local_domain_size ) ) )

    ! Tidy up time
    deallocate( local_undf, all_undfs )
    deallocate( dp_levels )
    deallocate( lat_data )
    deallocate( lon_data )
    deallocate( domain_index )
    nullify(domain_fs)

  end subroutine init_xios_ugrid_domain

  !> @brief   Initialises XIOS axis from function space and mesh
  !> @details Calculates vertical levels and coordinates from function space
  !>          and local mesh and passes information to XIOS axis object
  !>
  !> @param[in]  axis_id  The name of the XIOS axis
  !> @param[in]  mesh_id  The id of the partitioned mesh
  !> @param[in]  fs_id    The id of the function space corresponding to the
  !>                      domain
  !>
  subroutine init_xios_axis( axis_id, mesh_id, fs_id )

    implicit none

    character(len=*),       intent(in) :: axis_id
    integer(kind=i_def),    intent(in) :: mesh_id
    integer(kind=i_native), intent(in) :: fs_id

    type(function_space_type), pointer :: domain_fs => null()
    real(dp_xios), allocatable         :: dp_levels(:)
    integer(kind=i_def)                :: n_levels

    domain_fs => function_space_collection%get_fs( mesh_id,       &
                                                   element_order, &
                                                   fs_id )

    ! Get the function space levels information
    dp_levels = real( domain_fs%get_levels(), kind=dp_xios )
    n_levels = size(dp_levels)

    call xios_set_axis_attr( trim(axis_id),  &
                             n_glo=n_levels, &
                             value=dp_levels )

  end subroutine init_xios_axis

end module lfric_xios_io_mod
