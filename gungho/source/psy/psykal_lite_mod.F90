!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the PSy layer

!> @details Contains hand-rolled versions of the PSy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_mod

  use field_mod,                    only : field_type, field_proxy_type
  use scalar_mod,                   only : scalar_type
  use operator_mod,                 only : operator_type, operator_proxy_type
  use constants_mod,                only : r_def, i_def, cache_block
  use mesh_mod,                     only : mesh_type
  use function_space_mod,           only : BASIS, DIFF_BASIS

  use quadrature_xyoz_mod, only : quadrature_xyoz_type, &
                                  quadrature_xyoz_proxy_type
  use quadrature_face_mod, only : quadrature_face_type, &
                                  quadrature_face_proxy_type

  implicit none
  public

contains

!> Non pointwise Kernels

  !-------------------------------------------------------------------------------
  !> io_mod uses this routine. However, because io_mod is not a algorithm its currently
  !> not clear if it should call into PSy. #1253 will address this point and remove
  !> once decided.
  subroutine invoke_nodal_coordinates_kernel(nodal_coords, chi )
    use nodal_coordinates_kernel_mod, only: nodal_coordinates_code
    use mesh_mod,                     only: mesh_type
    implicit none

    type(field_type), intent(inout)      :: nodal_coords(3)
    type(field_type), intent(in)         :: chi(3)

    type(field_proxy_type) :: x_p(3), chi_p(3)

    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x
    integer                 :: undf_chi, undf_x
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf( )
    undf_x = x_p(1)%vspace%get_undf()
    nodes_x => x_p(1)%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    mesh => x_p(1)%vspace%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call nodal_coordinates_code(nlayers, &
                                   x_p(1)%data, &
                                   x_p(2)%data, &
                                   x_p(3)%data, &
                                   chi_p(1)%data, &
                                   chi_p(2)%data, &
                                   chi_p(3)%data, &
                                   ndf_x, undf_x, map_x, &
                                   ndf_chi, undf_chi, map_chi, &
                                   basis_chi &
                                  )
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_nodal_coordinates_kernel

  !-------------------------------------------------------------------------------
  !> io_mod uses this routine. However, because io_mod is not a algorithm its currently
  !> not clear if it should call into PSy. #1253 will address this point and remove
  !> once decided.
  subroutine invoke_nodal_xyz_coordinates_kernel(nodal_coords, chi, panel_id)
    use nodal_xyz_coordinates_kernel_mod, only: nodal_xyz_coordinates_code
    use mesh_mod,                     only: mesh_type
    implicit none

    type(field_type), intent(inout)      :: nodal_coords(3)
    type(field_type), intent(in)         :: chi(3)
    type(field_type), intent(in)         :: panel_id

    type(field_proxy_type) :: x_p(3), chi_p(3), panel_id_proxy

    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x, ndf_pid
    integer                 :: undf_chi, undf_x, undf_pid
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_pid(:) => null()
    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    panel_id_proxy = panel_id%get_proxy()

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf()
    undf_x = x_p(1)%vspace%get_undf()
    nodes_x => x_p(1)%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf()
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space()

    ndf_pid  = panel_id_proxy%vspace%get_ndf()
    undf_pid = panel_id_proxy%vspace%get_undf()

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    if (panel_id_proxy%is_dirty(depth=1)) then
      call panel_id_proxy%halo_exchange(depth=1)
    end if

    mesh => x_p(1)%vspace%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       map_pid => panel_id_proxy%vspace%get_cell_dofmap( cell )
       call nodal_xyz_coordinates_code(nlayers, &
                                       x_p(1)%data, &
                                       x_p(2)%data, &
                                       x_p(3)%data, &
                                       chi_p(1)%data, &
                                       chi_p(2)%data, &
                                       chi_p(3)%data, &
                                       panel_id_proxy%data, &
                                       ndf_x, undf_x, map_x, &
                                       ndf_chi, undf_chi, map_chi, &
                                       basis_chi, &
                                       ndf_pid, undf_pid, map_pid &
                                      )
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_nodal_xyz_coordinates_kernel

!-------------------------------------------------------------------------------
  subroutine invoke_compute_dof_level_kernel(level)

  use compute_dof_level_kernel_mod, only: compute_dof_level_code
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(inout) :: level
  type(field_proxy_type) :: l_p
  integer :: cell, ndf, undf
  real(kind=r_def), pointer :: nodes(:,:) => null()
  integer, pointer :: map(:) => null()
  type(mesh_type), pointer :: mesh => null()
  l_p = level%get_proxy()
  undf = l_p%vspace%get_undf()
  ndf  = l_p%vspace%get_ndf()
  nodes => l_p%vspace%get_nodes( )

  mesh => l_p%vspace%get_mesh()
  do cell = 1,mesh%get_last_halo_cell(1)
    map => l_p%vspace%get_cell_dofmap(cell)
    call compute_dof_level_code(l_p%vspace%get_nlayers(),                 &
                                l_p%data,                                 &
                                ndf,                                      &
                                undf,                                     &
                                map,                                      &
                                nodes                                     &
                               )
  end do
  call l_p%set_dirty()

  end subroutine invoke_compute_dof_level_kernel

!-------------------------------------------------------------------------------
!> invoke_subgrid_coeffs: Invoke the calculation of subgrid rho coefficients
subroutine invoke_subgrid_coeffs(a0,a1,a2,rho,cell_orientation,direction,rho_approximation_stencil_extent,halo_depth_to_compute)

    use flux_direction_mod,        only: x_direction, y_direction
    use stencil_dofmap_mod,        only: stencil_dofmap_type, &
                                         STENCIL_1DX,         &
                                         STENCIL_1DY
    use subgrid_coeffs_kernel_mod, only: subgrid_coeffs_code
    use subgrid_config_mod,        only: rho_approximation
    use mesh_mod,                  only: mesh_type
    use log_mod,                   only: log_event, LOG_LEVEL_ERROR

    implicit none

    type( field_type ), intent( inout ) :: a0
    type( field_type ), intent( inout ) :: a1
    type( field_type ), intent( inout ) :: a2
    type( field_type ), intent( in )    :: rho
    type( field_type ), intent( in )    :: cell_orientation
    integer, intent(in)                 :: direction
    integer, intent(in)                 :: rho_approximation_stencil_extent
    integer, intent(in)                 :: halo_depth_to_compute

    type( field_proxy_type )            :: rho_proxy
    type( field_proxy_type )            :: a0_proxy
    type( field_proxy_type )            :: a1_proxy
    type( field_proxy_type )            :: a2_proxy
    type( field_proxy_type )            :: cell_orientation_proxy

    type(stencil_dofmap_type), pointer  :: map_x_w3 => null()
    type(stencil_dofmap_type), pointer  :: map_y_w3 => null()
    integer, pointer                    :: map_w3(:) => null()
    integer, pointer                    :: stencil_map(:,:) => null()
    integer                             :: rho_stencil_size
    integer                 :: cell
    integer                 :: nlayers
    integer                 :: ndf_w3
    integer                 :: undf_w3
    type(mesh_type), pointer :: mesh => null()
    integer                  :: d
    logical                  :: swap
    integer                  :: ncells_to_iterate

    a0_proxy   = a0%get_proxy()
    a1_proxy   = a1%get_proxy()
    a2_proxy   = a2%get_proxy()
    rho_proxy  = rho%get_proxy()
    cell_orientation_proxy = cell_orientation%get_proxy()

    undf_w3 = rho_proxy%vspace%get_undf()
    ndf_w3  = rho_proxy%vspace%get_ndf()
    nlayers = rho_proxy%vspace%get_nlayers()

    ! Note stencil grid types are of the form:
    !                                   |5|
    !                                   |3|
    ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
    !                                   |2|
    !                                   |4|
    map_x_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_approximation_stencil_extent)
    map_y_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_approximation_stencil_extent)

    rho_stencil_size = map_x_w3%get_size()

    mesh => a0_proxy%vspace%get_mesh()

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_proxy%halo_exchange(depth=mesh%get_halo_depth())

    if (halo_depth_to_compute==0) then
      ncells_to_iterate = mesh%get_last_edge_cell()
    elseif (halo_depth_to_compute > 0) then
      ncells_to_iterate = mesh%get_last_halo_cell(halo_depth_to_compute)
    else
      call log_event( "Error: negative halo_depth_to_compute value in subgrid coeffs call", LOG_LEVEL_ERROR )
    endif

    !NOTE: The default looping limits for this type of field would be
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and some halo cells.
    do cell = 1, ncells_to_iterate

      map_w3 => rho_proxy%vspace%get_cell_dofmap(cell)

      if (direction == x_direction) then
        if (nint(cell_orientation_proxy%data(map_w3(1))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1))) == 4) then
          stencil_map => map_y_w3%get_dofmap(cell)
        else
          stencil_map => map_x_w3%get_dofmap(cell)
        end if
      elseif (direction == y_direction) then
        if (nint(cell_orientation_proxy%data(map_w3(1))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1))) == 4) then
          stencil_map => map_x_w3%get_dofmap(cell)
        else
          stencil_map => map_y_w3%get_dofmap(cell)
        end if
      end if

      call subgrid_coeffs_code( nlayers,                                  &
                                rho_approximation,                        &
                                undf_w3,                                  &
                                rho_proxy%data,                           &
                                cell_orientation_proxy%data,              &
                                ndf_w3,                                   &
                                rho_stencil_size,                         &
                                stencil_map,                              &
                                direction,                                &
                                a0_proxy%data,                            &
                                a1_proxy%data,                            &
                                a2_proxy%data                             &
                                )

    end do
    call a0_proxy%set_dirty()
    call a1_proxy%set_dirty()
    call a2_proxy%set_dirty()

  end subroutine invoke_subgrid_coeffs


!-------------------------------------------------------------------------------
!> invoke_subgrid_coeffs: Invoke the calculation of subgrid rho coefficients
!>                        The routine also includes a special type of halo
!>                        exchange where the values in the halos need to be
!>                        corrected. This is due to 1D direction updates of the
!>                        density field and the panels of the cubed-sphere having
!>                        different orientation.
  subroutine invoke_subgrid_coeffs_conservative( a0_x,                             &
                                                 a1_x,                             &
                                                 a2_x,                             &
                                                 a0_y,                             &
                                                 a1_y,                             &
                                                 a2_y,                             &
                                                 rho_x,                            &
                                                 rho_y,                            &
                                                 rho_x_halos_corrected,            &
                                                 rho_y_halos_corrected,            &
                                                 cell_orientation,                 &
                                                 rho_approximation_stencil_extent, &
                                                 dep_pt_stencil_extent  )

    use flux_direction_mod,               only: x_direction, y_direction
    use stencil_dofmap_mod,               only: stencil_dofmap_type, &
                                                STENCIL_1DX,         &
                                                STENCIL_1DY
    use subgrid_coeffs_kernel_mod,        only: subgrid_coeffs_code
    use subgrid_config_mod,               only: rho_approximation
    use mesh_mod,                         only: mesh_type
    use log_mod,                          only: log_event, LOG_LEVEL_ERROR
    use cosmic_halo_correct_x_kernel_mod, only: cosmic_halo_correct_x_code
    use cosmic_halo_correct_y_kernel_mod, only: cosmic_halo_correct_y_code

    implicit none

    type( field_type ), intent( inout ) :: a0_x
    type( field_type ), intent( inout ) :: a1_x
    type( field_type ), intent( inout ) :: a2_x
    type( field_type ), intent( inout ) :: a0_y
    type( field_type ), intent( inout ) :: a1_y
    type( field_type ), intent( inout ) :: a2_y
    type( field_type ), intent( in )    :: rho_x
    type( field_type ), intent( in )    :: rho_y
    type( field_type ), intent( inout ) :: rho_x_halos_corrected
    type( field_type ), intent( inout ) :: rho_y_halos_corrected
    type( field_type ), intent( in )    :: cell_orientation
    integer, intent(in)                 :: rho_approximation_stencil_extent
    integer, intent(in)                 :: dep_pt_stencil_extent

    type( field_proxy_type )            :: rho_x_proxy
    type( field_proxy_type )            :: rho_y_proxy
    type( field_proxy_type )            :: rho_x_halos_corrected_proxy
    type( field_proxy_type )            :: rho_y_halos_corrected_proxy
    type( field_proxy_type )            :: a0_x_proxy
    type( field_proxy_type )            :: a1_x_proxy
    type( field_proxy_type )            :: a2_x_proxy
    type( field_proxy_type )            :: a0_y_proxy
    type( field_proxy_type )            :: a1_y_proxy
    type( field_proxy_type )            :: a2_y_proxy
    type( field_proxy_type )            :: cell_orientation_proxy

    type(stencil_dofmap_type), pointer  :: map_x_w3 => null()
    type(stencil_dofmap_type), pointer  :: map_y_w3 => null()
    type(mesh_type), pointer            :: mesh => null()

    integer, pointer                    :: map_w3(:,:) => null()
    integer, pointer                    :: stencil_map(:,:) => null()
    integer                             :: rho_stencil_size
    integer                             :: cell
    integer                             :: nlayers
    integer                             :: ndf_w3
    integer                             :: undf_w3
    integer                             :: d
    logical                             :: swap
    integer                             :: ncells_to_iterate


    a0_x_proxy   = a0_x%get_proxy()
    a1_x_proxy   = a1_x%get_proxy()
    a2_x_proxy   = a2_x%get_proxy()
    a0_y_proxy   = a0_y%get_proxy()
    a1_y_proxy   = a1_y%get_proxy()
    a2_y_proxy   = a2_y%get_proxy()

    rho_x_proxy  = rho_x%get_proxy()
    rho_y_proxy  = rho_y%get_proxy()
    rho_x_halos_corrected_proxy = rho_x_halos_corrected%get_proxy()
    rho_y_halos_corrected_proxy = rho_y_halos_corrected%get_proxy()

    cell_orientation_proxy = cell_orientation%get_proxy()

    undf_w3 = rho_x_proxy%vspace%get_undf()
    ndf_w3  = rho_x_proxy%vspace%get_ndf()
    nlayers = rho_x_proxy%vspace%get_nlayers()

    mesh => rho_x_proxy%vspace%get_mesh()

    map_w3 => rho_x_proxy%vspace%get_whole_dofmap()

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_x_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_x_proxy%halo_exchange(depth=mesh%get_halo_depth())

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_y_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_y_proxy%halo_exchange(depth=mesh%get_halo_depth())

    ! Loop over all core and halo cells.
    do cell=1,mesh%get_ncells_2d()

      call cosmic_halo_correct_x_code(  nlayers,                            &
                                        rho_x_halos_corrected_proxy%data,   &
                                        rho_x_proxy%data,                   &
                                        rho_y_proxy%data,                   &
                                        cell_orientation_proxy%data,        &
                                        ndf_w3,                             &
                                        undf_w3,                            &
                                        map_w3(:,cell))
    end do

    ! Loop over all core and halo cells.
    do cell=1,mesh%get_ncells_2d()

      call cosmic_halo_correct_y_code(  nlayers,                            &
                                        rho_y_halos_corrected_proxy%data,   &
                                        rho_x_proxy%data,                   &
                                        rho_y_proxy%data,                   &
                                        cell_orientation_proxy%data,        &
                                        ndf_w3,                             &
                                        undf_w3,                            &
                                        map_w3(:,cell))
    end do

    map_x_w3 => rho_x_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_approximation_stencil_extent)
    map_y_w3 => rho_y_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_approximation_stencil_extent)

    rho_stencil_size = map_x_w3%get_size()

    if (dep_pt_stencil_extent==0) then
      ncells_to_iterate = mesh%get_last_edge_cell()
    elseif (dep_pt_stencil_extent > 0) then
      ncells_to_iterate = mesh%get_last_halo_cell(dep_pt_stencil_extent)
    else
      call log_event( "Error: negative dep_pt_stencil_extent value in subgrid coeffs call", LOG_LEVEL_ERROR )
    endif

    ! Note stencil grid types are of the form:
    !                                   |5|
    !                                   |3|
    ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
    !                                   |2|
    !                                   |4|

    ! NOTE: The default looping limits for this type of field would be
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and some halo cells.

    ! Calculate a0_y, a1_y, a2_y associated with rho_y but in the x_direction
    do cell = 1, ncells_to_iterate
      if (nint(cell_orientation_proxy%data(map_w3(1,cell))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1,cell))) == 4) then
        stencil_map => map_y_w3%get_dofmap(cell)
      else
        stencil_map => map_x_w3%get_dofmap(cell)
      end if
      call subgrid_coeffs_code( nlayers,                                    &
                                rho_approximation,                          &
                                undf_w3,                                    &
                                rho_y_halos_corrected_proxy%data,           &
                                cell_orientation_proxy%data,                &
                                ndf_w3,                                     &
                                rho_stencil_size,                           &
                                stencil_map,                                &
                                x_direction,                                &
                                a0_y_proxy%data,                            &
                                a1_y_proxy%data,                            &
                                a2_y_proxy%data                             &
                                )
    end do


    ! Calculate a0_x, a1_x, a2_x associated with rho_x but in the y_direction
    do cell = 1, ncells_to_iterate
      if (nint(cell_orientation_proxy%data(map_w3(1,cell))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1,cell))) == 4) then
        stencil_map => map_x_w3%get_dofmap(cell)
      else
        stencil_map => map_y_w3%get_dofmap(cell)
      end if
      call subgrid_coeffs_code( nlayers,                                  &
                                rho_approximation,                        &
                                undf_w3,                                  &
                                rho_x_halos_corrected_proxy%data,         &
                                cell_orientation_proxy%data,              &
                                ndf_w3,                                   &
                                rho_stencil_size,                         &
                                stencil_map,                              &
                                y_direction,                              &
                                a0_x_proxy%data,                          &
                                a1_x_proxy%data,                          &
                                a2_x_proxy%data                           &
                                )
    end do
    call a0_x_proxy%set_dirty()
    call a1_x_proxy%set_dirty()
    call a2_x_proxy%set_dirty()
    call a0_y_proxy%set_dirty()
    call a1_y_proxy%set_dirty()
    call a2_y_proxy%set_dirty()

  end subroutine invoke_subgrid_coeffs_conservative


!------------------------------------------------------------------------------
! One of the reasons (but not the only one) for this "light" implementation is
! passing the double precision deltaT value to fv_mass_flux_code. This should
! not be taken as a requirement, it is simply expedient to get the clock change
! on trunk. It is probably the wrong thing to be doing in the long run.
!
subroutine invoke_fv_mass_fluxes( rho,            &
                                  dep_pts,        &
                                  mass_flux,      &
                                  a0_coeffs,      &
                                  a1_coeffs,      &
                                  a2_coeffs,      &
                                  direction,      &
                                  stencil_extent )

  use fv_mass_flux_kernel_mod,      only: fv_mass_flux_code
  use flux_direction_mod,           only: x_direction, y_direction
  use stencil_dofmap_mod,           only: stencil_dofmap_type, &
                                          STENCIL_1DX, STENCIL_1DY
  use timestepping_config_mod,      only: dt
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(in)      :: rho
  type(field_type), intent(in)      :: dep_pts
  type(field_type), intent(inout)   :: mass_flux
  type(field_type), intent(in)      :: a0_coeffs
  type(field_type), intent(in)      :: a1_coeffs
  type(field_type), intent(in)      :: a2_coeffs
  integer, intent(in)               :: direction
  integer, intent(in)               :: stencil_extent

  type( field_proxy_type )  :: mass_flux_proxy, dep_pts_proxy, rho_proxy
  type( field_proxy_type )  :: a0_coeffs_proxy, a1_coeffs_proxy, a2_coeffs_proxy

  type(stencil_dofmap_type), pointer  :: map => null()

  integer, pointer :: map_rho(:) => null()
  integer, pointer :: map_w2(:) => null()
  integer, pointer :: stencil_map(:,:) => null()
  integer          :: stencil_size

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: cell
  integer :: nlayers
  type(mesh_type), pointer :: mesh => null()
  integer                  :: d
  logical                  :: swap

  rho_proxy     = rho%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()

  ndf_w3  = rho_proxy%vspace%get_ndf()
  undf_w3 = rho_proxy%vspace%get_undf()

  ndf_w2  = dep_pts_proxy%vspace%get_ndf()
  undf_w2 = dep_pts_proxy%vspace%get_undf()

  a0_coeffs_proxy = a0_coeffs%get_proxy()
  a1_coeffs_proxy = a1_coeffs%get_proxy()
  a2_coeffs_proxy = a2_coeffs%get_proxy()
  mass_flux_proxy = mass_flux%get_proxy()

  nlayers = rho_proxy%vspace%get_nlayers()

  ! Note stencil grid types are of the form:
  !                                   |5|
  !                                   |3|
  ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
  !                                   |2|
  !                                   |4|
  if (direction == x_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,stencil_extent)
  elseif (direction == y_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,stencil_extent)
  end if
  stencil_size = map%get_size()

  swap = .false.
  do d = 1,stencil_extent
    if (rho_proxy%is_dirty(depth=d)) swap = .true.
  end do
  if ( swap ) call rho_proxy%halo_exchange(depth=stencil_extent)

  mesh => rho_proxy%vspace%get_mesh()
  ! NOTE: The default looping limits for this type of field would be
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! in order to function correctly. See ticket #1058.
  ! The kernel loops over all core cells only.
  do cell = 1, mesh%get_last_edge_cell()
      map_rho => rho_proxy%vspace%get_cell_dofmap( cell )
      map_w2 => dep_pts_proxy%vspace%get_cell_dofmap( cell )

      stencil_map => map%get_dofmap(cell)

      call fv_mass_flux_code(  nlayers,                     &
                               undf_w3,                     &
                               ndf_w3,                      &
                               map_rho,                     &
                               rho_proxy%data,              &
                               a0_coeffs_proxy%data,        &
                               a1_coeffs_proxy%data,        &
                               a2_coeffs_proxy%data,        &
                               undf_w2,                     &
                               ndf_w2,                      &
                               map_w2,                      &
                               mass_flux_proxy%data,        &
                               dep_pts_proxy%data,          &
                               stencil_size,                &
                               stencil_map,                 &
                               direction,                   &
                               dt )

  end do
  call a0_coeffs_proxy%set_dirty()
  call a1_coeffs_proxy%set_dirty()
  call a2_coeffs_proxy%set_dirty()

end subroutine invoke_fv_mass_fluxes


!-------------------------------------------------------------------------------
!> invoke_calc_deppts: Invoke the calculation of departure points in 1D
subroutine invoke_calc_deppts(  u_n,                  &
                                u_np1,                &
                                dep_pts,              &
                                cell_orientation,     &
                                direction,            &
                                dep_pt_method,        &
                                dep_pt_stencil_extent )

  use calc_departure_point_kernel_mod,  only : calc_departure_point_code
  use stencil_dofmap_mod,               only : stencil_dofmap_type, &
                                               STENCIL_1DX, &
                                               STENCIL_1DY
  use flux_direction_mod,               only : x_direction, y_direction
  use mesh_mod,                         only : mesh_type
  implicit none

  type( field_type ), intent( in )    :: u_n
  type( field_type ), intent( in )    :: u_np1
  type( field_type ), intent( inout ) :: dep_pts
  type( field_type ), intent( in )    :: cell_orientation
  integer, intent(in)                 :: direction
  integer, intent(in)                 :: dep_pt_method
  integer, intent(in)                 :: dep_pt_stencil_extent

  type( field_proxy_type )        :: u_n_proxy
  type( field_proxy_type )        :: u_np1_proxy
  type( field_proxy_type )        :: dep_pts_proxy
  type( field_proxy_type )        :: cell_orientation_proxy
  type(stencil_dofmap_type), pointer  :: map=>null()
  type(stencil_dofmap_type), pointer  :: map_w3=>null()

  integer, pointer        :: stencil_map_w2(:,:) => null()
  integer, pointer        :: stencil_map_w3(:,:) => null()
  integer                 :: transport_stencil_size

  integer                 :: cell
  integer                 :: nlayers
  integer                 :: ndf_w2
  integer                 :: undf_w2
  integer                 :: ndf_w3
  integer                 :: undf_w3
  type(mesh_type), pointer :: mesh => null()

  u_n_proxy    = u_n%get_proxy()
  u_np1_proxy  = u_np1%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()
  cell_orientation_proxy = cell_orientation%get_proxy()

  ndf_w2  = u_n_proxy%vspace%get_ndf()
  undf_w2 = u_n_proxy%vspace%get_undf()

  ndf_w3  =   cell_orientation_proxy%vspace%get_ndf()
  undf_w3 =   cell_orientation_proxy%vspace%get_undf()

  if (direction == x_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,dep_pt_stencil_extent)
    map_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,dep_pt_stencil_extent)
  elseif (direction == y_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,dep_pt_stencil_extent)
    map_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,dep_pt_stencil_extent)
  endif
  transport_stencil_size = map%get_size()

  nlayers = u_n_proxy%vspace%get_nlayers()

  mesh => u_n_proxy%vspace%get_mesh()
  !NOTE: The default looping limits for this type of field would be
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! in order to function correctly. See ticket #1058.
  ! The kernel loops over all core cells only.
  do cell=1,mesh%get_last_edge_cell()

    stencil_map_w2 => map%get_dofmap(cell)
    stencil_map_w3 => map_w3%get_dofmap(cell)

    call calc_departure_point_code( nlayers,                      &
                                    dep_pts_proxy%data,           &
                                    transport_stencil_size,       &
                                    undf_w2,                      &
                                    ndf_w2,                       &
                                    stencil_map_w2,               &
                                    undf_w3,                      &
                                    ndf_w3,                       &
                                    stencil_map_w3,               &
                                    cell_orientation_proxy%data,  &
                                    u_n_proxy%data,               &
                                    u_np1_proxy%data,             &
                                    direction,                    &
                                    dep_pt_method )

  end do
  call dep_pts_proxy%set_dirty()

end subroutine invoke_calc_deppts

!-------------------------------------------------------------------------------
!> invoke_multiply_field_data:  z =  x * y
  subroutine invoke_multiply_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy,     &
                                          field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("PSy:multiply_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("PSy:multiply_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) * field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res_proxy%vspace%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_multiply_field_data

!-------------------------------------------------------------------------------
! Implemented in #965, kernel requires stencil support. Note that the w2_field is
! required to obtain the W2 stencil_cross which is used in determining
! orientation of cells in the halo
  subroutine invoke_mpi_calc_cell_orientation(w2_field,cell_orientation)

    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod,               only : stencil_dofmap_type, &
                                                 STENCIL_CROSS
    use calc_cell_orientation_kernel_mod, only : calc_cell_orientation_code

    implicit none

    type(field_type), intent(in)      :: w2_field
    type(field_type), intent(inout)   :: cell_orientation

    integer                 :: cell, nlayers
    integer                 :: ndf_w3
    integer                 :: undf_w3
    integer                 :: ndf_w2
    integer, pointer        :: map_w3(:) => null()

    type(field_proxy_type) :: cell_orientation_proxy
    type(field_proxy_type) :: w2_field_proxy

    type(mesh_type), pointer           :: mesh => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_w2 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size


    cell_orientation_proxy = cell_orientation%get_proxy()
    w2_field_proxy = w2_field%get_proxy()
    mesh => w2_field_proxy%vspace%get_mesh()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = w2_field_proxy%vspace%get_ndf( )

    ! Obtain the stencil for core cells only
    cross_stencil_w2 => w2_field_proxy%vspace%get_stencil_dofmap(             &
                                        STENCIL_CROSS, mesh%get_halo_depth())
    cross_stencil_w2_map => cross_stencil_w2%get_whole_dofmap()

    cross_stencil_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(     &
                                        STENCIL_CROSS, mesh%get_halo_depth())
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only
      map_w3 => cell_orientation_proxy%vspace%get_cell_dofmap(cell)

      call calc_cell_orientation_code(  nlayers,                              &
                                        cell_orientation_proxy%data,          &
                                        undf_w3,                              &
                                        ndf_w3,                               &
                                        map_w3,                               &
                                        ndf_w2,                               &
                                        cross_stencil_w3_size,                &
                                        cross_stencil_w2_map(:,:,cell),       &
                                        cross_stencil_w3_map(:,:,cell) )
    end do

  end subroutine invoke_mpi_calc_cell_orientation

  !-------------------------------------------------------------------------------
  !> This routine is called from psykal_lite due to the variable cell_orientation
  !> being passed into the kernel.
  !> The cell_orientation field should not be halo exchanged across panels as the
  !> orientation of cells is local to its own panel on the cubed-sphere.
  subroutine invoke_fv_divergence( mass_divergence,          &
                                   mass_flux_x,              &
                                   mass_flux_y,              &
                                   cell_orientation,         &
                                   direction )

    use mesh_mod,                         only : mesh_type
    use fv_divergence_kernel_mod,         only : fv_divergence_code
    use flux_direction_mod,               only : x_direction, y_direction

    implicit none

    type(field_type), intent(inout) :: mass_divergence
    type(field_type), intent(in)    :: mass_flux_x
    type(field_type), intent(in)    :: mass_flux_y
    type(field_type), intent(in)    :: cell_orientation
    integer, intent(in)             :: direction

    type(mesh_type), pointer        :: mesh => null()

    integer, pointer                :: map_w3(:,:) => null()
    integer, pointer                :: map_w2(:,:) => null()

    type(field_proxy_type)          :: cell_orientation_proxy
    type(field_proxy_type)          :: mass_flux_x_proxy, mass_flux_y_proxy
    type(field_proxy_type)          :: mass_divergence_proxy

    integer                         :: cell, nlayers
    integer                         :: ndf_w3, undf_w3
    integer                         :: ndf_w2, undf_w2

    cell_orientation_proxy           = cell_orientation%get_proxy()
    mass_divergence_proxy            = mass_divergence%get_proxy()
    mass_flux_x_proxy                = mass_flux_x%get_proxy()
    mass_flux_y_proxy                = mass_flux_y%get_proxy()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = mass_flux_x_proxy%vspace%get_ndf( )
    undf_w2 = mass_flux_x_proxy%vspace%get_undf()

    mesh   => mass_flux_x_proxy%vspace%get_mesh()
    map_w3 => cell_orientation_proxy%vspace%get_whole_dofmap()
    map_w2 => mass_flux_x_proxy%vspace%get_whole_dofmap()

    ! There is no automatic halo exchange on purpose at the moment. Since if there
    ! was then the x and y directional components would not be respected due to
    ! different panel orientations on the cubed-sphere.
    ! A ticket, #1147, has been created which addresses this issue of dealing
    ! with panel orientation when halo exchanging W2 fields.
    ! A similar implementation was made for invoke_subgrid_coeffs_conservative
    ! which dealt with W3 fields and ticket #1087 implemented this.
    if (direction == x_direction) then

      do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only

        call  fv_divergence_code(  nlayers,                             &
                                   mass_divergence_proxy%data,          &
                                   cell_orientation_proxy%data,         &
                                   mass_flux_x_proxy%data,              &
                                   direction,                           &
                                   ndf_w3,                              &
                                   undf_w3,                             &
                                   map_w3(:,cell),                      &
                                   ndf_w2,                              &
                                   undf_w2,                             &
                                   map_w2(:,cell) )

      end do

    elseif (direction == y_direction) then

      do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only

        call  fv_divergence_code(  nlayers,                             &
                                   mass_divergence_proxy%data,          &
                                   cell_orientation_proxy%data,         &
                                   mass_flux_y_proxy%data,              &
                                   direction,                           &
                                   ndf_w3,                              &
                                   undf_w3,                             &
                                   map_w3(:,cell),                      &
                                   ndf_w2,                              &
                                   undf_w2,                             &
                                   map_w2(:,cell) )

      end do

    end if

  end subroutine invoke_fv_divergence

  !-------------------------------------------------------------------------------
  !> This kernel routine is in psykal_lite due to the use of cell_orientation
  !> variable which cannot be halo exchanged and also the depth of the halos to
  !> loop over includes all halo depths, i.e. all cells in the core and halo.
  subroutine invoke_extract_xy(x_field_out,y_field_out,w2_field_in,cell_orientation)

    use extract_x_kernel_mod,        only : extract_x_code
    use extract_y_kernel_mod,        only : extract_y_code
    use mesh_mod,                    only : mesh_type

    implicit none

    type(field_type), intent(inout)    :: x_field_out
    type(field_type), intent(inout)    :: y_field_out
    type(field_type), intent(in)       :: w2_field_in
    type(field_type), intent(in)       :: cell_orientation

    type(field_proxy_type)             :: cell_orientation_proxy
    type(field_proxy_type)             :: x_field_out_proxy
    type(field_proxy_type)             :: y_field_out_proxy
    type(field_proxy_type)             :: w2_field_in_proxy

    integer                            :: cell, nlayers
    integer                            :: ndf_w3, undf_w3
    integer                            :: ndf_w2, undf_w2
    integer, pointer                   :: map_w3(:,:) => null()
    integer, pointer                   :: map_w2(:,:) => null()
    integer                            :: halo_depth

    type(mesh_type), pointer           :: mesh => null()

    cell_orientation_proxy = cell_orientation%get_proxy()
    x_field_out_proxy      = x_field_out%get_proxy()
    y_field_out_proxy      = y_field_out%get_proxy()
    w2_field_in_proxy      = w2_field_in%get_proxy()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = x_field_out_proxy%vspace%get_ndf( )
    undf_w2 = x_field_out_proxy%vspace%get_undf()

    mesh   => x_field_out_proxy%vspace%get_mesh()
    map_w2 => x_field_out_proxy%vspace%get_whole_dofmap()
    map_w3 => cell_orientation_proxy%vspace%get_whole_dofmap()

    halo_depth = mesh%get_halo_depth()
    call w2_field_in_proxy%halo_exchange(depth=halo_depth)

    ! Extract the x-component of the W2 field
    do cell = 1, mesh%get_ncells_2d() ! Loop over core and halo cells

      call extract_x_code(  nlayers,                             &
                            cell_orientation_proxy%data,         &
                            w2_field_in_proxy%data,              &
                            x_field_out_proxy%data,              &
                            undf_w3,                             &
                            ndf_w3,                              &
                            map_w3(:,cell),                      &
                            undf_w2,                             &
                            ndf_w2,                              &
                            map_w2(:,cell) )
    end do

    ! Extract the y-component of the W2 field
    do cell = 1, mesh%get_ncells_2d() ! Loop over core and halo cells

      call extract_y_code(  nlayers,                             &
                            cell_orientation_proxy%data,         &
                            w2_field_in_proxy%data,              &
                            y_field_out_proxy%data,              &
                            undf_w3,                             &
                            ndf_w3,                              &
                            map_w3(:,cell),                      &
                            undf_w2,                             &
                            ndf_w2,                              &
                            map_w2(:,cell) )
    end do

  end subroutine invoke_extract_xy

  !-------------------------------------------------------------------------------
  ! Ticket #1156. Stephen Pring
  ! This code is implemented in psykal-lite because the cells to
  ! iterate over include all core cells and all halo cells. At present, the default
  ! iteration is over core cells and a halo depth of 1. The cosmic transport scheme
  ! uses a larger halo depth and this routine requires iteration over all values
  ! in the halo as well.
  subroutine invoke_cosmic_departure_wind(dep_wind_x,dep_wind_y,u_piola_x,u_piola_y,detj_at_w2,direction)
    use cosmic_departure_wind_kernel_mod, only: cosmic_departure_wind_code
    use mesh_mod,                         only: mesh_type
    use flux_direction_mod,               only: x_direction, y_direction
    use log_mod,                          only: log_event, LOG_LEVEL_ERROR

    implicit none

    type(field_type), intent(inout)      :: dep_wind_x, dep_wind_y
    type(field_type), intent(in)         :: u_piola_x, u_piola_y
    integer, intent(in)                  :: direction
    type(field_type), intent(in)         :: detj_at_w2

    type(field_proxy_type) :: dep_wind_x_p, dep_wind_y_p
    type(field_proxy_type) :: u_piola_x_p, u_piola_y_p
    type(field_proxy_type) :: detj_at_w2_p

    integer                 :: cell, nlayers
    integer                 :: ndf_w2
    integer                 :: undf_w2
    integer, pointer        :: map(:) => null()

    type(mesh_type), pointer :: mesh => null()
    integer :: halo_depth

    u_piola_x_p = u_piola_x%get_proxy()
    u_piola_y_p = u_piola_y%get_proxy()
    dep_wind_x_p = dep_wind_x%get_proxy()
    dep_wind_y_p = dep_wind_y%get_proxy()
    detj_at_w2_p = detj_at_w2%get_proxy()

    mesh => u_piola_x_p%vspace%get_mesh()
    halo_depth = mesh%get_halo_depth()
    call detj_at_w2_p%halo_exchange(depth=halo_depth)

    nlayers = u_piola_x_p%vspace%get_nlayers()

    ndf_w2  = u_piola_x_p%vspace%get_ndf()
    undf_w2 = u_piola_x_p%vspace%get_undf()


    ! NOTE: The default looping limits for this type of field would be
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and all halo cells.
    if (direction == x_direction) then
      do cell = 1,mesh%get_ncells_2d()
         map     => u_piola_x_p%vspace%get_cell_dofmap( cell )
         call cosmic_departure_wind_code( nlayers,                                  &
                                          dep_wind_x_p%data,                        &
                                          u_piola_x_p%data,                         &
                                          detj_at_w2_p%data,                        &
                                          ndf_w2, undf_w2, map,                     &
                                          direction                                 &
                                           )
      end do
    elseif (direction == y_direction) then
      do cell = 1,mesh%get_ncells_2d()
         map     => u_piola_y_p%vspace%get_cell_dofmap( cell )
         call cosmic_departure_wind_code( nlayers,                                  &
                                          dep_wind_y_p%data,                        &
                                          u_piola_y_p%data,                         &
                                          detj_at_w2_p%data,                        &
                                          ndf_w2, undf_w2, map,                     &
                                          direction                                 &
                                           )
      end do
    else
      call log_event("Direction incorrectly specified in invoke_cosmic_departure_wind",LOG_LEVEL_ERROR)
    end if

    call dep_wind_x_p%set_dirty()
    call dep_wind_y_p%set_dirty()

  end subroutine invoke_cosmic_departure_wind


  !-------------------------------------------------------------------------------
  ! Ticket #1156. Stephen Pring
  ! This code is implemented in psykal-lite because the cells to
  ! iterate over include all core cells and all halo cells. At present, the default
  ! iteration is over core cells and a halo depth of 1. The cosmic transport scheme
  ! uses a larger halo depth and this routine requires iteration over all values
  ! in the halo as well.
  subroutine invoke_correct_cosmic_wind(wind_x_out,                   &
                                        wind_y_out,                   &
                                        departure_wind_x_in,          &
                                        departure_wind_y_in,          &
                                        orientation_of_cells,         &
                                        direction)

    use correct_cosmic_wind_kernel_mod, only: correct_cosmic_wind_code
    use flux_direction_mod,             only: x_direction, y_direction
    use mesh_mod,                       only: mesh_type
    use log_mod,                        only: log_event, LOG_LEVEL_ERROR

    implicit none

    type(field_type), intent(inout) :: wind_x_out, wind_y_out
    type(field_type), intent(in)    :: departure_wind_x_in,departure_wind_y_in,orientation_of_cells
    integer,          intent(in)    :: direction

    type(field_proxy_type) :: wind_x_in_proxy, wind_y_in_proxy
    type(field_proxy_type) :: wind_x_out_proxy, wind_y_out_proxy
    type(field_proxy_type) :: orientation_proxy

    integer                 :: cell, nlayers
    integer                 :: ndf_w2, ndf_w3
    integer                 :: undf_w2, undf_w3
    integer, pointer        :: map_w3(:) => null()
    integer, pointer        :: map_w2(:) => null()

    type(mesh_type), pointer :: mesh => null()

    wind_x_in_proxy = departure_wind_x_in%get_proxy()
    wind_y_in_proxy = departure_wind_y_in%get_proxy()
    wind_x_out_proxy = wind_x_out%get_proxy()
    wind_y_out_proxy = wind_y_out%get_proxy()
    orientation_proxy = orientation_of_cells%get_proxy()

    nlayers = orientation_proxy%vspace%get_nlayers()
    ndf_w3  = orientation_proxy%vspace%get_ndf()
    undf_w3 = orientation_proxy%vspace%get_undf()

    ndf_w2  = wind_x_in_proxy%vspace%get_ndf( )
    undf_w2 = wind_x_in_proxy%vspace%get_undf()

    mesh => orientation_proxy%vspace%get_mesh()


    if (direction == x_direction) then
      do cell = 1, mesh%get_ncells_2d()
        map_w3 => orientation_proxy%vspace%get_cell_dofmap(cell)
        map_w2 => wind_x_in_proxy%vspace%get_cell_dofmap(cell)

        call correct_cosmic_wind_code(  nlayers,                        &
                                        wind_x_out_proxy%data,          &
                                        wind_x_in_proxy%data,           &
                                        orientation_proxy%data,         &
                                        undf_w2,                        &
                                        ndf_w2,                         &
                                        map_w2,                         &
                                        undf_w3,                        &
                                        ndf_w3,                         &
                                        map_w3,                         &
                                        direction )

      end do
    elseif (direction == y_direction) then
      do cell = 1, mesh%get_ncells_2d()
        map_w3 => orientation_proxy%vspace%get_cell_dofmap(cell)
        map_w2 => wind_x_in_proxy%vspace%get_cell_dofmap(cell)

        call correct_cosmic_wind_code(  nlayers,                        &
                                        wind_y_out_proxy%data,          &
                                        wind_y_in_proxy%data,           &
                                        orientation_proxy%data,         &
                                        undf_w2,                        &
                                        ndf_w2,                         &
                                        map_w2,                         &
                                        undf_w3,                        &
                                        ndf_w3,                         &
                                        map_w3,                         &
                                        direction )

      end do
    else
      call log_event("Direction incorrectly specified in invoke_correct_cosmic_wind",LOG_LEVEL_ERROR)
    end if

  end subroutine invoke_correct_cosmic_wind

  !----------------------------------------------------------------------------
  !> Handles passing double precision deltaT to the vertical_flux kernel.
  !>
  !> It is probably not the correct thing to do but it is expedient for getting
  !> the clock change on trunk.
  !>
  subroutine invoke_vertical_flux_kernel( mass_flux_z, &
                                          dep_pts,     &
                                          rho,         &
                                          a0, a1, a2,  &
                                          dt )

    use mesh_mod,                 only: mesh_type
    use vertical_flux_kernel_mod, only: vertical_flux_code

    implicit none

    real(kind=r_def), intent(in)    :: dt
    type(field_type), intent(inout) :: mass_flux_z
    type(field_type), intent(in)    :: dep_pts, rho, a0, a1, a2
    integer cell
    integer nlayers
    type(field_proxy_type) :: mass_flux_z_proxy, &
                              dep_pts_proxy,     &
                              rho_proxy,         &
                              a0_proxy, a1_proxy, a2_proxy
    integer, pointer :: map_w2(:,:) => null(), map_w3(:,:) => null()
    integer ndf_w2, undf_w2, ndf_w3, undf_w3
    type(mesh_type), pointer :: mesh => null()
    !
    ! initialise field and/or operator proxies
    !
    mass_flux_z_proxy = mass_flux_z%get_proxy()
    dep_pts_proxy = dep_pts%get_proxy()
    rho_proxy = rho%get_proxy()
    a0_proxy = a0%get_proxy()
    a1_proxy = a1%get_proxy()
    a2_proxy = a2%get_proxy()
    !
    ! initialise number of layers
    !
    nlayers = mass_flux_z_proxy%vspace%get_nlayers()
    !
    ! create a mesh object
    !
    mesh => mass_flux_z_proxy%vspace%get_mesh()
    !
    ! look-up dofmaps for each function space
    !
    map_w2 => mass_flux_z_proxy%vspace%get_whole_dofmap()
    map_w3 => rho_proxy%vspace%get_whole_dofmap()
    !
    ! initialise number of dofs for w2
    !
    ndf_w2 = mass_flux_z_proxy%vspace%get_ndf()
    undf_w2 = mass_flux_z_proxy%vspace%get_undf()
    !
    ! initialise number of dofs for w3
    !
    ndf_w3 = rho_proxy%vspace%get_ndf()
    undf_w3 = rho_proxy%vspace%get_undf()
    !
    ! call kernels and communication routines
    !
    if (mass_flux_z_proxy%is_dirty(depth=1)) then
      call mass_flux_z_proxy%halo_exchange(depth=1)
    end if
    !
    if (dep_pts_proxy%is_dirty(depth=1)) then
      call dep_pts_proxy%halo_exchange(depth=1)
    end if
    !
    if (rho_proxy%is_dirty(depth=1)) then
      call rho_proxy%halo_exchange(depth=1)
    end if
    !
    if (a0_proxy%is_dirty(depth=1)) then
      call a0_proxy%halo_exchange(depth=1)
    end if
    !
    if (a1_proxy%is_dirty(depth=1)) then
      call a1_proxy%halo_exchange(depth=1)
    end if
    !
    if (a2_proxy%is_dirty(depth=1)) then
      call a2_proxy%halo_exchange(depth=1)
    end if
    !
    do cell=1,mesh%get_last_halo_cell(1)
      !
      call vertical_flux_code(nlayers, mass_flux_z_proxy%data, dep_pts_proxy%data, rho_proxy%data, a0_proxy%data, a1_proxy%data, &
  &a2_proxy%data, dt, ndf_w2, undf_w2, map_w2(:,cell), ndf_w3, undf_w3, map_w3(:,cell))
    end do
    !
    ! set halos dirty/clean for fields modified in the above loop
    !
    call mass_flux_z_proxy%set_dirty()

  end subroutine invoke_vertical_flux_kernel

  !----------------------------------------------------------------------------
  ! This requires a stencil of horizontal cells for the operators
  ! see PSyclone #1103: https://github.com/stfc/PSyclone/issues/1103
  ! The LFRic infrastructure for this will be introduced in #2532
  subroutine invoke_helmholtz_operator_kernel_type(helmholtz_operator, hb_lumped_inv, stencil_depth, u_normalisation, div_star, &
                                                   t_normalisation, ptheta2v, compound_div, m3_exner_star, p3theta, m3_inv, &
                                                   w2_mask)
    use helmholtz_operator_kernel_mod, only: helmholtz_operator_code
    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod, only: stencil_cross
    use stencil_dofmap_mod, only: stencil_dofmap_type

    implicit none

    type(field_type), intent(in) :: helmholtz_operator(9), hb_lumped_inv, u_normalisation, t_normalisation, w2_mask
    type(operator_type), intent(in) :: div_star, ptheta2v, compound_div, m3_exner_star, p3theta, m3_inv
    integer(kind=i_def), intent(in) :: stencil_depth
    integer(kind=i_def) :: stencil_size
    integer(kind=i_def) cell
    integer(kind=i_def) nlayers
    type(operator_proxy_type) div_star_proxy, ptheta2v_proxy, compound_div_proxy, m3_exner_star_proxy, p3theta_proxy, m3_inv_proxy
    type(field_proxy_type) helmholtz_operator_proxy(9), hb_lumped_inv_proxy, u_normalisation_proxy, t_normalisation_proxy, &
                           w2_mask_proxy
    integer(kind=i_def), pointer :: map_w2(:,:) => null(), map_w3(:,:) => null(), map_wtheta(:,:) => null()
    integer(kind=i_def) ndf_w3, undf_w3, ndf_w2, undf_w2, ndf_wtheta, undf_wtheta
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) hb_lumped_inv_stencil_size
    integer(kind=i_def), pointer :: hb_lumped_inv_stencil_dofmap(:,:,:) => null()
    type(stencil_dofmap_type), pointer :: hb_lumped_inv_stencil_map => null()
    integer(kind=i_def) :: i
    integer(kind=i_def), allocatable :: cell_stencil(:)
    !
    ! Initialise field and/or operator proxies
    !
    helmholtz_operator_proxy(1) = helmholtz_operator(1)%get_proxy()
    helmholtz_operator_proxy(2) = helmholtz_operator(2)%get_proxy()
    helmholtz_operator_proxy(3) = helmholtz_operator(3)%get_proxy()
    helmholtz_operator_proxy(4) = helmholtz_operator(4)%get_proxy()
    helmholtz_operator_proxy(5) = helmholtz_operator(5)%get_proxy()
    helmholtz_operator_proxy(6) = helmholtz_operator(6)%get_proxy()
    helmholtz_operator_proxy(7) = helmholtz_operator(7)%get_proxy()
    helmholtz_operator_proxy(8) = helmholtz_operator(8)%get_proxy()
    helmholtz_operator_proxy(9) = helmholtz_operator(9)%get_proxy()
    hb_lumped_inv_proxy = hb_lumped_inv%get_proxy()
    u_normalisation_proxy = u_normalisation%get_proxy()
    div_star_proxy = div_star%get_proxy()
    t_normalisation_proxy = t_normalisation%get_proxy()
    ptheta2v_proxy = ptheta2v%get_proxy()
    compound_div_proxy = compound_div%get_proxy()
    m3_exner_star_proxy = m3_exner_star%get_proxy()
    p3theta_proxy = p3theta%get_proxy()
    m3_inv_proxy = m3_inv%get_proxy()
    w2_mask_proxy = w2_mask%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = helmholtz_operator_proxy(1)%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => helmholtz_operator_proxy(1)%vspace%get_mesh()
    !
    ! Initialise stencil dofmaps
    !
    hb_lumped_inv_stencil_map => hb_lumped_inv_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS,stencil_depth)
    hb_lumped_inv_stencil_dofmap => hb_lumped_inv_stencil_map%get_whole_dofmap()
    hb_lumped_inv_stencil_size = hb_lumped_inv_stencil_map%get_size()
    !
    ! Look-up dofmaps for each function space
    !
    map_w3 => helmholtz_operator_proxy(1)%vspace%get_whole_dofmap()
    map_w2 => hb_lumped_inv_proxy%vspace%get_whole_dofmap()
    map_wtheta => t_normalisation_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for w3
    !
    ndf_w3 = helmholtz_operator_proxy(1)%vspace%get_ndf()
    undf_w3 = helmholtz_operator_proxy(1)%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = hb_lumped_inv_proxy%vspace%get_ndf()
    undf_w2 = hb_lumped_inv_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for wtheta
    !
    ndf_wtheta = t_normalisation_proxy%vspace%get_ndf()
    undf_wtheta = t_normalisation_proxy%vspace%get_undf()
    !
    ! Call kernels and communication routines
    !
    if (hb_lumped_inv_proxy%is_dirty(depth=1)) then
      call hb_lumped_inv_proxy%halo_exchange(depth=1)
    end if
    if (u_normalisation_proxy%is_dirty(depth=1)) then
      call u_normalisation_proxy%halo_exchange(depth=1)
    end if
    if (w2_mask_proxy%is_dirty(depth=1)) then
      call w2_mask_proxy%halo_exchange(depth=1)
    end if
    !
    ! Create cell stencil of the correct size
    stencil_size =  1 + 4*stencil_depth
    allocate( cell_stencil( stencil_size ) )

    !$omp parallel default(shared), private(cell, cell_stencil, i)
    !$omp do schedule(static)
    do cell=1,mesh%get_last_edge_cell()
      !
      ! Populate cell_stencil array used for operators
      ! (this is the id of each cell in the stencil)
      cell_stencil(1) = cell
      do i = 1,4
        cell_stencil(i+1) = mesh%get_cell_next(i, cell)
      end do
      call helmholtz_operator_code(stencil_size,                     &
                                   cell_stencil, nlayers,            &
                                   helmholtz_operator_proxy(1)%data, &
                                   helmholtz_operator_proxy(2)%data, &
                                   helmholtz_operator_proxy(3)%data, &
                                   helmholtz_operator_proxy(4)%data, &
                                   helmholtz_operator_proxy(5)%data, &
                                   helmholtz_operator_proxy(6)%data, &
                                   helmholtz_operator_proxy(7)%data, &
                                   helmholtz_operator_proxy(8)%data, &
                                   helmholtz_operator_proxy(9)%data, &
                                   hb_lumped_inv_proxy%data, &
                                   hb_lumped_inv_stencil_size, &
                                   hb_lumped_inv_stencil_dofmap(:,:,cell), &
                                   u_normalisation_proxy%data, &
                                   div_star_proxy%ncell_3d, &
                                   div_star_proxy%local_stencil, &
                                   t_normalisation_proxy%data, &
                                   ptheta2v_proxy%ncell_3d, &
                                   ptheta2v_proxy%local_stencil, &
                                   compound_div_proxy%ncell_3d, &
                                   compound_div_proxy%local_stencil, &
                                   m3_exner_star_proxy%ncell_3d, &
                                   m3_exner_star_proxy%local_stencil, &
                                   p3theta_proxy%ncell_3d, &
                                   p3theta_proxy%local_stencil, &
                                   m3_inv_proxy%ncell_3d, &
                                   m3_inv_proxy%local_stencil, &
                                   w2_mask_proxy%data, &
                                   ndf_w3, undf_w3, map_w3(:,cell), &
                                   ndf_w2, undf_w2, map_w2(:,cell), &
                                   ndf_wtheta, undf_wtheta, map_wtheta(:,cell))
    end do
    !$omp end do
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    !$omp master
    call helmholtz_operator_proxy(1)%set_dirty()
    call helmholtz_operator_proxy(2)%set_dirty()
    call helmholtz_operator_proxy(3)%set_dirty()
    call helmholtz_operator_proxy(4)%set_dirty()
    call helmholtz_operator_proxy(5)%set_dirty()
    call helmholtz_operator_proxy(6)%set_dirty()
    call helmholtz_operator_proxy(7)%set_dirty()
    call helmholtz_operator_proxy(8)%set_dirty()
    call helmholtz_operator_proxy(9)%set_dirty()
    !$omp end master
    !
    !$omp end parallel
    !
  end subroutine invoke_helmholtz_operator_kernel_type

end module psykal_lite_mod
