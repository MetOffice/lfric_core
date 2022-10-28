!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the flux in x using 1D PPM.
!> @details This kernel computes the flux in the x direction. PPM is used
!!          to compute the subgrid reconstruction of the form a0 + a1 x + a2 x^2,
!!          and this is integrated between the flux point and the departure point.
!!          For CFL > 1 the field values are summed between the flux point and
!!          the departure cell. As this is used for the final steps of the FFSL
!!          transport scheme the flux is computed using field_y. At cubed sphere
!!          panel edges the correct direction must be used, hence field_x is also
!!          an input.
!!
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module ffsl_flux_final_x_kernel_mod

  use argument_mod,       only : arg_type,                 &
                                 GH_FIELD, GH_REAL,        &
                                 CELL_COLUMN, GH_WRITE,    &
                                 GH_READ, GH_SCALAR,       &
                                 STENCIL, X1D, GH_INTEGER, &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_def, i_def
  use fs_continuity_mod,  only : W3, W2
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_final_x_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                                      &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2),                                      & ! flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),                        & ! field_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),                        & ! field_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(X1D)), & ! panel_id
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                                      & ! dep_pts
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                      & ! order
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                      & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                                       & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_final_x_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_flux_final_x_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] flux              The output flux in x
  !> @param[in]     field_x           Field from x direction
  !> @param[in]     stencil_size_x    Local length of field_x W3 stencil
  !> @param[in]     stencil_map_x     Dofmap for the field_x stencil
  !> @param[in]     field_y           Field from y direction
  !> @param[in]     stencil_size_y    Local length of field_y W3 stencil
  !> @param[in]     stencil_map_y     Dofmap for the field_y stencil
  !> @param[in]     panel_id          Panel ID of cell
  !> @param[in]     stencil_size_p    Local length of panel ID stencil
  !> @param[in]     stencil_map_p     Dofmap for the panel ID stencil
  !> @param[in]     dep_pts           Departure points in x
  !> @param[in]     order             Order of reconstruction
  !> @param[in]     extent_size       Stencil extent needed for the LAM edge
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_w2            Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2           Number of unique degrees of freedom for W2
  !> @param[in]     map_w2            Map for W2
  !> @param[in]     ndf_w3            Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3           Number of unique degrees of freedom for W3
  !> @param[in]     map_w3            Map for W3
  !> @param[in]     ndf_wp            Number of degrees of freedom for panel ID
  !!                                  function space per cell
  !> @param[in]     undf_wp           Number of unique degrees of freedom for
  !!                                  panel ID function space
  !> @param[in]     map_wp            Map for panel ID function space

  subroutine ffsl_flux_final_x_code( nlayers,        &
                                     flux,           &
                                     field_x,        &
                                     stencil_size_x, &
                                     stencil_map_x,  &
                                     field_y,        &
                                     stencil_size_y, &
                                     stencil_map_y,  &
                                     panel_id,       &
                                     stencil_size_p, &
                                     stencil_map_p,  &
                                     dep_pts,        &
                                     order,          &
                                     extent_size,    &
                                     dt,             &
                                     ndf_w2,         &
                                     undf_w2,        &
                                     map_w2,         &
                                     ndf_w3,         &
                                     undf_w3,        &
                                     map_w3,         &
                                     ndf_wp,         &
                                     undf_wp,        &
                                     map_wp )

    use subgrid_rho_mod,            only: second_order_coeffs
    use cosmic_flux_mod,            only: frac_and_int_part,       &
                                          calc_integration_limits, &
                                          return_part_mass,        &
                                          get_index_negative,      &
                                          get_index_positive
    use ffsl_cubed_sphere_edge_mod, only: get_local_rho_y

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2
    integer(kind=i_def), intent(in) :: ndf_w2
    integer(kind=i_def), intent(in) :: undf_wp
    integer(kind=i_def), intent(in) :: ndf_wp
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_y
    integer(kind=i_def), intent(in) :: stencil_size_p

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_wp), intent(in) :: map_wp
    integer(kind=i_def), dimension(ndf_w3,stencil_size_x), intent(in) :: stencil_map_x
    integer(kind=i_def), dimension(ndf_w3,stencil_size_y), intent(in) :: stencil_map_y
    integer(kind=i_def), dimension(ndf_wp,stencil_size_p), intent(in) :: stencil_map_p

    ! Arguments: Fields
    real(kind=r_def), dimension(undf_w2), intent(inout) :: flux
    real(kind=r_def), dimension(undf_w3), intent(in)    :: field_x
    real(kind=r_def), dimension(undf_w3), intent(in)    :: field_y
    real(kind=r_def), dimension(undf_wp), intent(in)    :: panel_id
    real(kind=r_def), dimension(undf_w2), intent(in)    :: dep_pts
    integer(kind=i_def), intent(in)                     :: order
    integer(kind=i_def), intent(in)                     :: extent_size
    real(kind=r_def), intent(in)                        :: dt

    ! Variables for flux calculation
    real(kind=r_def) :: mass_total
    real(kind=r_def) :: departure_dist
    real(kind=r_def) :: fractional_distance
    real(kind=r_def) :: mass_frac
    real(kind=r_def) :: mass_from_whole_cells
    real(kind=r_def) :: left_integration_limit
    real(kind=r_def) :: right_integration_limit

    ! Local fields
    real(kind=r_def)    :: field_local(1:stencil_size_y)
    real(kind=r_def)    :: field_x_local(1:stencil_size_x)
    real(kind=r_def)    :: field_y_local(1:stencil_size_y)
    integer(kind=i_def) :: ipanel_local(1:stencil_size_p)

    ! PPM coefficients
    real(kind=r_def)    :: coeffs(1:3)

    ! DOFs
    integer(kind=i_def) :: local_dofs(1:2)
    integer(kind=i_def) :: dof_iterator

    ! Indices
    integer(kind=i_def) :: n_cells_to_sum
    integer(kind=i_def) :: ind_lo, ind_hi
    integer(kind=i_def) :: k, ii, jj

    ! Stencils
    integer(kind=i_def) :: stencil_half, stencil_size, lam_edge_size

    ! Stencil has order e.g.        | 5 | 4 | 3 | 2 | 1 | 6 | 7 | 8 | 9 | for extent 4
    ! Local fields have order e.g.  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | for extent 4
    ! Fluxes calculated for centre cell, e.g. cell 1 for stencil, cell 5 for local

    ! x-direction
    local_dofs = (/ 1, 3 /)

    ! Use stencil_size_y as each stencil size should be equal
    stencil_size = stencil_size_y
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 2_i_def*extent_size+1_i_def

    if (lam_edge_size > stencil_size) then

      ! At edge of LAM, so set output to zero
      do k = 0,nlayers-1
        do dof_iterator = 1,2
         flux( map_w2(local_dofs(dof_iterator)) + k ) = 0.0_r_def
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes

      ! Set up local panel ID
      do jj = 1, stencil_half
        ipanel_local(jj) = int(panel_id(stencil_map_p(1,stencil_half+1-jj)), i_def)
      end do
      do jj = stencil_half+1, stencil_size
        ipanel_local(jj) = int(panel_id(stencil_map_p(1,jj)), i_def)
      end do
      ! Initialise field_local to zero
      field_local(1:stencil_size) = 0.0_r_def
      coeffs(1:3) = 0.0_r_def

      do k = 0,nlayers-1

        do dof_iterator = 1,2
          ! Loop over the x direction dofs to compute flux at each dof

          ! Get the departure distance
          departure_dist = dep_pts( map_w2(local_dofs(dof_iterator)) + k )

          ! Calculates number of cells of interest and fraction of a cell to add.
          call frac_and_int_part(departure_dist,n_cells_to_sum,fractional_distance)

          ! Get local field values - this will depend on panel ID at cubed sphere edges
          do jj = 1, stencil_half
            field_y_local(jj) = field_y(stencil_map_y(1,stencil_half+1-jj) + k)
            field_x_local(jj) = field_x(stencil_map_x(1,stencil_half+1-jj) + k)
          end do
          do jj = stencil_half+1, stencil_size
            field_y_local(jj) = field_y(stencil_map_y(1,jj) + k)
            field_x_local(jj) = field_x(stencil_map_x(1,jj) + k)
          end do
          call get_local_rho_y(field_local,field_x_local,field_y_local,ipanel_local,stencil_size,stencil_half)

          ! Get a0, a1, a2 in the required cell and build up whole cell part
          mass_from_whole_cells = 0.0_r_def
          if (departure_dist >= 0.0_r_def ) then
            call get_index_positive(ind_lo,ind_hi,n_cells_to_sum,dof_iterator,stencil_size,stencil_half)
            do ii = 1, n_cells_to_sum-1
              mass_from_whole_cells = mass_from_whole_cells + field_local(stencil_half - (2-dof_iterator) - (ii-1) )
            end do
          else
            call get_index_negative(ind_lo,ind_hi,n_cells_to_sum,dof_iterator,stencil_size,stencil_half)
            do ii = 1, n_cells_to_sum-1
              mass_from_whole_cells = mass_from_whole_cells + field_local(stencil_half + (dof_iterator-1) + (ii-1) )
            end do
          end if
          if ( order == 0 ) then
            ! Piecewise constant reconstruction
            coeffs(1) = field_local(ind_lo+2)
          else
            ! Piecewise parabolic reconstruction
            call second_order_coeffs( field_local(ind_lo:ind_hi), coeffs, .false., .false.)
          end if

          ! Calculates the left and right integration limits for the fractional cell.
          call calc_integration_limits( departure_dist,         &
                                        fractional_distance,    &
                                        left_integration_limit, &
                                        right_integration_limit )

          ! Compute fractional flux
          mass_frac = return_part_mass(3,coeffs,left_integration_limit,right_integration_limit)

          ! Get total flux, i.e. fractional part + whole cell part
          mass_total = mass_from_whole_cells + mass_frac

          ! Assign to flux variable and divide by dt to get the correct form
          flux(map_w2(local_dofs(dof_iterator)) + k) =  sign(1.0_r_def,departure_dist) * mass_total / dt

        end do

      end do

    end if

  end subroutine ffsl_flux_final_x_code

end module ffsl_flux_final_x_kernel_mod
