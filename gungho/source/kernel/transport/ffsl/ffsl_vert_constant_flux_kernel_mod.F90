!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the mass flux for FFSL using PCM in the z direction.
!> @details The flux form semi-Lagrangian (FFSL) scheme updates density in
!!          the x, y and z directions separately.
!!          This code calculates the mass flux in the z direction (vertical) using
!!          the piecewise constant method (PCM).
!!          The scheme calculates the total mass which has flowed through a cell
!!          wall and divides by dt (timestep length) to obtain a rate of flow.
!!          The scheme is able to handle Courant numbers larger than 1.

module ffsl_vert_constant_flux_kernel_mod

use argument_mod,      only : arg_type,            &
                              GH_FIELD, GH_SCALAR, &
                              GH_READ, GH_WRITE,   &
                              GH_REAL, CELL_COLUMN
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_vert_constant_flux_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/              &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2), & ! flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2), & ! dep_pts
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3), & ! rho
       arg_type(GH_SCALAR, GH_REAL, GH_READ)       & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_vert_constant_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_vert_constant_flux_code

contains

!> @brief Kernel which computes the mass flux for FFSL in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      Flux values which are calculated
!> @param[in]     dep_pts   Departure points
!> @param[in]     rho       Density values in W3
!> @param[in]     deltaT    Length of a timestep
!> @param[in]     ndf_w2    Number of degrees of freedom per cell
!> @param[in]     undf_w2   Number of unique degrees of freedom
!> @param[in]     map_w2    Dofmap for the cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom
!> @param[in]     map_w3    Dofmap for the cell at the base of the column
subroutine ffsl_vert_constant_flux_code( nlayers, &
                                         flux,    &
                                         dep_pts, &
                                         rho,     &
                                         deltaT,  &
                                         ndf_w2,  &
                                         undf_w2, &
                                         map_w2,  &
                                         ndf_w3,  &
                                         undf_w3, &
                                         map_w3 )

  use cosmic_flux_mod,    only : frac_and_int_part,       &
                                 calc_integration_limits, &
                                 return_part_mass,        &
                                 calc_local_vertical_index

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                       :: nlayers
  integer(kind=i_def), intent(in)                       :: ndf_w3
  integer(kind=i_def), intent(in)                       :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)    :: map_w3
  real(kind=r_def), dimension(undf_w3), intent(in)      :: rho
  integer(kind=i_def), intent(in)                       :: ndf_w2
  integer(kind=i_def), intent(in)                       :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(inout)   :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)      :: dep_pts
  real(kind=r_def),                     intent(in)      :: deltaT

  ! Internal variables
  integer(kind=i_def) :: k
  integer(kind=i_def) :: ii
  integer(kind=i_def) :: edge
  integer(kind=i_def) :: n_cells_to_sum
  integer(kind=i_def) :: df
  integer(kind=i_def) :: local_density_index(nlayers)

  real(kind=r_def)    :: mass_total
  real(kind=r_def)    :: departure_dist
  real(kind=r_def)    :: fractional_distance
  real(kind=r_def)    :: mass_frac
  real(kind=r_def)    :: mass_from_whole_cells
  real(kind=r_def)    :: left_integration_limit
  real(kind=r_def)    :: right_integration_limit
  real(kind=r_def)    :: subgrid_coeffs(3)
  real(kind=r_def)    :: rho_local(nlayers)

  ! Flux boundary conditions
  k = 0
  df = 5
  ! Bottom boundary condition, zero flux
  flux(map_w2(df)+k) = 0.0_r_def
  k = nlayers-1
  df = 6
  ! Top boundary condition, zero flux
  flux(map_w2(df)+k) = 0.0_r_def

  ! Initialise to zero
  local_density_index = 0_i_def
  rho_local = 0.0_r_def

  do k=0,nlayers-2
    departure_dist = dep_pts(map_w2(df)+k)

    ! Calculates number of cells of interest and fraction of a cell to add.
    ! There is a test within this function to detect if the departure points
    ! go out of bounds.
    call frac_and_int_part(departure_dist,n_cells_to_sum,fractional_distance)

    edge = 1
    call calc_local_vertical_index(local_density_index,departure_dist,n_cells_to_sum,k,nlayers,edge)

    do ii=1,n_cells_to_sum
      rho_local(ii) = rho(map_w3(1)+local_density_index(ii))
    end do

    ! Calculates the left and right integration limits for the fractional cell
    call calc_integration_limits( departure_dist,         &
                                  fractional_distance,    &
                                  left_integration_limit, &
                                  right_integration_limit )

    mass_from_whole_cells = sum(rho_local(1:n_cells_to_sum-1))

    ! Constant reconstruction, therefore set a0 = rho and a1=a2=0
    subgrid_coeffs = (/ rho_local(n_cells_to_sum), &
                        0.0_r_def,                 &
                        0.0_r_def /)

    mass_frac = return_part_mass(3,subgrid_coeffs,left_integration_limit,right_integration_limit)

    ! Calculate the total mass, which, depending on the departure distance
    ! is the sum of mass from whole cells, i.e. the integer part of the
    ! departure distance, and the fractional part of the departure distance.
    mass_total = mass_from_whole_cells + mass_frac

    flux(map_w2(df)+k) = sign(1.0_r_def,departure_dist)*mass_total/deltaT

  end do

end subroutine ffsl_vert_constant_flux_code

end module ffsl_vert_constant_flux_kernel_mod
