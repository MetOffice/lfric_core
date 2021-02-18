!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes vertical fluxes through fitting a high order 1D
!>        upwind reconstruction.
!> @details Computes the flux for a tracer density field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used
!>          for the polynomial is centred on the upwind cell.
!>          Near the boundaries the order of reconstruction may be reduced
!>          if there are not enough points to compute desired order.
!>          This method is only valid for lowest order elements.
module poly1d_vert_adv_kernel_mod

use argument_mod,         only : arg_type, func_type,   &
                                 GH_FIELD, GH_INTEGER,  &
                                 GH_READWRITE, GH_READ, &
                                 GH_BASIS, CELLS, GH_EVALUATOR
use constants_mod,        only : r_def, i_def
use fs_continuity_mod,    only : W2, Wtheta
use kernel_mod,           only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_READWRITE, Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      W2),                         &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: poly1d_vert_adv_code
end type

type, public, extends(kernel_type) :: poly1d_vert_adv_old_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_READWRITE, Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      W2),                         &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_FIELD,   GH_READ,      Wtheta),                     &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: poly1d_vert_adv_old_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_adv_code
public poly1d_vert_adv_old_code

contains

!> @brief Computes the vertical fluxes for a tracer density.
!! @param[in]  nlayers Number of layers
!! @param[in,out] advective Advective update to increment
!! @param[in]  wind Wind field
!! @param[in]  tracer Tracer field to advect
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  map_wt Cell dofmaps for the tracer space
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux &
!!                     wind fields
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at wt nodes
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  nfaces_v Number of vertical faces (used by PSyclone to size
!!                      coeff array)
!! @param[in]  logspace If true, then perform interpolation in log space
subroutine poly1d_vert_adv_code( nlayers,              &
                                 advective,            &
                                 wind,                 &
                                 tracer,               &
                                 coeff,                &
                                 ndf_wt,               &
                                 undf_wt,              &
                                 map_wt,               &
                                 ndf_w2,               &
                                 undf_w2,              &
                                 map_w2,               &
                                 global_order,         &
                                 nfaces_v,             &
                                 logspace )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_v

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer

  real(kind=r_def), dimension(global_order+1, nfaces_v, undf_wt), intent(in) :: coeff

  logical, intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, ij, p, stencil, order, &
                                    m, ijkp, direction, boundary_offset, &
                                    vertical_order, km

  real(kind=r_def)                         :: w, tracer_p, tracer_m

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  ij = map_wt(1)

  ! Reconstruct upwind tracer at cell centres with stencil
  ! direction determined by w at Wtheta points
  do k = 1, nlayers - 1

    ! Check if this is the upwind cell
    w = sign(1.0_r_def,wind(map_w2(5) + k ))
    ! w > 0, direction = 0
    ! w < 0, direction = 1
    direction = int( 0.5_r_def*(abs(w)-w),i_def )

    ! Compute value at W3 point in cell below this cell
    km = k - 1
    order = min(vertical_order, min(2*(km+1), 2*(nlayers-1 - (km-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( km == 0 .and. direction == 0 ) order = 1

    if (logspace) then
      ! Interpolate log(tracer)
      ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
      !                 = tracer_1**c_1*tracer_2**c_2...
      ! Note that we further take the absolute value before raising to the
      ! fractional power. This code should only be used for a positive
      ! quantity, but adding in the abs ensures no errors are thrown
      ! if negative numbers are passed through in redundant calculations
      ! in the haloes
      tracer_m = 1.0_r_def

      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        tracer_m = tracer_m  * abs(tracer( ijkp ))**coeff( p, direction+1, ij+km )
      end do
    else
      tracer_m = 0.0_r_def

      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        tracer_m = tracer_m  + tracer( ijkp )*coeff( p, direction+1, ij+km )
      end do
    end if

    ! Compute value at W3 point in this cell
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( k == nlayers - 1 .and. direction == 1 ) order = 1


    ! Offset at the top boundary to ensure we use the correct entries
    boundary_offset = 0
    if ( k == nlayers - 1 .and. direction == 1)  boundary_offset = -1

    if (logspace) then
      ! Interpolate log(tracer)
      ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
      !                 = tracer_1**c_1*tracer_2**c_2...
      tracer_p = 1.0_r_def

      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        tracer_p = tracer_p  * abs(tracer( ijkp ))**coeff( p, direction+1, ij+k )
      end do
    else
      tracer_p = 0.0_r_def

      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        tracer_p = tracer_p  + tracer( ijkp )*coeff( p, direction+1, ij+k )
      end do
    end if

    ! Compute advective increment
    advective(map_wt(1) + k ) = advective(map_wt(1) + k ) &
                              + wind(map_w2(5) + k )      &
                              *(tracer_p - tracer_m)
  end do

  deallocate( smap )
end subroutine poly1d_vert_adv_code

!> @brief Computes the vertical fluxes for a tracer density.
!>        Uses the old method where the averaged velocity
!>        at cell centres is used to determine the upwind
!>        direction.
!! @param[in]  nlayers Number of layers
!! @param[in,out] advective Advective update to increment
!! @param[in]  wind Wind field
!! @param[in]  tracer Tracer field to advect
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  map_wt Cell dofmaps for the tracer space
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux &
!!                     wind fields
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at wt nodes
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  nfaces_v Number of vertical faces (used by PSyclone to size
!!                      coeff array)
!! @param[in]  logspace If true, then perform interpolation in log space
subroutine poly1d_vert_adv_old_code( nlayers,              &
                                     advective,            &
                                     wind,                 &
                                     tracer,               &
                                     coeff,                &
                                     ndf_wt,               &
                                     undf_wt,              &
                                     map_wt,               &
                                     ndf_w2,               &
                                     undf_w2,              &
                                     map_w2,               &
                                     global_order,         &
                                     nfaces_v,             &
                                     logspace )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_v

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer

  real(kind=r_def), dimension(global_order+1, nfaces_v, undf_wt), intent(in) :: coeff

  logical, intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, ij, p, stencil, order, &
                                    m, ijkp, direction, boundary_offset, &
                                    vertical_order

  real(kind=r_def)                         :: w
  real(kind=r_def), dimension(0:nlayers-1) :: polynomial_tracer

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  ij = map_wt(1)

  ! Reconstruct upwind tracer at cell centres with stencil
  ! direction determined by w at W3 points
  do k = 0, nlayers - 1

    ! Check if this is the upwind cell
    w = sign(1.0_r_def,wind(map_w2(5) + k ) + wind(map_w2(6) + k ))
    ! w > 0, direction = 0
    ! w < 0, direction = 1
    direction = int( 0.5_r_def*(abs(w)-w),i_def )

    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( ( k == 0           .and. direction == 0) .or. &
         ( k == nlayers - 1 .and. direction == 1) ) order = 1

    ! Offset at the top boundary to ensure we use the correct entries
    boundary_offset = 0
    if ( k == nlayers - 1 .and. direction == 1)  boundary_offset = -1

    if (logspace) then
      ! Interpolate log(tracer)
      ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
      !                 = tracer_1**c_1*tracer_2**c_2...
      polynomial_tracer(k) = 1.0_r_def

      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        polynomial_tracer(k) = polynomial_tracer(k) &
                             * abs(tracer( ijkp ))**coeff( p, direction+1, ij+k )
      end do
    else
      polynomial_tracer(k) = 0.0_r_def

      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        polynomial_tracer(k) = polynomial_tracer(k) &
                             + tracer( ijkp )*coeff( p, direction+1, ij+k )
      end do
    end if
  end do

  do k = 1, nlayers - 1
    advective(map_wt(1) + k ) = advective(map_wt(1) + k ) &
                              + wind(map_w2(5) + k )      &
                              *(polynomial_tracer(k)      &
                              - polynomial_tracer(k-1))
  end do

  deallocate( smap )
end subroutine poly1d_vert_adv_old_code

end module poly1d_vert_adv_kernel_mod
