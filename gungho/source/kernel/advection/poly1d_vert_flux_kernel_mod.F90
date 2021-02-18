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
!>          This method is only valid for the lowest order elements based on
!>          a reference cube.
module poly1d_vert_flux_kernel_mod

use argument_mod,      only : arg_type, func_type,           &
                              reference_element_data_type,   &
                              GH_FIELD, GH_INTEGER,          &
                              GH_INC, GH_READ,               &
                              GH_BASIS, CELLS, GH_EVALUATOR, &
                              outward_normals_to_vertical_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_flux_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                   &
       arg_type(GH_FIELD,   GH_INC,   W2),                              &
       arg_type(GH_FIELD,   GH_READ,  W2),                              &
       arg_type(GH_FIELD,   GH_READ,  W3),                              &
       arg_type(GH_FIELD,   GH_READ,  W3),                              &
       arg_type(GH_INTEGER, GH_READ),                                   &
       arg_type(GH_INTEGER, GH_READ)                                    &
       /)
  type(func_type) :: meta_funcs(1) = (/                                 &
       func_type(W2, GH_BASIS)                                          &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/   &
       reference_element_data_type( outward_normals_to_vertical_faces ) &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: poly1d_vert_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_flux_code

contains

!> @brief Computes the vertical fluxes for a tracer density.
!! @param[in]  nlayers Number of layers
!! @param[out] flux Mass flux field to compute
!! @param[in]  wind Wind field
!! @param[in]  density Tracer density
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux &
!!                     wind fields
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at w2 nodes
!! @param[in]  ndf_w3 Number of degrees of freedom per cell
!! @param[in]  undf_w3 Number of unique degrees of freedom for the density field
!! @param[in]  map_w3 Cell dofmaps for the density space
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  nfaces_re_v Number of vertical faces (used by PSyclone to size
!!                         coeff array)
!! @param[in]  outward_normals_to_vertical_faces Vector of normals to the
!!                                               reference element vertical
!!                                               "outward faces"
!! @param[in]  logspace If true, then perform interpolation in log space
subroutine poly1d_vert_flux_code( nlayers,                           &
                                  flux,                              &
                                  wind,                              &
                                  density,                           &
                                  coeff,                             &
                                  ndf_w2,                            &
                                  undf_w2,                           &
                                  map_w2,                            &
                                  basis_w2,                          &
                                  ndf_w3,                            &
                                  undf_w3,                           &
                                  map_w3,                            &
                                  global_order,                      &
                                  nfaces_re_v,                       &
                                  outward_normals_to_vertical_faces, &
                                  logspace )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_re_v

  real(kind=r_def), dimension(undf_w2), intent(out)  :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)   :: density

  real(kind=r_def), dimension(global_order+1, nfaces_re_v, undf_w3), intent(in) :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  real(kind=r_def), intent(in)  :: outward_normals_to_vertical_faces(:,:)

  logical, intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, df, ij, p, stencil, order, id, &
                                    m, ijkp, boundary_offset, vertical_order
  integer(kind=i_def)            :: vert_offset
  real(kind=r_def)               :: direction
  real(kind=r_def), dimension(2) :: v_dot_n
  real(kind=r_def)               :: polynomial_density

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! TODO: Set the offset for vertical faces for the reference cube. This should
  ! not be hardcoded as it depends on element type (cube vs prism). Instead it
  ! should really be passed from the PSy layer for the kernel to be
  ! element-agnostic
  vert_offset = 4_i_def

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  do id = 1,nfaces_re_v
    df = id + vert_offset
    v_dot_n(id) =  dot_product(basis_w2(:,df,df), outward_normals_to_vertical_faces(:,id))
  end do

  ij = map_w3(1)

  ! Vertical flux computation
  do k = 0, nlayers - 1
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))

    boundary_offset = 0
    if ( order > 0 ) then
      if ( k == 0 )           boundary_offset =  1
      if ( k == nlayers - 1 ) boundary_offset = -1
    end if

    do id = 1,nfaces_re_v

      df = id + vert_offset
      ! Check if this is the upwind cell
      direction = wind(map_w2(df) + k )*v_dot_n(id)
      if ( direction >= 0.0_r_def ) then
        if (logspace) then
          ! Interpolate log(tracer)
          ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
          !                 = tracer_1**c_1*tracer_2**c_2...
          ! Note that we further take the absolute value before raising to the
          ! fractional power. This code should only be used for a positive
          ! quantity, but adding in the abs ensures no errors are thrown
          ! if negative numbers are passed through in redundant calculations
          ! in the haloes
          polynomial_density = 1.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            polynomial_density = polynomial_density &
                               * abs(density( ijkp ))**coeff( p, id, ij+k )
          end do
        else
          polynomial_density = 0.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            polynomial_density = polynomial_density &
                               + density( ijkp )*coeff( p, id, ij+k )
          end do
        end if
        flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
      end if
    end do
  end do

  ! Boundary conditions
  ! When computing a flux the wind = 0 so flux is zero
  ! For averaging the wind = +/-1 so the resulting average is non-zero
  ! Make sure we always compute the boundary values
  order = min(vertical_order, 2)
  boundary_offset = 0

  ! The boundary values only have a cell on one side of them so the
  ! reconstruction above may not have performed. Therefore we ensure
  ! that the values at the boundaries are computed

  k = 0
  id = 1
  if ( order > 0 ) boundary_offset = 1
  df = id + vert_offset
  if (logspace) then
    ! Interpolate log(tracer)
    ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
    !                 = tracer_1**c_1*tracer_2**c_2...
    polynomial_density = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      polynomial_density = polynomial_density &
                         * abs(density( ijkp ))**coeff( p, id, ij+k )
    end do
  else
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      polynomial_density = polynomial_density &
                         + density( ijkp )*coeff( p, id, ij+k )
    end do
  end if
  flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density

  k = nlayers-1
  id = 2
  if ( order > 0 ) boundary_offset =  -1
  df = id + vert_offset
  if (logspace) then
    ! Interpolate log(tracer)
    ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
    !                 = tracer_1**c_1*tracer_2**c_2...
    polynomial_density = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      polynomial_density = polynomial_density &
                         * abs(density( ijkp ))**coeff( p, id, ij+k )
    end do
  else
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      polynomial_density = polynomial_density &
                         + density( ijkp )*coeff( p, id, ij+k )
    end do
  end if
  flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density

  deallocate( smap )

end subroutine poly1d_vert_flux_code

end module poly1d_vert_flux_kernel_mod
