!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the fluxes through fitting a high order 1D
!>        upwind reconstruction
!> @details Compute the flux for a tracer density field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used for the
!>          polynomial is centred, with an upwind bias if an odd number of 
!>          points are used, therefore for an upwind shceme an even ordered 
!>          polynomial should be used (this means odd ordered for the tracer
!>          integral. In the vertical the order is reduced 
!>          near the boundaries depending on the number of points available.
!>          This method is only valid for lowest order elements 
module sample_poly_flux_kernel_mod

use argument_mod,  only : arg_type, func_type,                  &
                          GH_FIELD, GH_WRITE, GH_READ,          &
                          W2, W3, GH_BASIS, CELLS,              &
                          GH_EVALUATOR, EVALUATOR
use constants_mod, only : r_def, i_def
use kernel_mod,    only : kernel_type
use reference_element_mod, only: out_face_normal


implicit none

! Precomputed operators, these are the same for all model columns
real(kind=r_def), allocatable,    private :: coeff_matrix(:,:)
real(kind=r_def), allocatable,    private :: coeff(:), density_stencil(:)
integer(kind=i_def), allocatable, private :: stencil(:,:)
integer(kind=i_def), allocatable, private :: np_v(:,:)
integer(kind=i_def),              private :: np
real(kind=r_def),                 private :: x0
real(kind=r_def), allocatable,    private :: x_to_p(:)
real(kind=r_def), allocatable,    private :: coeff_matrix_v(:,:,:)

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_poly_flux_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
!       arg_type(GH_FIELD,   GH_READ,  W3, STENCIL(cross))             &
       arg_type(GH_FIELD,   GH_READ,  W3)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = EVALUATOR
contains
  procedure, nopass ::sample_poly_flux_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface sample_poly_flux_kernel_type
   module procedure sample_poly_flux_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public sample_poly_flux_code
public sample_poly_flux_init
contains

type(sample_poly_flux_kernel_type) function sample_poly_flux_kernel_constructor() result(self)
  return
end function sample_poly_flux_kernel_constructor

!> @brief Computes the fluxes for a tracer density
!! @param[in]  nlayers Number of layers
!! @param[out] flux Mass flux field to compute 
!! @param[in]  wind Wind field
!! @param[in]  density Tracer density
!! @param[in]  undf_w2 Number of unique degrees of freedom for the flux & wind fields
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at w2 nodes
!! @param[in]  undf_w3 Number of unique degrees of freedom for the density field
!! @param[in]  ndf_w3 Number of degrees of freedom per cell
!! @param[in]  stencil_size Size of the stencil (number of cells)
!! @param[in]  stencil_map Dofmaps for the stencil
subroutine sample_poly_flux_code( nlayers,              &
                                  flux,                 &
                                  wind,                 &
                                  density,              &
                                  ndf_w2,               &
                                  undf_w2,              &
                                  map_w2,               &
                                  basis_w2,             &
                                  ndf_w3,               &
                                  undf_w3,              &
                                  stencil_size,         &
                                  stencil_map           &
                                  )

  implicit none

  !Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_w2), intent(out)  :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)   :: density

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  integer(kind=i_def),                                 intent(in) :: stencil_size
  integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map

  !Internal variables
  integer(kind=i_def) :: k, df, p, i, nv
  real(kind=r_def)    :: direction, polynomial_density 
  real(kind=r_def)    :: z0

  ! Horizontal terms
  do k = 0, nlayers - 1
    do df = 1,4
      ! Check if this is the upwind cell
      direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                               out_face_normal(:,df))

      if ( direction >= 0.0_r_def ) then         
        do p = 1,np       
          density_stencil(p) = density( stencil_map(1,stencil(p,df)) + k )
        end do
        coeff(:) = matmul(coeff_matrix,density_stencil)
        polynomial_density = 0.0_r_def
        do p = 1,np
          polynomial_density = polynomial_density + coeff(p)*x_to_p(p)
        end do
        flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
      end if
    end do
  end do

  ! Vertical terms
  ! Bottom cell - apply zero flux on bottom boundary (df=5)
  !             - apply first order upwind on top boundary (df=6) 
  k = 0
  df = 5
  flux(map_w2(df)) = 0.0_r_def
  df = 6
  direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                           out_face_normal(:,df))
  if ( direction >= 0.0_r_def ) &
    flux(map_w2(df)) = wind(map_w2(df))*density(stencil_map(1,1))

  ! Internal cells
  do k = 1,nlayers - 2
    df = 5
    direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                             out_face_normal(:,df))
    if ( direction >= 0.0_r_def ) then
      nv = np_v(1,k)
      do p = 1,nv
        i = stencil_map(1,1) + k + (nv-1)/2 ! Location of first upwind point
        density_stencil(p) = density( i - (p-1)  )
      end do
    ! Use appropriate inverse matrix for the order of this point
      coeff(1:nv) = matmul(coeff_matrix_v(1:nv,1:nv,nv), &
                           density_stencil(1:nv))
      polynomial_density = 0.0_r_def
      z0 = real(1+nv/2,r_def)
      do p = 1,nv
        polynomial_density = polynomial_density + coeff(p)*z0**(p-1)
      end do
      flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density      
    end if

    df = 6
    direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                             out_face_normal(:,df))
    if ( direction >= 0.0_r_def ) then
      nv = np_v(2,k)
      do p = 1,nv
        i =  stencil_map(1,1) + k - (nv-1)/2 ! Location of first upwind point
        density_stencil(p) = density( i + (p-1) )
      end do
    ! Use appropriate inverse matrix for the order of this point
      coeff(1:nv) = matmul(coeff_matrix_v(1:nv,1:nv,nv), &
                           density_stencil(1:nv))
      polynomial_density = 0.0_r_def
      z0 = real(1+nv/2,r_def)
      do p = 1,nv
        polynomial_density = polynomial_density + coeff(p)*z0**(p-1)
      end do
      flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density      
    end if
  end do

  ! Top cell - apply zero flux on top boundary (df=6)
  !          - apply first order upwind on bottom boundary (df=5) 
  k = nlayers - 1
  df = 5
  direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                           out_face_normal(:,df))
  if ( direction >= 0.0_r_def ) &
    flux(map_w2(df)+k) = wind(map_w2(df)+k)*density(stencil_map(1,1)+k)
  df = 6
  flux(map_w2(df)+k) = 0.0_r_def

end subroutine sample_poly_flux_code

!=============================================================================!
!>@brief Initialise various quantities needed for sample_poly_flux_code
!>@param[in] order Polynomial order for flux computations
!>@param[in] nlayers Number of vertical layers
subroutine sample_poly_flux_init(order, nlayers)

  use matrix_invert_mod, only: matrix_invert

  implicit none
  
  integer(kind=i_def), intent(in) :: order, nlayers
  integer(kind=i_def)             :: p, i, j, k
  integer(kind=i_def)             :: nupwindcells, ndownwindcells
  real(kind=r_def), allocatable   :: inv_coeff_matrix(:,:)

  ! Order p polynomial has p+1 coefficients
  np = order + 1_i_def
 
  ! Compute the stencil table
  ! This give the cell id's in the stencilmap in a different order so that
  ! in each direction the first entry is the most upwind cell and the last
  ! the most downwind cell.
  ! e.g for a 1depth stencil map
  !     ---
  !     |5|
  !   -------
  !   |2|1|4|
  !   -------
  !     |3|
  !     ---
  !
  ! The stencil array will be:
  ! stencil = ( 4,1,2
  !             5,1,3
  !             2,1,4
  !             3,1,5) 
  ! that is ( i+1, i, i-1
  !           j+1, j, j-1
  !           i-1, i, i+1
  !           j-1, j, j+1) 
  allocate(stencil(np,4))

  ! Index of first upwind cell is  (j + (nupwindcells-1)*4)
  ! where j = [2,3,4,5] for [W,S,E,N] directions
  nupwindcells = int(real(np,r_def)/2.0_r_def)
  i = 1
  do p = 1,nupwindcells
    stencil(i,1) = int((4 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,2) = int((5 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,3) = int((2 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,4) = int((3 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    i = i + 1
  end do
  !index of centre cell in stencil is always 1
  stencil(i,:) = 1
  i = i + 1
  ! Index of first downwind cell is j
  ! where j = [5,4,3,2] for [W,S,E,N] directions
  ndownwindcells = int(real(np-1,r_def)/2.0_r_def)
  do p = 1,ndownwindcells
    stencil(i,1) = int(2 + (p-1)*4,i_def)
    stencil(i,2) = int(3 + (p-1)*4,i_def)
    stencil(i,3) = int(4 + (p-1)*4,i_def)
    stencil(i,4) = int(5 + (p-1)*4,i_def)
    i = i + 1
  end do

  ! Build Coefficient matrix of arbritray order
  allocate(coeff_matrix(np,np), inv_coeff_matrix(np,np))
  do i = 1,np
    do j = 1,np
      ! int(x^(j-1),x=i-1..i)
      inv_coeff_matrix(i,j) = (real(i,r_def)**j - real(i-1,r_def)**j)/real(j,r_def)
    end do
  end do
  call matrix_invert(inv_coeff_matrix,coeff_matrix,np)
  
  ! Find sampling point,
  ! first tracer point in stencil is at x = 0, 
  ! flux points are at boundaries so first flux point is x = 1/2
  ! x0 then depends upon the number of upwind cells used
  x0 = real(1+int(real(np)/2.0),r_def)

  ! Compute monomials at sampling point
  allocate( x_to_p(np) )
  do p = 1,np
    x_to_p(p) = x0**(p-1)
  end do

  allocate(coeff(np), density_stencil(np) )

  ! For vertical terms we may need to reduce orders near the boundaries
  ! due to lack of points so for each cell compute order and 
  ! number of points in stencil of vertical terms
  allocate( np_v(2,nlayers) )

  ! max polynomial order p = min(2#up -1,2#down)
  ! np = p + 1
  ! first index of array is for bottom boundary of cell 
  ! and second for the top boundary
  np_v(:,1) = 1_i_def
  do k = 1,nlayers-2  
    ! #down = k, #up = nl - k
    np_v(1,k+1) = int(min(np,min(2*(nlayers - k)-1, 2*k) + 1),i_def) 
    ! #down = nl - k - 1, #up = k + 1
    np_v(2,k+1) = int(min(np,min(2*(k+1) -1  , 2*(nlayers-1-k)) + 1),i_def)
  end do
  np_v(:,nlayers) = 1_i_def

  ! Compute inverses for vertical terms
  ! Each lower order matrix is just a subset of the high order matrix
  allocate(coeff_matrix_v(1:np,1:np,np) )
  do p = 1,np
    call matrix_invert(inv_coeff_matrix(1:p,1:p),coeff_matrix_v(1:p,1:p,p),p)
  end do

end subroutine sample_poly_flux_init

end module sample_poly_flux_kernel_mod
