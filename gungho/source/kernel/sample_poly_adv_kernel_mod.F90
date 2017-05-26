!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the advective update u.grad(t) through fitting a high order 
!>        upwind 1D reconstruction
!> @details Compute the advective update for a tracer field using a high order
!>          polynomial fit to the tracer values. The stencil used for the
!>          polynomial is centred, with an upwind bias if an even number of 
!>          points are used, therefore for an upwind shceme an odd ordered 
!>          polynomial should be used. In the vertical the order is reduced 
!>          near the boundaries depending on the number of points available.
!>          This method is only valid for lowest order elements 

module sample_poly_adv_kernel_mod

use argument_mod,  only : arg_type, func_type,                  &
                          GH_FIELD, GH_WRITE, GH_READ,          &
                          W2, Wtheta, GH_BASIS, CELLS,          &
                          GH_EVALUATOR, EVALUATOR
use constants_mod, only : r_def, i_def
use kernel_mod,    only : kernel_type
use reference_element_mod, only: W, E, N, S


implicit none

! Precomputed operators, these are the same for all model columns
real(kind=r_def), allocatable,    private :: coeff_matrix(:,:)
real(kind=r_def), allocatable,    private :: coeff(:), tracer_stencil(:)
integer(kind=i_def), allocatable, private :: stencil(:,:)
integer(kind=i_def), allocatable, private :: np_v(:)
integer(kind=i_def),              private :: np
real(kind=r_def),                 private :: x0
real(kind=r_def), allocatable,    private :: dx0(:)
real(kind=r_def), allocatable,    private :: coeff_matrix_v(:,:,:)

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_poly_adv_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, Wtheta),                         &
       arg_type(GH_FIELD,   GH_READ,  Wtheta),                         &
       arg_type(GH_FIELD,   GH_READ,  W2)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = EVALUATOR
contains
  procedure, nopass ::sample_poly_adv_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface sample_poly_adv_kernel_type
   module procedure sample_poly_adv_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public sample_poly_adv_code
public sample_poly_adv_init
contains

type(sample_poly_adv_kernel_type) function sample_poly_adv_kernel_constructor() result(self)
  return
end function sample_poly_adv_kernel_constructor

!> @brief Computes the advective update
!! @param[in]  nlayers Number of layers
!! @param[out] advection Advective update to compute 
!! @param[in]  wind Wind field
!! @param[in]  tracer Tracer field
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the wind field
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at wt nodes
!! @param[in]  stencil_size Size of the stencil (number of cells)
!! @param[in]  stencil_map Dofmaps for the stencil
subroutine sample_poly_adv_code( nlayers,              &
                                 advection,            &
                                 tracer,               &
                                 wind,                 &
                                 ndf_wt,               &
                                 undf_wt,              &
                                 ndf_w2,               &
                                 undf_w2,              &
                                 map_w2,               &
                                 basis_w2,             &
                                 stencil_size,         &
                                 stencil_map           &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_wt), intent(out)  :: advection
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)   :: tracer

  real(kind=r_def), dimension(3,ndf_w2,ndf_wt), intent(in) :: basis_w2

  integer(kind=i_def),                                 intent(in) :: stencil_size
  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map

  ! Internal variables
  integer(kind=i_def) :: k, df, dft, p, dir, id
  real(kind=r_def)    :: u(3,nlayers+1)
  real(kind=r_def)    :: polynomial_tracer, advection_update, z0 

  ! Compute wind at tracer points
  u(:,:) = 0.0_r_def
  do k = 0, nlayers - 1
    do df = 1,ndf_w2
      do dft = 1,ndf_wt
        u(:,k + dft) = u(:,k + dft) + wind(map_w2(df) + k)*basis_w2(:,df,dft)
      end do
    end do   
  end do
  ! Average velocities across cell boundaries due to double counting
  ! the winds at the surface and lid are also reduced to mimic
  ! the weighting in the fem mass matrix
  ! which for the surface term is int(1-z,z=0..1) = 1/2
  ! whilst for interior terms it is  int(1-z,z=0..1) + int(z,z=0..1) = 1
  do k = 1,nlayers+1
    u(:,k) = u(:,k)*0.5_r_def   
  end do

  ! tracer data is stored contiguously and so can just use dft = 1
  ! and increment by 1 to get each subsequent level up to dft + nlayers 
  ! for the model lid
  do k = 0, nlayers
    dft = 1     
    ! Compute x stencil
    ! stencil_dir is East if u < 0 otherwise it is West
    if (  u(1,k+dft) >= 0.0_r_def ) then
      dir = W
    else
      dir = E
    end if
    do p = 1,np
      tracer_stencil(p) = tracer( stencil_map(dft,stencil(p,dir)) + k )
    end do
    coeff(:) = matmul(coeff_matrix,tracer_stencil)
    polynomial_tracer = 0.0_r_def
    do p = 1,np-1
      polynomial_tracer = polynomial_tracer + coeff(p+1)*dx0(p)
    end do
    ! Polynomial tracer contains the directional derivative
    ! ( dP/dx for u > 0, -dP/dx for u < 0) so need to cancel
    ! out sign of wind - hence presense of abs(u)
    advection_update = abs(u(1,k+dft))*polynomial_tracer

    ! Compute y stencil
    ! stencil_dir is north if v < 0 otherwise it is south
    if (  u(2,k+dft) >= 0.0_r_def ) then
      dir = S
    else
      dir = N
    end if
    do p = 1,np
      tracer_stencil(p) = tracer( stencil_map(dft,stencil(p,dir)) + k )
    end do
    coeff(:) = matmul(coeff_matrix,tracer_stencil)
    polynomial_tracer = 0.0_r_def
    do p = 1,np-1
      polynomial_tracer = polynomial_tracer + coeff(p+1)*dx0(p)
    end do
    ! Polynomial tracer contains the directional derivative
    ! ( dP/dy for v > 0, -dP/dy for v < 0) so need to cancel
    ! out sign of wind - hence presense of abs(v)
    advection_update = advection_update + abs(u(2,k+dft))*polynomial_tracer
   
    ! Compute z stencil, dir < 0 if w > 0
    if (  u(3,k+dft) >= 0.0_r_def ) then
      dir = -1
    else
      dir = 1
    end if
    do p = 1,np_v(k) 
      id = stencil_map(dft,1) + k + dir*(np_v(k)/2 - (p-1))
      tracer_stencil(p) = tracer( id )
    end do

    ! Use appropriate inverse matrix for the order of this point
    coeff(1:np_v(k)) = matmul(coeff_matrix_v(1:np_v(k),1:np_v(k),np_v(k)), &
                              tracer_stencil(1:np_v(k)))
    z0 = real(np_v(k)/2,r_def)
    polynomial_tracer = 0.0_r_def
    do p = 1,np_v(k)-1
      polynomial_tracer = polynomial_tracer + coeff(p+1)*real(p,r_def)*z0**(p-1)
    end do
    ! Polynomial tracer contains the directional derivative
    ! ( dP/dz for w > 0, -dP/dz for w < 0) so need to cancel
    ! out sign of wind - hence presense of abs(w)
    advection_update = advection_update + abs(u(3,k+dft))*polynomial_tracer
    
    advection(stencil_map(dft,1)+k) = advection_update
  end do
end subroutine sample_poly_adv_code

!=============================================================================!
!>@brief Initialise various quantities needed for sample_poly_adv_code
!>@param[in] order Polynomial order for advective computations
!>@param[in] nlayers Number of vertical layers
subroutine sample_poly_adv_init(order, nlayers)

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
  ! stencil = ( 2,1,4
  !             3,1,5
  !             4,1,2
  !             5,1,3)
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
    stencil(i,1) = int((2 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,2) = int((3 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,3) = int((4 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    stencil(i,4) = int((5 + (nupwindcells-1)*4) - (p-1)*4,i_def)
    i = i + 1
  end do
  !index of centre cell in stencil is always 1
  stencil(i,:) = 1
  i = i + 1
  ! Index of first downwind cell is j
  ! where j = [5,4,3,2] for [W,S,E,N] directions
  ndownwindcells = int(real(np-1,r_def)/2.0_r_def)
  do p = 1,ndownwindcells
    stencil(i,1) = int(4 + (p-1)*4,i_def)
    stencil(i,2) = int(5 + (p-1)*4,i_def)
    stencil(i,3) = int(2 + (p-1)*4,i_def)
    stencil(i,4) = int(3 + (p-1)*4,i_def)
    i = i + 1
  end do
  
  ! Build Coefficient matrix of arbritray order
  allocate(coeff_matrix(np,np), inv_coeff_matrix(np,np))
  do i = 1,np
    do j = 1,np
      ! (i-1)**(j-1)
      inv_coeff_matrix(i,j) = real(i-1,r_def)**(j-1)
    end do
  end do
  call matrix_invert(inv_coeff_matrix,coeff_matrix,np)

  ! Find sampling point,
  ! first tracer point in stencil is at x = 0, 
  ! x0 then depends upon the number of upwind cells used
  x0 = floor(real(np,r_def)/2.0_r_def)

  ! compute derivative of x0^p
  allocate( dx0(np-1) )
  do p = 1,np-1
    dx0(p) = real(p,r_def)*x0**(p-1)
  end do

  allocate(coeff(np), tracer_stencil(np) )
  
  ! For vertical terms we may need to reduce orders near the boundaries
  ! due to lack of points so for each cell compute order and 
  ! number of points in stencil of vertical terms
  allocate( np_v(0:nlayers) )

  np_v(0) = 1_i_def
  do k = 1,nlayers-1
    np_v(k) = int(min(np,min(2*k+1,2*(nlayers-k)+1)),i_def)
  end do
  np_v(nlayers) = 1_i_def
  
! Compute inverses for vertical
! Each lower order matrix is just a subset of the high order matrix
  allocate(coeff_matrix_v(1:np,1:np,np) )
  do p = 1,np
    call matrix_invert(inv_coeff_matrix(1:p,1:p),coeff_matrix_v(1:p,1:p,p),p)
  end do

end subroutine sample_poly_adv_init

end module sample_poly_adv_kernel_mod
