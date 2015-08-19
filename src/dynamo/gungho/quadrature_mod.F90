!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Contains the routines used for (Gaussian) quadrature.

!> @details This module has a type for the (Gaussian) quadrature and a static
!> copy of the quadrature that is used throughout the model. The first
!> time the quadrature is required, it is created and a pointer to it
!> returned. Subsequent times, the pointer to the already created 
!> quadrature is returned.

module quadrature_mod
use constants_mod, only: r_def, PI, EPS
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public :: quadrature_type
  private
  !> allocatable arrays which holds the values of the gaussian quadrature
  real(kind=r_def), allocatable :: xqp(:), xqp_h(:,:), wqp(:), wqp_h(:)

  !> enumerated integer representing this instance of the quadrature rule
  integer :: qr

  !> integer Number of quadrature points in the horizontal
  integer :: nqp_h

    !> integer Number of quadrature points in the vertical
  integer :: nqp_v

contains
  !> @brief Returns a pointer to the quadrature. If a quadrature had not
  !! yet been created, it creates one before returning the pointer to it.
  !> @param[in] quadrature Integer, quadrature rule
  !> @param[in] nqp_h Integer, number of quadrature points in the horizontal
  !> @param[in] nqp_v Integer, number of quadrature points in the vertical
  procedure, nopass :: get_instance

  !final     :: final_gauss
  !> @brief Writes out an answer for a test
  !> @param self The calling gaussian quadrature
  procedure :: test_integrate

  !> @brief Computes 3D quadrature integration of function f
  !> @param[in] self The calling quadrature rule
  !> @param f Real 3D array each of size ngp which holds the sample values of the
  !! function to be integrated on the quadrature points
  !> @return real The value of the function thus integrated
  procedure :: integrate

  !> @brief Returns the 2-d array of quadrature points in the horizontal
  !> @param[in] self The calling quadrature rule
  !> @return xgp_h The array to copy the quadrature points into
  procedure :: get_xqp_h

  !> @brief Returns the 1-d array of vertical quadrature points
  !> @param[in] self The calling quadrature rule
  !> @return xqp_v The array to copy the quadrature points into
  procedure :: get_xqp_v

  !> @brief Returns the enumerated integer for the gaussian_quadrature
  !! which is this gaussian_quadrature
  !> @param[in] self The calling quadrature rule
  !> @return qr The quadrature rule
  procedure :: which

  !> @brief Returns the 1-d array of horizontal quadrature weights
  !> @param[in] self The calling quadrature rule
  !> @return wqp_h The pointer to the horizontal quadrature weights
  procedure :: get_wqp_h

  !> @brief Returns the 1-d array of vertical quadrature weights
  !> @param[in] self The calling quadrature rule
  !> @return wqp_v The pointer to the vertical quadrature weights
  procedure :: get_wqp_v

  !> @brief Returns the number of quadrature points in the vertical
  !> @param[in] self The calling quadrature rule
  !> @return nqp_v Number of quadrature points in the vertical
  procedure :: get_nqp_v

  !> @brief Returns the number of quadrature points in the horizontal
  !> @param[in] self The calling quadrature rule
  !> @return nqp_h Number of quadrature points in the horizontal
  procedure :: get_nqp_h

end type

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
!> Integer that defines the type of Gaussian quadrature required
integer, public, parameter      :: QR3 = 1001, QRnewton = 1002

!> All fields are integrated onto a fixed Guassian quadrature.
!> This is a static copy of that Gaussian quadrature object 
type(quadrature_type), target, allocatable, save :: qr_3, qr_newton

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Returns a pointer to the quadrature
!-------------------------------------------------------------------------------
function get_instance(quadrature,nqp_h,nqp_v) result(instance)

  use log_mod, only : log_event, LOG_LEVEL_ERROR

  implicit none

  integer :: quadrature
  type(quadrature_type), pointer :: instance
  integer, intent(in) :: nqp_h, nqp_v

  select case (quadrature)
  case (QR3)
    if(.not.allocated(qr_3)) then
      allocate(qr_3)
      call init_quadrature(qr_3, QR3, nqp_h, nqp_v)
      call create_gaussian_quadrature(qr_3)
    end if
    instance => qr_3
  case (QRnewton)
    if(.not.allocated(qr_newton)) then
      allocate(qr_newton)
      call init_quadrature(qr_newton, QRnewton, nqp_h, nqp_v)
      call create_newton_cotes_quadrature(qr_newton)
    end if
    instance => qr_newton
  case default
    instance => null() ! To keep GCC's uninitialised use checker happy.
    ! Not a recognised  quadrature. Logging an event with severity:
    ! LOG_LEVEL_ERROR will cause execution to abort
    call log_event( 'Quadrature type not recognised in '// &
                    'quadrature%get_instance', LOG_LEVEL_ERROR )
  end select

  return
end function get_instance

!-------------------------------------------------------------------------------
! Initialises the quadrature rule
!-------------------------------------------------------------------------------
!> @brief Initialises the quadrature rule 
!> @param[in] nqp_h Integer, number of quadrature points in the horizontal
!> @param[in] nqp_v Integer, number of quadrature points in the vertical
!> @param[in] qr Integer, quadrature rule
!-------------------------------------------------------------------------------
subroutine init_quadrature(self, qr, nqp_h, nqp_v)

  implicit none

  class(quadrature_type) :: self
  integer,intent(in)  :: nqp_h,nqp_v
  integer, intent(in) :: qr

  self%nqp_h = nqp_h
  self%nqp_v = nqp_v

  allocate( self%xqp(nqp_v) )
  allocate( self%wqp(nqp_v) ) 
  allocate( self%xqp_h(nqp_h,2) ) 
  allocate( self%wqp_h(nqp_h) )

  self%xqp(:) = 0.0_r_def
  self%wqp(:) = 0.0_r_def
  self%xqp_h(:,:) = 0.0_r_def
  self%wqp_h(:) = 0.0_r_def
  self%qr = qr
end subroutine init_quadrature 

!-----------------------------------------------------------------------------
! Computes Gaussian quadrature points and weights 
!-----------------------------------------------------------------------------
!> @brief Computes Gaussian quadrature points (xqp) and weights (wqp)
!> @param[in] self The calling quadrature rule
!-----------------------------------------------------------------------------
subroutine create_gaussian_quadrature(self)

  implicit none

  class(quadrature_type) :: self 
  integer             :: i, j, m
  real(kind=r_def)    :: p1, p2, p3, pp, z, z1
  real(kind=r_def), parameter :: DOMAIN_CHANGE_FACTOR = 0.5_r_def
  real(kind=r_def), parameter :: DOMAIN_SHIFT_FACTOR  = 1.0_r_def

  z1 = 0.0_r_def
  m = (self%nqp_v + 1) / 2

  !Roots are symmetric in the interval - so only need to find half of them

  do i = 1, m ! Loop over the desired roots

    z = cos( PI * (i - 0.25_r_def) / (self%nqp_v + 0.5_r_def) )

    !Starting with the above approximation to the ith root, we enter the main
    !loop of refinement by NEWTON'S method
    pp = 1.0_r_def
    do while ( abs(z-z1) > eps )
      p1 = 1.0_r_def
      p2 = 0.0_r_def

      !Loop up the recurrence relation to get the Legendre polynomial evaluated
      !at z                 
      do j = 1, self%nqp_v
        p3 = p2
        p2 = p1
        p1 = ((2.0_r_def * j - 1.0_r_def) * z * p2 - (j - 1.0_r_def) * p3) / j
      end do

      !p1 is now the desired Legendre polynomial. We next compute pp, its
      !derivative, by a standard relation involving also p2, the polynomial of
      ! one lower order.
      pp = self%nqp_v * (z * p1 - p2)/(z*z - 1.0_r_def)
      z1 = z
      z = z1 - p1/pp             ! Newton's Method  
    end do

    self%xqp(i) =  - z                                  ! Roots will be bewteen -1.0 & 1.0 
    self%xqp(self%nqp_v+1-i) =  + z                     ! and symmetric about the origin  
    self%wqp(i) = 2.0_r_def/((1.0_r_def - z*z) * pp*pp) ! Compute the wgpht and its       
    self%wqp(self%nqp_v+1-i) = self%wqp(i)              ! symmetric counterpart         

  end do ! i loop

  !Shift quad points from [-1,1] to [0,1]
  do i=1,self%nqp_v
    self%xqp(i) = DOMAIN_CHANGE_FACTOR*(self%xqp(i) + DOMAIN_SHIFT_FACTOR)
    self%wqp(i) = DOMAIN_CHANGE_FACTOR*self%wqp(i)
  end do

  ! This is correct for quads (will need modification for hexes/triangles)
  m = 1
  do i=1,self%nqp_v
    do j=1,self%nqp_v 
      self%xqp_h(m,1) = self%xqp(i)
      self%xqp_h(m,2) = self%xqp(j)
      self%wqp_h(m) = self%wqp(i)*self%wqp(j)

      m = m + 1
    end do
  end do

  return
end subroutine create_gaussian_quadrature

!-----------------------------------------------------------------------------
! Computes Newton Cotes quadrature points and weights
!-----------------------------------------------------------------------------
!> @brief Computes Newton Cotes quadrature points (xqp) and weights (wqp)
!> @param[in] self The calling quadrature rule
!-----------------------------------------------------------------------------
subroutine create_newton_cotes_quadrature(self)

  use matrix_invert_mod, only: matrix_invert
 
  implicit none

  class(quadrature_type) :: self 
  integer :: i, j, ij
  real(kind=r_def), allocatable :: A(:,:), Ainv(:,:), b(:)

  allocate( A   (self%nqp_v,self%nqp_v),  &
            Ainv(self%nqp_v,self%nqp_v),  &
            b   (self%nqp_v) )

  if ( self%nqp_v == 1 ) then
    self%xqp(1) = 0.5_r_def
  else
    do i = 1,self%nqp_v
     self%xqp(i) = real(i-1)/real(self%nqp_v-1)
    end do
  end if
  ij = 1
  do i = 1,self%nqp_v
    do j = 1,self%nqp_v
      self%xqp_h(ij,1) =  self%xqp(i)
      self%xqp_h(ij,2) =  self%xqp(j)
      ij = ij + 1
    end do    
  end do

  ! Compute weights
  do i = 1,self%nqp_v
    b(i) = (real(self%nqp_v)**i-1.0_r_def)/real(i)
  end do
  ! Compute coefficient matrix A
  do i = 1,self%nqp_v
    do j = 1,self%nqp_v
      A(i,j) = real(j)**(i-1)
    end do
  end do
  call matrix_invert(A,Ainv,self%nqp_v)
  self%wqp(:) = matmul(Ainv,b)
  ij = 1
  do i = 1,self%nqp_v
    do j = 1,self%nqp_v
      self%wqp(ij) = self%wqp(i)*self%wqp(j)
    end do
  end do
  
end subroutine create_newton_cotes_quadrature

!-----------------------------------------------------------------------------
! Writes out an answer for a test
!-----------------------------------------------------------------------------
subroutine test_integrate(self)

  use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  class(quadrature_type) :: self

  integer          :: i, k
  real(kind=r_def) :: func(self%nqp_h, self%nqp_v)
  real(kind=r_def) :: answer

  do i=1,self%nqp_h
    do k=1,self%nqp_v
      func(i,k) = self%xqp_h(i,1)*self%xqp_h(i,2)*1.0_r_def*1.0_r_def
    end do
  end do

  answer = self%integrate(func)
  write( log_scratch_space, '(A,F0.0)') 'int(x^2,x=0..1,y=0..1,z=0..1) = ', &
                                        answer
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  return
end subroutine test_integrate
  
!-----------------------------------------------------------------------------
! Computes 3D quadrature integration of function f  
!-----------------------------------------------------------------------------  
function integrate(self,f)
  implicit none

  class(quadrature_type), intent(in) :: self

  real(kind=r_def), intent(in) :: f(self%nqp_h,self%nqp_v)
  real(kind=r_def)             :: integrate

  integer :: i,k

  integrate = 0.0_r_def
  do k=1,self%nqp_v 
    do i=1,self%nqp_h
      integrate = integrate + self%wqp_h(i)*self%wqp(k)*f(i,k)
    end do
  end do
  
  return
end function integrate

!-----------------------------------------------------------------------------
! Returns the quadrature points in the horizontal
!-----------------------------------------------------------------------------
function get_xqp_h(self) result(xqp_h)
  implicit none
  class(quadrature_type), target, intent(in) :: self
  real(kind=r_def), pointer :: xqp_h(:,:)

  xqp_h => self%xqp_h
  return
end function get_xqp_h

!-----------------------------------------------------------------------------
! Returns the quadrature points in the vertical
!-----------------------------------------------------------------------------
function get_xqp_v(self) result(xqp_v)
  implicit none
  class(quadrature_type), target, intent(in) :: self
  real(kind=r_def), pointer :: xqp_v(:)

  xqp_v => self%xqp
  return
end function get_xqp_v

!-----------------------------------------------------------------------------
! Returns the enumerated integer for quadrature rule
!-----------------------------------------------------------------------------
function which(self) result(qr)
  implicit none
  class(quadrature_type),  intent(in) :: self
  integer :: qr

  qr = self%qr
  return
end function which

!-----------------------------------------------------------------------------
! Returns the number of quadrature points in the vertical
!-----------------------------------------------------------------------------
function get_nqp_v(self) result(nqp_v)
  implicit none
  class(quadrature_type), intent(in) :: self
  integer :: nqp_v

  nqp_v = self%nqp_v
  return
end function get_nqp_v

!-----------------------------------------------------------------------------
! Returns the number of quadrature points in the horizontal
!-----------------------------------------------------------------------------
function get_nqp_h(self) result(nqp_h)
  implicit none
  class(quadrature_type), intent(in) :: self
  integer :: nqp_h

  nqp_h = self%nqp_h
  return
end function get_nqp_h

!-----------------------------------------------------------------------------
! Returns the quadrature weights in the horizontal
!-----------------------------------------------------------------------------
function get_wqp_h(self) result(wqp_h)
  implicit none
  class(quadrature_type), target, intent(in) :: self
  real(kind=r_def), pointer :: wqp_h(:) 

  wqp_h => self%wqp_h
  return
end function get_wqp_h 

!-----------------------------------------------------------------------------
! Returns the quadrature weights in the vertical
!-----------------------------------------------------------------------------
function get_wqp_v(self) result(wqp_v)
  implicit none
  class(quadrature_type), target, intent(in) :: self
  real(kind=r_def), pointer :: wqp_v(:) 

  wqp_v => self%wqp
  return
end function get_wqp_v 


end module quadrature_mod
