!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Contains the routines used for Gaussian quadrature.

!> @details This module has a type for the Gaussian quadrature and a static
!> copy of the Gaussian quadrature that is used throughout the model. The first
!> time the Gaussian quadrature is required, it is created and a pointer to it
!> returned. Subsequent times, the pointer to the already created Guassian
!> quadrature is returned.

module gaussian_quadrature_mod
use constants_mod, only: r_def, pi, eps
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public :: gaussian_quadrature_type
  private
  !> allocatable arrays which holds the values of the gaussian quadrature
  real(kind=r_def), allocatable :: xgp(:), xgp_h(:,:), wgp(:), wgp_h(:)

  !> enumerated integer representing this instance of the gaussian_quadrature
  integer :: gq
contains
  !> Function returns a pointer to the Gaussian quadrature. If a Gaussian quadrature
  !> quadrature had not yet been created, it creates one before returning the pointer
  !> to it
  procedure, nopass :: get_instance
  !final     :: final_gauss
  !> Subroutine writes out an answer for a test
  !! @param self the calling gaussian quadrature
  procedure :: test_integrate

  !> Function quassian quadrature integration of a function f 
  !! @param self the calling gp type
  !! @param f real 3D array each of size ngp which holds the sample values of the
  !! function to be integrated
  !! @return real the value of the function thus integrated
  procedure :: integrate

  !> function returns the 2-d array of horizontal quadrature points
  procedure :: get_xgp_h
  
  !> function returns the 1-d array of vertical quadrature points
  procedure :: get_xgp_v

  !> function returns the enumerated integer for the gaussian_quadrature
  !! which is this gaussian_quadrature
  procedure :: which

end type

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
!> integer that defines the type of Gaussian quadrature required
integer, public, parameter      :: GQ3 = 203

!> integer The number of gaussian quadrature points in the vertical
integer, public, parameter      :: ngp_v = 3
!> integer The number of gaussian quadrature points in the horizontal
!! nqp_h=ngp_v*ngp_v for quads. They can be different (triangles or hexes)
!! but there is no setup code for this
integer, public, parameter      :: ngp_h = 9
!> All fields are integrated onto a fixed Guassian quadrature.
!> This is a static copy of that Gaussian quadrature object 
type(gaussian_quadrature_type), target, allocatable, save :: gq_3

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

function get_instance(gaussian_quadrature) result(instance)

  use log_mod, only : log_event, LOG_LEVEL_ERROR

  implicit none

  integer :: gaussian_quadrature
  type(gaussian_quadrature_type), pointer :: instance

  select case (gaussian_quadrature)
  case (GQ3)
    if(.not.allocated(gq_3)) then
      allocate(gq_3)
      call init_gauss(gq_3, GQ3) 
    end if
    instance => gq_3
  case default
    ! Not a recognised Gaussian quadrature. Logging an event with severity:
    ! LOG_LEVEL_ERROR will cause execution to abort
    call log_event( 'Gaussian quadrature type not recognised in '// &
                    'gaussian_quadrature%get_instance', LOG_LEVEL_ERROR )
  end select

  return
end function get_instance
 
subroutine init_gauss(self, gq)
  !-----------------------------------------------------------------------------
  ! Subroutine to compute the Gaussian points (xgp) and (wgp) wgphts 
  !-----------------------------------------------------------------------------
  implicit none

  class(gaussian_quadrature_type) :: self

  integer             :: i, j, m
  real(kind=r_def)    :: p1, p2, p3, pp, z, z1
  integer, intent(in) :: gq
  
  allocate( self%xgp(ngp_v) )
  allocate( self%wgp(ngp_v) ) 
  allocate( self%xgp_h(ngp_h,2) ) 
  allocate( self%wgp_h(ngp_h) ) 

  z1 = 0.0_r_def
  m = (ngp_v + 1) / 2

  !Roots are symmetric in the interval - so only need to find half of them  

  do i = 1, m ! Loop over the desired roots 

    z = cos( pi * (i - 0.25_r_def) / (ngp_v + 0.5_r_def) )

    !Starting with the above approximation to the ith root, we enter the main
    !loop of refinement by NEWTON'S method   
    do while ( abs(z-z1) > eps )
      p1 = 1.0_r_def
      p2 = 0.0_r_def

      !Loop up the recurrence relation to get the Legendre polynomial evaluated
      !at z                 
      do j = 1, ngp_v
        p3 = p2
        p2 = p1
        p1 = ((2.0_r_def * j - 1.0_r_def) * z * p2 - (j - 1.0_r_def) * p3) / j
      end do

      !p1 is now the desired Legendre polynomial. We next compute pp, its
      !derivative, by a standard relation involving also p2, the polynomial of one
      !lower order.      
      pp = ngp_v * (z * p1 - p2)/(z*z - 1.0_r_def)
      z1 = z
      z = z1 - p1/pp             ! Newton's Method  
    end do

    self%xgp(i) =  - z                                  ! Roots will be bewteen -1.0 & 1.0 
    self%xgp(ngp_v+1-i) =  + z                          ! and symmetric about the origin  
    self%wgp(i) = 2.0_r_def/((1.0_r_def - z*z) * pp*pp) ! Compute the wgpht and its       
    self%wgp(ngp_v+1-i) = self%wgp(i)                   ! symmetric counterpart         

  end do     ! i loop
      
  !Shift quad points from [-1,1] to [0,1]
  do i=1,ngp_v
    self%xgp(i) = 0.5_r_def*(self%xgp(i) + 1.0_r_def)
  end do

! This is correct for quads (will need modification for hexes/triangles)
  m = 1
  do i=1,ngp_v
    do j=1,ngp_v 
      self%xgp_h(m,1) = self%xgp(i)
      self%xgp_h(m,2) = self%xgp(j)
      self%wgp_h(m) = self%wgp(i)*self%wgp(j)
      
      m = m + 1
    end do
  end do

  self%gq = gq

  return
end subroutine init_gauss

subroutine test_integrate(self)
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------

  use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  class(gaussian_quadrature_type) :: self

  integer          :: i, k
  real(kind=r_def) :: func(ngp_v*ngp_v, ngp_v)
  real(kind=r_def) :: answer

  do i=1,ngp_h
    do k=1,ngp_v
      func(i,k) = self%xgp_h(i,1)*self%xgp_h(i,2)*1.0_r_def*1.0_r_def
    end do
  end do

  answer = self%integrate(func)
  write( log_scratch_space, '(A,F0.0)') 'int(x^2,x=0..1,y=0..1,z=0..1) = ', &
                                        answer
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  return
end subroutine test_integrate
  
!-----------------------------------------------------------------------------
! Compute 3D Gaussian integration of function f  
!-----------------------------------------------------------------------------  
!> Function to integrate a function f
!> @param[in] self the calling quadrature rule
!> @param[in] f the function to be integrated evaluated on the quadrature points
function integrate(self,f)
  implicit none

  class(gaussian_quadrature_type), intent(in) :: self

  real(kind=r_def), intent(in) :: f(ngp_h,ngp_v)
  real(kind=r_def)             :: integrate

  integer :: i,k

  integrate = 0.0_r_def
  do k=1,ngp_v 
    do i=1,ngp_h
      integrate = integrate + self%wgp_h(i)*self%wgp(k)*f(i,k)
    end do
  end do
  
  integrate = 0.125_r_def*integrate

  return
end function integrate

!-----------------------------------------------------------------------------
! Return Gaussian quadrature points
!-----------------------------------------------------------------------------
!> Function to return the quadrature points in the horizontal
!> @param[in] self the calling quadrature rule
!> @param[in] xgp_h the array to copy the quadrature points into
function get_xgp_h(self) result(xgp_h)
  implicit none
  class(gaussian_quadrature_type), intent(in) :: self
  real(kind=r_def) :: xgp_h(ngp_h,2)

  xgp_h(:,:) = self%xgp_h(:,:)
  return
end function get_xgp_h 

!> Function to return the quadrature points in the vertical
!> @param[in] self the calling quadrature rule
!> @param[in] xgp_v the array to copy the quadrature points into
function get_xgp_v(self) result(xgp_v)
  implicit none
  class(gaussian_quadrature_type), intent(in) :: self
  real(kind=r_def) :: xgp_v(ngp_v)

  xgp_v(:) = self%xgp(:)
  return
end function get_xgp_v

function which(self) result(gq)
  implicit none
  class(gaussian_quadrature_type),  intent(in) :: self
  integer :: gq
  
  gq = self%gq
  return
end function which

end module gaussian_quadrature_mod
