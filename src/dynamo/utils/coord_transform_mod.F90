!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2012.
! However, it has been created by John Thuburn.
!-------------------------------------------------------------------------------
!>  @brief   Routines for coordinate transformations.
!!
!!  @details Contains routines for conversion of e.g. lat-long to Cartesian XYZ.
!-------------------------------------------------------------------------------
module coord_transform_mod

use constants_mod, only : r_def

implicit none

private

! Public subroutines
public :: ll2xyz
public :: llr2xyz
public :: xyz2ll
public :: xyz2llr
public :: starea2
public :: spdist
public :: cartesian_distance
public :: sphere2cart_vector
public :: cart2sphere_vector

!--------------------------------------------------------------------------------
! Contained functions / subroutines
!--------------------------------------------------------------------------------
contains

!--------------------------------------------------------------------------------
!>  @brief  Converts longitude and latitude to Cartesian coordinates on the unit
!!          sphere.
!!
!!  @param[in]   long  The longitude to convert.
!!  @param[in]   lat   The latitude  to convert.
!!  @param[out]  x     The x coordinate.
!!  @param[out]  y     The y coordinate.
!!  @param[out]  z     The z coordinate.
!--------------------------------------------------------------------------------
subroutine ll2xyz(long,lat,x,y,z)
  implicit none
 
  !Arguments
  real(kind=r_def), intent(in)  :: long,lat
  real(kind=r_def), intent(out) :: x,y,z 
 
  !Internal variables
  real(kind=r_def) :: cos_long, sin_long, cos_lat, sin_lat
 
  sin_long = sin(long)
  cos_long = cos(long)
  sin_lat  = sin(lat)
  cos_lat  = cos(lat)
 
  x = cos_long * cos_lat
  y = sin_long * cos_lat
  z = sin_lat
 
  return
end subroutine ll2xyz

!--------------------------------------------------------------------------------
!>  @brief  Converts longitude and latitude to Cartesian coordinates on 
!!          a sphere with some specified radius.
!!
!!  @param[in]   long    The longitude to convert.
!!  @param[in]   lat     The latitude  to convert.
!!  @param[in]   radius  The radius of the sphere.
!!  @param[out]  x       The x coordinate.
!!  @param[out]  y       The y coordinate.
!!  @param[out]  z       The z coordinate.
!--------------------------------------------------------------------------------
subroutine llr2xyz(long,lat,radius,x,y,z)
  implicit none
 
  !Arguments
  real(kind=r_def), intent(in)  :: long,lat,radius
  real(kind=r_def), intent(out) :: x,y,z 
 
  !Internal variables
  real(kind=r_def) :: cos_long, sin_long, cos_lat, sin_lat
 
  sin_long = sin(long)
  cos_long = cos(long)
  sin_lat  = sin(lat)
  cos_lat  = cos(lat)
 
  x = radius * cos_long * cos_lat
  y = radius * sin_long * cos_lat
  z = radius * sin_lat
 
  return
end subroutine llr2xyz

!--------------------------------------------------------------------------------
!>  @brief  Converts Cartesian coordinates to longitude and latitude.
!!
!!  @param[in]   x     The x coordinate to convert.
!!  @param[in]   y     The y coordinate to convert.
!!  @param[in]   z     The z coordinate to convert.
!!  @param[out]  long  The longitude.
!!  @param[out]  lat   The latitude.
!--------------------------------------------------------------------------------
subroutine xyz2ll(x,y,z,long,lat)
  use constants_mod, only : PI
  implicit none

  !Arguments
  real(kind=r_def), intent(in)  :: x,y,z
  real(kind=r_def), intent(out) :: long, lat

  !Internal variables
  real(kind=r_def) :: tan_long, tan_lat, radius

  if (x == 0.0_r_def) then
    if (y >= 0.0_r_def) then
      long = 0.5_r_def*PI
    else
      long = 1.5_r_def*PI
    end if
  else
    tan_long=y/x
    long=atan(tan_long)
    if (x < 0.0_r_def) then
      long=long+PI
    endif
    if (long < 0.0_r_def) then
      long=long+2.0_r_def*PI
    endif
  end if

  radius = sqrt(x*x+y*y)
  if (radius == 0.0_r_def) then
    if (z > 0.0_r_def) then
      lat=0.5_r_def*pi
    else
      lat=-0.5_r_def*pi
    end if
  else
    tan_lat = z / radius
    lat = atan(tan_lat)
  end if

  return
end subroutine xyz2ll

!-------------------------------------------------------------------------------
!>  @brief  Converts Cartesian coordinates to longitude, latitude and radius
!!  @param[in]   x     The x coordinate to convert.
!!  @param[in]   y     The y coordinate to convert.
!!  @param[in]   z     The z coordinate to convert.
!!  @param[out]  long  The longitude.
!!  @param[out]  lat   The latitude.
!!  @param[out]  radius  The radius of the sphere..
!-------------------------------------------------------------------------------     
subroutine xyz2llr(x,y,z,long,lat,r)
  use constants_mod, only: PI
  implicit none

  real(kind=r_def), intent(in)  :: x, y, z
  real(kind=r_def), intent(out) :: long, lat, r
  real(kind=r_def)              :: tln, tlt
  real(kind=r_def)              :: tol = 10.0e-8_r_def

  if ( abs(x) < tol ) then
    if ( y >= 0.0_r_def ) then
      long = 0.5_r_def*PI
    else
      long = 1.5_r_def*PI
    end if
  else
    tln = y/x
    long = atan(tln)
    if ( x < 0.0_r_def ) then
      long = long + PI
    end if
    if ( long < 0.0_r_def ) then
      long = long + 2.0_r_def*PI
    end if
  end if

  r = sqrt(x*x + y*y)
  if ( abs(r) < tol ) then
    if (z > 0.0_r_def ) then
      lat =  0.5_r_def*PI
    else
      lat = -0.5_r_def*PI
    end if
  else
    tlt = z/r
    lat = atan(tlt)
  end if
  
  r = sqrt(x*x + y*y + z*z)  

end subroutine xyz2llr


!-------------------------------------------------------------------------------
!>  @brief  Calculates the area of a spherical triangle.
!!
!!  @details Calculates the area of the spherical triangle whose corners have
!!           Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2).
!!           The formula below is more robust to roundoff error than the better
!!           known sum of angle - PI formula
!!
!!  @param[in]   x0    Coordinate of a triangle corner.
!!  @param[in]   y0    Coordinate of a triangle corner.
!!  @param[in]   z0    Coordinate of a triangle corner.
!!  @param[in]   x1    Coordinate of a triangle corner.
!!  @param[in]   y1    Coordinate of a triangle corner.
!!  @param[in]   z1    Coordinate of a triangle corner.
!!  @param[in]   x2    Coordinate of a triangle corner.
!!  @param[in]   y2    Coordinate of a triangle corner.
!!  @param[in]   z2    Coordinate of a triangle corner.
!!  @param[out]  area  Area of the spherical triangle.
!-------------------------------------------------------------------------------
subroutine starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)
  implicit none

  !Arguments
  real(kind=r_def), intent(in)  :: x0, y0, z0   
  real(kind=r_def), intent(in)  :: x1, y1, z1
  real(kind=r_def), intent(in)  :: x2, y2, z2
  real(kind=r_def), intent(out) :: area

  !Internal variables
  real(kind=r_def) :: d0,d1,d2,s,t0,t1,t2,t3

  !Distances between pairs of points
  call spdist(x0,y0,z0,x1,y1,z1,d2)
  call spdist(x1,y1,z1,x2,y2,z2,d0)
  call spdist(x2,y2,z2,x0,y0,z0,d1)

  !Half perimeter
  s=0.5_r_def*(d0+d1+d2)

  !Tangents
  t0 = tan(0.5_r_def*(s-d0))
  t1 = tan(0.5_r_def*(s-d1))
  t2 = tan(0.5_r_def*(s-d2))
  t3 = tan(0.5_r_def*s)

  !Area
  area = 4.0_r_def*atan(sqrt(t0*t1*t2*t3))

  return
end subroutine starea2

!-------------------------------------------------------------------------------
!> @brief  Calculates the spherical distance between two points.
!!
!! @details  Calculates the spherical distance s between two points with 
!!           Cartesian coordinates (x1,y1,z1), (x2,y2,z2) on the unit sphere.
!!
!! @param[in]  x1  First Cartesian coordinate.
!! @param[in]  x2  First Cartesian coordinate.
!! @param[in]  x3  First Cartesian coordinate.
!! @param[in]  y1  Second Cartesian coordinate.
!! @param[in]  y2  Second Cartesian coordinate.
!! @param[in]  y3  Second Cartesian coordinate.
!! @param[out] s   Spherical distance between the points.
!-------------------------------------------------------------------------------
subroutine spdist(x1,y1,z1,x2,y2,z2,s)
  implicit none

  !Arguments
  real(kind=r_def), intent(in)  :: x1, y1, z1, x2, y2, z2
  real(kind=r_def), intent(out) :: s

  !Internal variables
  real(kind=r_def) :: dx, dy, dz
  real(kind=r_def) :: ad

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  ad = sqrt(dx*dx + dy*dy + dz*dz)
  s = 2.0_r_def*asin(0.5_r_def*ad)

  return
end subroutine spdist

!-------------------------------------------------------------------------------
!> @brief  Calculates the Cartesian distance between two points.
!!
!! @details  Calculates the Cartesian distance s between two points with 
!!           Cartesian coordinates (x1,y1,z1), (x2,y2,z2) 
!!
!! @param[in]  x(3)  First point in Cartesian coordinates.
!! @param[in]  y(3)  Second point in Cartesian coordinates.
!! @param[out] s   Cartesian distance between the points.
!-------------------------------------------------------------------------------
pure function cartesian_distance(x,y) result( s )
  implicit none

  !Arguments
  real(kind=r_def), intent(in)  :: x(3), y(3)
  real(kind=r_def)              :: s

  !Internal variables
  real(kind=r_def) :: dx, dy, dz

  dx = y(1) - x(1)
  dy = y(2) - x(2)
  dz = y(3) - x(3)

  s = sqrt(dx*dx + dy*dy + dz*dz)
  return
end function cartesian_distance


!-------------------------------------------------------------------------------
!> @brief Converts a vector in sperical coordinates to one in
!> Cartesian coodinates
!> @param[in]  llr location in spherical coodinates
!> @param[in]  dlambda input vector in spherical coordinates
!> @param[return] dx output vector in Cartesian coordinates
!-------------------------------------------------------------------------------
function sphere2cart_vector( dlambda, llr ) result ( dx )
  use constants_mod,     only: r_def
  use matrix_invert_mod, only: matrix_invert_3x3
  implicit none

  real(kind=r_def), intent(in)  :: dlambda(3)
  real(kind=r_def), intent(in)  :: llr(3)
  real(kind=r_def)              :: dx(3)

  real(kind=r_def)              :: A(3,3), A_inv(3,3)

! Form transformation matrix
  A(1,:) = (/ -sin(llr(1)),              cos(llr(1)),              0.0_r_def   /)
  A(2,:) = (/ -sin(llr(2))*cos(llr(1)), -sin(llr(2))*sin(llr(1)),  cos(llr(2)) /)
  A(3,:) = (/  cos(llr(2))*cos(llr(1)),  cos(llr(2))*sin(llr(1)),  sin(llr(2)) /)
! form inverse
  A_inv = matrix_invert_3x3(A)
  dx(:) = matmul(A_inv, dlambda)
  return

end function sphere2cart_vector

!-------------------------------------------------------------------------------
!> @brief Converts a flux vector in Cartesian coordinates (x,y,z) 
!> to one in spherical coodinates (lambda,phi,r)
!> @detail Converts 3d Cartesian velocity u = (u1,u2,u3) in Cartesian coordinates at (x,y,z) 
!> to spherical velocity v = (v1,v2,v3) (in m/s) in spherical coordinates at (lambda,phi,r)
!> @param[in]  x_vec vector (x,y,z) location in Cartesian coodinates
!> @param[in]  Cartesian_vec components of a flux vector (u,v,w) in Cartesian coordinates
!> @param[out] spherical_vec components of a flux vector (u,v,w) on spherical coordinates
!-------------------------------------------------------------------------------
pure function cart2sphere_vector(x_vec, cartesian_vec) result ( spherical_vec )
     use constants_mod, only: r_def, PI

     implicit none

     real(kind=r_def), intent(in)  :: x_vec(3)
     real(kind=r_def), intent(in)  :: cartesian_vec(3)

     real(kind=r_def)              :: spherical_vec(3)

     real(kind=r_def)              :: r, t, phi

     t = x_vec(1)**2 + x_vec(2)**2
     r = sqrt(t + x_vec(3)**2)

     spherical_vec(1) = (- x_vec(2)*cartesian_vec(1) &
                         + x_vec(1)*cartesian_vec(2) ) / t
     spherical_vec(2) = (-x_vec(1)*x_vec(3)*cartesian_vec(1) &
                         -x_vec(2)*x_vec(3)*cartesian_vec(2) &
                                        + t*cartesian_vec(3))/(r*r*sqrt(t))
     spherical_vec(3) = (x_vec(1)*cartesian_vec(1) &
                       + x_vec(2)*cartesian_vec(2) &
                       + x_vec(3)*cartesian_vec(3)) / r

! convert from (dlambda/dt,dphi/dt,dr/dt) to (u,v,w) in m/s
     phi = 0.5_r_def*PI - acos(x_vec(3)/r)
     spherical_vec(1) = spherical_vec(1)*r*cos(phi)
     spherical_vec(2) = spherical_vec(2)*r

end function cart2sphere_vector

end module coord_transform_mod

