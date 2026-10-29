!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Module contains utility routines for mesh rotation
!> @details The rotations are based on the assumption that the coords in the
!>          unrotated frame of reference are rotated such that the TRUE NORTH
!>          moves to the target_north_pole.
module rotation_mod

  use constants_mod,       only: r_def, i_def, str_def, PI
  use log_mod,             only: log_event, LOG_LEVEL_ERROR
  use coord_transform_mod, only: rebase_longitude_range

  implicit none

  private

  ! Set Parameters for true north pole / null island
  real(r_def), parameter, public :: TRUE_NORTH_POLE_XYZ(3) = &
                                    (/ 0.0_r_def,  0.0_r_def, 1.0_r_def /)
  real(r_def), parameter, public :: TRUE_NORTH_POLE_LL(2)  = &
                                    (/ 0.0_r_def, 90.0_r_def /)
  real(r_def), parameter, public :: TRUE_NULL_ISLAND_LL(2) = &
                                    (/ 0.0_r_def,  0.0_r_def /)

  public :: rotate_mesh_coords
  public :: get_target_north_pole
  public :: get_target_null_island

  contains
  !-------------------------------------------------------------------------------
  !> @brief   Rotates mesh coordinates
  !> @details Rotates the coordinates to give a new pole.
  !>
  !> The rotation is achieved by rotating about the vector rot_vec and
  !> by an angle alpha_vec.
  !>
  !> @param[in,out] vert_coords       Array of (lon,lat) coordinates to be rotated
  !> @param[in]     target_north_pole The target north pole (lon,lat) to rotate to.
  !-------------------------------------------------------------------------------
  subroutine rotate_mesh_coords(vert_coords, target_north_pole)

    use coord_transform_mod, only: xyz2ll, ll2xyz, rodrigues_rotation
    use cross_product_mod,   only: cross_product

    implicit none

    real(r_def), intent(inout) :: vert_coords(:,:)
    real(r_def), intent(in)    :: target_north_pole(2)

    real(r_def)    :: x_vec(3)     ! Cartesian vector of points
    real(r_def)    :: rot_vec(3)   ! Cartesian axis of rotation
    real(r_def)    :: new_pole(3)  ! Cartesian vector for new pole
    real(r_def)    :: alpha_rot    ! Angle of rotation about rot_vec
    integer(i_def) :: icell, nverts
    real(r_def)    :: rotation_angle

    logical        :: rotate_the_pole    ! To determine if we want a new north
    logical        :: rotate_about_pole  ! To determine if we want to rotate about north

    ! Set up axis of rotation.
    call ll2xyz(target_north_pole(1), &
                target_north_pole(2), &
                new_pole(1),          &
                new_pole(2),          &
                new_pole(3))

    ! Rotate from true north to the new pole.
    rot_vec   = cross_product(true_north_pole_xyz, new_pole)

    alpha_rot = acos(dot_product(true_north_pole_xyz, new_pole) /   &
       sqrt(dot_product(new_pole,   new_pole) *                     &
       dot_product(true_north_pole_xyz, true_north_pole_xyz)))

    ! If angle is less than ~ 0.0001 degrees (i.e. norm2(rot_vec) < 1.e-6)
    ! then don't change the axis.
    rotate_the_pole = .true.
    if (norm2(rot_vec) < 1.e-6_r_def) rotate_the_pole = .false.

    ! The rotation_angle is chosen to from 0 degrees to pole_lon.
    rotation_angle = target_north_pole(1) + PI

    ! If rotation_angle < 0.0001 degrees, then don't rotate.
    rotate_about_pole = .true.
    if (abs(target_north_pole(1)) < 0.0001_r_def) rotate_about_pole = .false.

    if (rotate_the_pole .or. rotate_about_pole) then

      nverts = size( vert_coords, dim=2 )

      do icell = 1, nverts

        call ll2xyz(vert_coords(1,icell),         &
                    vert_coords(2,icell),         &
                    x_vec(1), x_vec(2), x_vec(3))

        ! First rotate pole to a new longitude - rotate around true_north_pole_xyz by an
        ! rotation_angle.
        if (rotate_about_pole) x_vec = rodrigues_rotation(x_vec, true_north_pole_xyz,  &
           rotation_angle)

        ! Next rotate pole to a new latitude - keep the pole longitude fixed -
        ! rotate about rot_vec by the angle alpha_rot.
        if (rotate_the_pole) x_vec = rodrigues_rotation(x_vec, rot_vec,         &
           alpha_rot)

        ! Convert back to lat, lon and send back to vert_coords.
        call xyz2ll(x_vec(1), x_vec(2), x_vec(3), &
                    vert_coords(1,icell),         &
                    vert_coords(2,icell))

      end do

    end if

  end subroutine rotate_mesh_coords

  !-----------------------------------------------------------------------------
  !> @brief   Return the rotated north pole using the specified centre of the domain
  !> @details Determines the appropriate rotated north pole based on a target
  !>          rotated domain centre (Null Island).  This constrains the
  !>          rotation such that zero longitude coincides with the
  !>          (unrotated) north pole.
  !>          This is done in degrees (i.e. not radians)
  !>          Note that this should always return a pole in the northern
  !>          hemisphere.
  !>
  !> @return    north_pole  Location of target pole as a
  !>                        [longitude, latitude] array.
  !>
  !> @param[in] null_island Location of target centre (Null Island)
  !>                        [longitude, latitude] array.
  !-----------------------------------------------------------------------------
  function get_target_north_pole(null_island)

    implicit none

    real(kind=r_def), intent(in)  :: null_island(2)
    real(kind=r_def) :: north_pole(2)

    real(kind=r_def) :: get_target_north_pole(2)

    ! Set the target pole
    if (null_island(2) > 0.0_r_def) then
      ! If the null island is in the northern hemisphere, then the
      ! target pole will be offset by 180 degrees longitude
      north_pole(1) = 180.0_r_def + null_island(1)
      north_pole(2) = 90.0_r_def - null_island(2)
    else
      ! If the null island is in the southern hemisphere, then the
      ! target pole will in the northern hemisphere, but on the
      ! same line of longitude
      north_pole(1) = null_island(1)
      north_pole(2) = 90.0_r_def + null_island(2)
    end if

    north_pole(1) = rebase_longitude_range( north_pole(1), -180.0_r_def )

    get_target_north_pole = north_pole

  end function get_target_north_pole

  !-----------------------------------------------------------------------------
  !> @brief   Return the rotated null island from the target rotated pole
  !> @details Determines the appropriate rotated null island based on a target
  !>          rotated north pole).
  !>          This is done in degrees (i.e. not radians)
  !>          Note that this requires the input pole to be in the northern
  !>          hemisphere and should always return a null island in the northern
  !>          hemisphere.
  !>
  !> @return    null_island Location of target null island as a
  !>                        [longitude, latitude] array.
  !>
  !> @param[in] north_pole Location of target pole as a
  !>                        [longitude, latitude] array.
  !-----------------------------------------------------------------------------
  function get_target_null_island(north_pole)

    implicit none

    real(kind=r_def), intent(in)  :: north_pole(2)
    real(kind=r_def) :: null_island(2)

    real(kind=r_def) :: get_target_null_island(2)

    ! Set the target pole
    if (north_pole(2) >= 0.0_r_def) then
      ! null island will be offset by 180 degrees longitude
      null_island(1) = north_pole(1) - 180.0_r_def
      null_island(2) = 90.0_r_def - north_pole(2)
    else
      call log_event( "Target pole must be in the Northern Hemisphere", &
         LOG_LEVEL_ERROR )
    end if

    null_island(1) = rebase_longitude_range( null_island(1), -180.0_r_def )

    get_target_null_island = null_island

  end function get_target_null_island
  !-------------------------------------------------------------------------------

end module rotation_mod
