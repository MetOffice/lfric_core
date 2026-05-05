!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!>  @brief   Routines for transforming the chi coordinate fields
!!
!!  @details Contains routines for conversion of chi coordinate fields. These
!!           are accessed through the chi2ABC interface functions, so that
!!           which coord_system chi is in, it will convert the
!!           coordinates to the ABC system
!------------------------------------------------------------------------------
module sci_chi_transform_mod

use constants_mod,             only : r_def, i_def, l_def,     &
                                      str_def, EPS, PI, rmdi
use coord_transform_mod,       only : alphabetar2xyz,          &
                                      alphabetar2llr,          &
                                      xyz2alphabetar,          &
                                      llr2xyz, xyz2llr,        &
                                      xyz2ll,                  &
                                      mesh_rotation_matrix,    &
                                      schmidt_transform_xyz,   &
                                      inverse_schmidt_transform_xyz
use log_mod,                   only : log_event,               &
                                      log_scratch_space,       &
                                      LOG_LEVEL_ERROR,         &
                                      LOG_LEVEL_DEBUG,         &
                                      LOG_LEVEL_WARNING
use matrix_invert_mod,         only : matrix_invert_3x3

use base_mesh_config_mod,      only : geometry,                &
                                      geometry_spherical,      &
                                      geometry_planar,         &
                                      topology,                &
                                      topology_fully_periodic
use finite_element_config_mod, only : coord_system,            &
                                      coord_system_xyz,        &
                                      coord_system_native
use planet_config_mod,         only : scaled_radius

implicit none

private

! ---------------------------------------------------------------------------- !
! Private matrices or values that need computing once
! ---------------------------------------------------------------------------- !

real(kind=r_def)    :: chi2xyz_rot_mat(3,3)
real(kind=r_def)    :: xyz2chi_rot_mat(3,3)
real(kind=r_def)    :: stretch_factor
logical(kind=l_def) :: to_rotate
logical(kind=l_def) :: to_stretch

! ---------------------------------------------------------------------------- !
! Public subroutines
! ---------------------------------------------------------------------------- !
public :: init_chi_transforms
public :: final_chi_transforms
public :: chi2xyz
public :: chi2abr
public :: chi2llr
public :: chir2xyz
public :: get_mesh_rotation_matrix
public :: get_inverse_mesh_rotation_matrix
public :: get_stretch_factor
public :: get_to_rotate
public :: get_to_stretch

!------------------------------------------------------------------------------
! Contained functions / subroutines
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!> @brief  Initialise the coordinate transform information
!!
!> @param[in] mesh_collection    Optional: a collection of meshes, which contain
!!                               metadata used to determine the rotation matrix
!!                               and stretching factors.
!> @param[in] north_pole_arg     Optional: target north pole, used to generate
!!                               the rotation matrix. This is incompatible with
!!                               the mesh_collection argument, and ideally
!!                               should only be used for unit-testing.
!> @param[in] equator_lat_arg    Optional: Latitude of the equator of the mesh,
!!                               allowing a stretching to be described.
!!                               This is incompatible with the mesh_collection
!!                               argument, and ideally should only be used for
!!                               unit-testing.
!------------------------------------------------------------------------------
subroutine init_chi_transforms( geometry, topology, &
                                mesh_collection,    &
                                north_pole_arg, equator_lat_arg )

  use local_mesh_mod,            only : local_mesh_type
  use mesh_collection_mod,       only : mesh_collection_type
  use mesh_mod,                  only : mesh_type

  implicit none

  integer(i_def), intent(in) :: geometry
  integer(i_def), intent(in) :: topology

  type(mesh_collection_type), optional, intent(in) :: mesh_collection
  real(kind=r_def),           optional, intent(in) :: north_pole_arg(2)
  real(kind=r_def),           optional, intent(in) :: equator_lat_arg

  type(mesh_type),         pointer :: mesh
  type(local_mesh_type),   pointer :: local_mesh
  character(str_def),  allocatable :: all_mesh_names(:)

  real(kind=r_def) :: north_pole(2)
  real(kind=r_def) :: null_island(2)
  real(kind=r_def) :: equatorial_latitude

  ! -------------------------------------------------------------------------- !
  ! Extract stretching and rotation information from mesh
  ! -------------------------------------------------------------------------- !

  ! Begin by assuming no stretching and no rotation
  to_stretch = .false.
  to_rotate = .false.
  north_pole(1) = PI
  north_pole(2) = PI/2.0_r_def
  null_island(1) = 0.0_r_def
  null_island(2) = 0.0_r_def
  equatorial_latitude = 0.0_r_def

  if ( present(mesh_collection) .and.                                          &
       (present(equator_lat_arg) .or. present(north_pole_arg)) ) then
    call log_event(                                                            &
      'init_chi_transform: mesh_compatible argument cannot be passed with ' // &
      'another argument', LOG_LEVEL_ERROR                                      &
    )
  end if

  if (present(mesh_collection)) then
    ! NB:
    ! At this stage, we will assume that the stretching and rotation are the same
    ! for all meshes. If they weren't, we would need to extract a different factor
    ! and different rotation matrix for each mesh. The chi2*** transforms would
    ! also need to take mesh_id as an argument, which would be a major API change
    ! since it would need passing through each kernel.
    ! Therefore, extract first mesh from collection ...
    all_mesh_names = mesh_collection%get_mesh_names()
    if (SIZE(all_mesh_names) > 0) then
      mesh => mesh_collection%get_mesh(all_mesh_names(1))
    else
      call log_event(                                                          &
        'init_chi_transform: unable to determine mesh rotation and ' //        &
        'stretching because there are no meshes!', LOG_LEVEL_ERROR             &
      )
    end if

    ! Extract rotation and stretching information from global mesh
    local_mesh => mesh%get_local_mesh()
    north_pole = local_mesh%get_north_pole()
    null_island = local_mesh%get_null_island()
    equatorial_latitude = local_mesh%get_equatorial_latitude()

    ! If any variables are unset, set them to defaults here --------------------
    if ( abs(north_pole(1) - rmdi) < EPS                                       &
         .or. abs(north_pole(2) - rmdi) < EPS ) then
      north_pole(1) = 0.0_r_def
      north_pole(2) = PI/2.0_r_def
      call log_event(                                                          &
        'Mesh North Pole not set, so using (lon=0, lat=pi/2) as default',      &
         LOG_LEVEL_WARNING                                                     &
      )
    end if
    if ( abs(null_island(1) - rmdi) < EPS                                      &
         .or. abs(null_island(2) - rmdi) < EPS ) then
      null_island(1) = 0.0_r_def
      null_island(2) = 0.0_r_def
      call log_event(                                                          &
        'Mesh Null Island not set, so using (lon=0, lat=0) as default',        &
         LOG_LEVEL_WARNING                                                     &
      )
    end if
    if ( abs(equatorial_latitude - rmdi) < EPS .or.                            &
         geometry == geometry_planar .or. topology /= topology_fully_periodic ) then
      equatorial_latitude = 0.0_r_def
      call log_event(                                                          &
        'Equatorial latitude for mesh not set, so using 0.0 as default',       &
         LOG_LEVEL_WARNING                                                     &
      )
    end if
  end if

  if (present(north_pole_arg)) north_pole = north_pole_arg
  if (present(equator_lat_arg)) equatorial_latitude = equator_lat_arg


  ! Now that parameters have been read in, determine if stretching or rotation
  ! are actually happening
  to_stretch = abs(equatorial_latitude) > EPS
  ! It's probably safer to check both the null island and the north pole here
  to_rotate = ( abs(north_pole(2) - PI/2.0_r_def) > EPS                        &
                .or. abs(null_island(1)) > EPS .or. abs(null_island(2)) > EPS )

  ! Compute Schmidt stretch factor ---------------------------------------------
  stretch_factor = sqrt( (1.0_r_def - sin(equatorial_latitude))                &
                         / (1.0_r_def + sin(equatorial_latitude)) )

  ! Compute rotation matrix ----------------------------------------------------
  chi2xyz_rot_mat = mesh_rotation_matrix(north_pole)

  ! Compute inverse rotation matrix --------------------------------------------
  xyz2chi_rot_mat = matrix_invert_3x3(chi2xyz_rot_mat)

  write(log_scratch_space,'(A,L6,A,2E16.8)')                                   &
    'Mesh rotation: ', to_rotate, ' north pole: ', north_pole(1), north_pole(2)
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
  write(log_scratch_space,'(A,L6,A,E16.8,A,E16.8)')                            &
    'Mesh stretching: ', to_stretch, ' stretching factor: ', stretch_factor,   &
    '   latitude of equator: ', equatorial_latitude
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

end subroutine init_chi_transforms

!------------------------------------------------------------------------------
!>  @brief  Nullify the coordinate transform values
!------------------------------------------------------------------------------
subroutine final_chi_transforms()

  implicit none

  to_stretch = .false.
  to_rotate = .false.
  stretch_factor = rmdi
  chi2xyz_rot_mat(:,:) = 0.0_r_def
  xyz2chi_rot_mat(:,:) = 0.0_r_def

end subroutine final_chi_transforms


!-------------------------------------------------------------------------------
!> @brief Transforms a coordinate field chi from any system into global
!>        Cartesian (X,Y,Z) coordinates. If chi is in a spherical coordinate
!>        system, the third coordinate should be height, and the scaled_radius
!>        will be added to the height to give the radius before the coordinates
!>        are transformed to (X,Y,Z) coordinates.
!!
!! @param[in]   chi_1      The first coordinate field in
!! @param[in]   chi_2      The second coordinate field in
!! @param[in]   chi_3      The third coordinate field in
!! @param[in]   panel_id   The mesh panel ID
!! @param[out]  x          The first coordinate field out (global Cartesian X)
!! @param[out]  y          The second coordinate field out (global Cartesian Y)
!! @param[out]  z          The third coordinate field out (global Cartesian Z)
!-------------------------------------------------------------------------------
subroutine chi2xyz(chi_1, chi_2, chi_3, panel_id, x, y, z)

  implicit none

  integer(kind=i_def), intent(in)  :: panel_id
  real(kind=r_def),    intent(in)  :: chi_1, chi_2, chi_3
  real(kind=r_def),    intent(out) :: x, y, z

  real(kind=r_def) :: xyz(3)

  if (geometry == geometry_planar .or. coord_system == coord_system_xyz) then
    ! chi already uses (geocentric) Cartesian coordinates
    x = chi_1
    y = chi_2
    z = chi_3

  else if (topology /= topology_fully_periodic) then
    ! domain is a spherical LAM, using (lon,lat,z) coordinates
    call llr2xyz(chi_1, chi_2, chi_3+scaled_radius, x, y, z)

    if (to_rotate) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = matmul(chi2xyz_rot_mat, xyz)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if

  else
    ! cubed-sphere coordinates
    ! transform to native (X,Y,Z) coordinates
    call alphabetar2xyz(chi_1, chi_2, chi_3+scaled_radius, panel_id, x, y, z)

    ! stretch, if necessary
    if (to_stretch) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = schmidt_transform_xyz(xyz, stretch_factor)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if

    ! rotate, if necessary
    if (to_rotate) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = matmul(chi2xyz_rot_mat, xyz)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if
  end if

end subroutine chi2xyz


!-------------------------------------------------------------------------------
!> @brief Transforms a coordinate field chi from any system into global
!>        Cartesian (X,Y,Z) coordinates. If chi is in a spherical coordinate
!>        system, the third coordinate should be radius (distinguishing this
!>        function from chi2xyz above). Therefore this will not add the
!>        scaled_radius to transform.
!!
!! @param[in]   chi_1      The first coordinate field in
!! @param[in]   chi_2      The second coordinate field in
!! @param[in]   chi_3      The third coordinate field in
!! @param[in]   panel_id   The mesh panel ID
!! @param[out]  x          The first coordinate field out (global Cartesian X)
!! @param[out]  y          The second coordinate field out (global Cartesian Y)
!! @param[out]  z          The third coordinate field out (global Cartesian Z)
!-------------------------------------------------------------------------------
subroutine chir2xyz(chi_1, chi_2, chi_3, panel_id, x, y, z)

  implicit none

  integer(kind=i_def), intent(in)  :: panel_id
  real(kind=r_def),    intent(in)  :: chi_1, chi_2, chi_3
  real(kind=r_def),    intent(out) :: x, y, z

  real(kind=r_def) :: xyz(3)

  if (coord_system == coord_system_xyz .or. geometry == geometry_planar) then
    ! chi already uses (geocentric) Cartesian coordinates
    x = chi_1
    y = chi_2
    z = chi_3

  else if (topology /= topology_fully_periodic) then
    ! domain is a spherical LAM, using (lon,lat,z) coordinates
    call llr2xyz(chi_1, chi_2, chi_3, x, y, z)

    if (to_rotate) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = matmul(chi2xyz_rot_mat, xyz)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if

  else
    ! cubed-sphere coordinates
    ! transform to native (X,Y,Z) coordinates
    call alphabetar2xyz(chi_1, chi_2, chi_3, panel_id, x, y, z)

    ! stretch, if necessary
    if (to_stretch) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = schmidt_transform_xyz(xyz, stretch_factor)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if

    ! rotate, if necessary
    if (to_rotate) then
      xyz(1) = x
      xyz(2) = y
      xyz(3) = z

      xyz = matmul(chi2xyz_rot_mat, xyz)

      x = xyz(1)
      y = xyz(2)
      z = xyz(3)
    end if
  end if

end subroutine chir2xyz


!-------------------------------------------------------------------------------
!> @brief Transforms a coordinate field chi from any system into spherical polar
!>        (longitude, latitude, radius) coordinates
!!
!! @param[in]   chi_1      The first coordinate field in
!! @param[in]   chi_2      The second coordinate field in
!! @param[in]   chi_3      The third coordinate field in
!! @param[in]   panel_id   The mesh panel ID
!! @param[out]  longitude  The first coordinate field out (longitude)
!! @param[out]  latitude   The second coordinate field out (latitude)
!! @param[out]  radius     The third coordinate field out (radius)
!-------------------------------------------------------------------------------
subroutine chi2llr(chi_1, chi_2, chi_3, panel_id, lon, lat, radius)

  implicit none

  integer(kind=i_def), intent(in)  :: panel_id
  real(kind=r_def),    intent(in)  :: chi_1, chi_2, chi_3
  real(kind=r_def),    intent(out) :: lon, lat, radius

  real(kind=r_def) :: xyz(3)

  if (geometry == geometry_planar .or. coord_system == coord_system_xyz) then
    ! chi uses (geocentric) Cartesian coordinates
    call xyz2llr(chi_1, chi_2, chi_3, lon, lat, radius)

  else if (topology /= topology_fully_periodic) then
    ! domain is a spherical LAM, already using (lon,lat,z) coordinates
    ! may need to rotate these to the physical (lon,lat) coordinates

    ! avoid conversions in computing radius
    radius = chi_3 + scaled_radius

    if (to_rotate) then
      call llr2xyz(chi_1, chi_2, radius, xyz(1), xyz(2), xyz(3))
      xyz = matmul(chi2xyz_rot_mat, xyz)
      call xyz2ll(xyz(1), xyz(2), xyz(3), lon, lat)
    else
      lon = chi_1
      lat = chi_2
    end if

  else
    ! cubed-sphere coordinates
    ! transform to native (X,Y,Z) coordinates
    radius = chi_3 + scaled_radius

    if (to_stretch .or. to_rotate) then
      call alphabetar2xyz(chi_1, chi_2, radius, panel_id, xyz(1), xyz(2), xyz(3))

      ! stretch, if necessary
      if (to_stretch) then
        xyz = schmidt_transform_xyz(xyz, stretch_factor)
      end if

      ! rotate, if necessary
      if (to_rotate) then
        xyz = matmul(chi2xyz_rot_mat, xyz)
      end if

      ! convert to spherical polar coordinates
      call xyz2ll(xyz(1), xyz(2), xyz(3), lon, lat)

    else
      call alphabetar2llr(chi_1, chi_2, radius, panel_id, lon, lat)
    end if

  end if

end subroutine chi2llr


!-------------------------------------------------------------------------------
!> @brief Transforms a coordinate field chi from any system into *native*
!!        equiangular cubed sphere (alpha,beta,radius) coordinates
!!
!! @param[in]   chi_1      The first coordinate field in
!! @param[in]   chi_2      The second coordinate field in
!! @param[in]   chi_3      The third coordinate field in
!! @param[in]   panel_id   The mesh panel ID
!! @param[out]  alpha      The first coordinate field out (alpha)
!! @param[out]  beta       The second coordinate field out (beta)
!! @param[out]  radius     The third coordinate field out (radius)
!-------------------------------------------------------------------------------
subroutine chi2abr(chi_1, chi_2, chi_3, panel_id, alpha, beta, radius)

  implicit none

  integer(kind=i_def), intent(in)  :: panel_id
  real(kind=r_def),    intent(in)  :: chi_1, chi_2, chi_3
  real(kind=r_def),    intent(out) :: alpha, beta, radius

  real(kind=r_def) :: xyz(3)

  if (topology /= topology_fully_periodic .or. geometry /= geometry_spherical) then
    call log_event(                                                            &
      'chi2abr can only be used on cubed-sphere meshes', LOG_LEVEL_ERROR       &
    )

  else if (coord_system == coord_system_native) then
    alpha = chi_1
    beta = chi_2
    radius = chi_3 + scaled_radius

  else
    ! geocentric Cartesian coordinates
    xyz(1) = chi_1
    xyz(2) = chi_2
    xyz(3) = chi_3

    ! un-rotate, if necessary
    if (to_rotate) then
      xyz = matmul(xyz2chi_rot_mat, xyz)
    end if

    ! un-stretch, if necessary
    if (to_stretch) then
      xyz = inverse_schmidt_transform_xyz(xyz, stretch_factor)
    end if

    ! transform to equiangular cubed-sphere coordinates
    call xyz2alphabetar(xyz(1), xyz(2), xyz(3), panel_id, alpha, beta, radius)
  end if

end subroutine chi2abr

!-------------------------------------------------------------------------------
!> @brief Returns a pointer to the rotation matrix for transforming from the
!!        native Cartesian coordinates to the physical Cartesian coordinates
!-------------------------------------------------------------------------------
function get_mesh_rotation_matrix() result(rot_mat)
  implicit none
  real(kind=r_def) :: rot_mat(3,3)

  rot_mat = chi2xyz_rot_mat

end function get_mesh_rotation_matrix

!-------------------------------------------------------------------------------
!> @brief Returns a pointer to the inverse rotation matrix, transforming from
!!        physical Cartesian coordinates to native Cartesian coordinates
!-------------------------------------------------------------------------------
function get_inverse_mesh_rotation_matrix() result(rot_mat)
  implicit none
  real(kind=r_def) :: rot_mat(3,3)

  rot_mat = xyz2chi_rot_mat

end function get_inverse_mesh_rotation_matrix

!-------------------------------------------------------------------------------
!> @brief Returns the Schmidt transform stretch factor
!-------------------------------------------------------------------------------
function get_stretch_factor() result(stretch_factor_out)
  implicit none
  real(kind=r_def) :: stretch_factor_out

  stretch_factor_out = stretch_factor

end function get_stretch_factor

!-------------------------------------------------------------------------------
!> @brief Returns whether coordinates are rotated
!-------------------------------------------------------------------------------
function get_to_rotate() result(to_rotate_out)
  implicit none
  logical(kind=l_def) :: to_rotate_out

  to_rotate_out = to_rotate

end function get_to_rotate

!-------------------------------------------------------------------------------
!> @brief Returns whether coordinates are stretched
!-------------------------------------------------------------------------------
function get_to_stretch() result(to_stretch_out)
  implicit none
  logical(kind=l_def) :: to_stretch_out

  to_stretch_out = to_stretch

end function get_to_stretch

end module sci_chi_transform_mod

