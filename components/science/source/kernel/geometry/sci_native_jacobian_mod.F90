!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Support routines for computing Jacobians with native coordinate fields
!> @details Contains optimised routines for computing Jacobian matrices, using
!!          native coordinate fields. The first two coordinate fields of the
!!          native coordinates are identical for a whole column, which means
!!          only the lowest values need to be used for a whole column.
!!          This gives a data access optimisation.
module sci_native_jacobian_mod

  use constants_mod,             only: l_def, i_def, r_def, r_single
  use coord_transform_mod,       only: PANEL_ROT_MATRIX, &
                                       alphabetar2xyz,   &
                                       xyz2llr,          &
                                       xyz2ll,           &
                                       llr2xyz,          &
                                       schmidt_transform_lat
  use sci_chi_transform_mod,     only: get_mesh_rotation_matrix, &
                                       get_to_stretch,           &
                                       get_to_rotate,            &
                                       get_stretch_factor

  use finite_element_config_mod, only: coord_system_xyz, &
                                       coord_system_native
  use base_mesh_config_mod,      only: geometry_planar, &
                                       topology_fully_periodic

  implicit none

  private

  public  :: native_jacobian
  ! Public for unit-testing
  public  :: jacobian_stretched
  private :: jacobian_abr2XYZ
  private :: jacobian_llr2XYZ
  private :: jacobian_XYZ2llr

contains

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------

  !> @brief Compute the Jacobian matrices at a 1D array of points (e.g. DoFs)
  !!        for a whole column, using the native coordinates of the mesh
  !! @param[in] coord_system   Finite-element coordinate system enumeration.
  !! @param[in] geometry       Mesh geometry enumeration.
  !! @param[in] topology       Mesh topology enumeration.
  !! @param[in] scaled_radius  Scaled planetary radius.
  !> @param[in] ndf_chi        Num DoFs per cell for coordinate fields
  !> @param[in] nlayers        Number of layers in the mesh
  !> @param[in] chi_1          First native coord field, for a single cell
  !> @param[in] chi_2          Second native coord field, for a single cell
  !> @param[in] chi_3          Third native coord field, for the whole column
  !> @param[in] panel_id       Mesh panel ID value for the column
  !> @param[in] basis          Wchi basis, evaluated at a 1D array of points
  !> @param[in] diff_basis     Derivatives of Wchi basis functions, evaluated at
  !!                           a 1D array of points
  !> @param[in,out] jac        Array of Jacobian matrices to be calculated for
  !!                           a whole column
  !> @param[in,out] dj         Jacobian determinants for the whole column
  subroutine native_jacobian(coord_system, geometry, topology, scaled_radius, &
                             ndf_chi, nlayers, chi_1, chi_2, chi_3, panel_id, &
                             basis, diff_basis, jac, dj)
    implicit none

    integer(kind=i_def),  intent(in) :: coord_system
    integer(kind=i_def),  intent(in) :: geometry
    integer(kind=i_def),  intent(in) :: topology
    real(kind=r_def),     intent(in) :: scaled_radius

    integer(kind=i_def), intent(in)  :: ndf_chi
    integer(kind=i_def), intent(in)  :: nlayers
    integer(kind=i_def), intent(in)  :: panel_id

    real(kind=r_def),    intent(in)    :: chi_1(ndf_chi), chi_2(ndf_chi)
    real(kind=r_def),    intent(in)    :: chi_3(nlayers,ndf_chi)
    real(kind=r_def),    intent(in)    :: basis(1,ndf_chi)
    real(kind=r_def),    intent(in)    :: diff_basis(3,ndf_chi)
    real(kind=r_def),    intent(inout) :: jac(nlayers,3,3)
    real(kind=r_def),    intent(inout) :: dj(nlayers)

    ! Local variables
    real(kind=r_def) :: chi_1_df, chi_2_df
    real(kind=r_def) :: chi_3_df(nlayers)
    real(kind=r_def) :: jac_ref2sph(nlayers,3,3)
    real(kind=r_def) :: jac_sph2XYZ(nlayers,3,3)
    real(kind=r_def) :: lowest_radius
    real(kind=r_def) :: radius
    real(kind=r_def) :: rotation_matrix(3,3)
    real(kind=r_def) :: jac_S(nlayers,3,3)
    real(kind=r_def) :: stretch_factor
    real(kind=r_def) :: native_x, native_y, native_z
    real(kind=r_def) :: native_lon, native_lat

    integer(kind=i_def) :: i, j, k
    logical(kind=l_def) :: to_rotate
    logical(kind=l_def) :: to_stretch

    integer(kind=i_def) :: df, dir

    ! Jacobian from reference element to native coords -------------------------
    jac_ref2sph(:,:,:) = 0.0_r_def
    chi_1_df = 0.0_r_def
    chi_2_df = 0.0_r_def
    chi_3_df(:) = 0.0_r_def
    do dir = 1, 3
      do df = 1, ndf_chi
        jac_ref2sph(:,1,dir) = jac_ref2sph(:,1,dir) + chi_1(df)*diff_basis(dir,df)
        jac_ref2sph(:,2,dir) = jac_ref2sph(:,2,dir) + chi_2(df)*diff_basis(dir,df)
        jac_ref2sph(:,3,dir) = jac_ref2sph(:,3,dir) + chi_3(:,df)*diff_basis(dir,df)
      end do
    end do

    do df = 1, ndf_chi
      chi_1_df = chi_1_df + chi_1(df)*basis(1,df)
      chi_2_df = chi_2_df + chi_2(df)*basis(1,df)
      chi_3_df(:) = chi_3_df(:) + chi_3(:,df)*basis(1,df)
    end do

    ! Jacobian from native to (native) Cartesian coordinates -------------------
    if (coord_system == coord_system_xyz .or. geometry == geometry_planar) then
      ! Using (X,Y,Z) coordinates or on a plane
      jac = jac_ref2sph

    else if (topology == topology_fully_periodic) then
      radius = real(scaled_radius, kind=r_def)
      jac_sph2XYZ = jacobian_abr2XYZ(nlayers, chi_1_df, chi_2_df, chi_3_df+radius, panel_id)

      ! Matrix multiplication of jac_sph2XYZ and jac_ref2sph
      jac(:,:,:) = 0.0_r_def
      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            jac(:,i,j) = jac(:,i,j) + jac_sph2XYZ(:,i,k) * jac_ref2sph(:,k,j)
          end do
        end do
      end do

      ! Apply stretching by Schmidt transform ----------------------------------
      to_stretch = get_to_stretch()
      if (to_stretch) then
        ! Convert chi to spherical polar (un-stretched) coordinates
        lowest_radius = chi_3_df(1)+radius
        call alphabetar2xyz(chi_1_df, chi_2_df, lowest_radius, panel_id,       &
                            native_x, native_y, native_z)
        call xyz2ll(native_x, native_y, native_z,                              &
                    native_lon, native_lat)
        stretch_factor = real(get_stretch_factor(), r_def)
        jac_S = jacobian_stretched(nlayers, native_lon, native_lat, chi_3_df+radius, stretch_factor)

        ! Matrix multiplication of jac_S and jac
        jac_sph2XYZ = jac   ! Store current Jacobian in jac2sphXYZ
        jac(:,:,:) = 0.0_r_def
        do i = 1, 3
          do j = 1, 3
            do k = 1, 3
              jac(:,i,j) = jac(:,i,j) + jac_S(:,i,k) * jac_sph2XYZ(:,k,j)
            end do
          end do
        end do
      end if

      ! Apply rotation ---------------------------------------------------------
      to_rotate = get_to_rotate()
      if (to_rotate) then
        rotation_matrix = get_mesh_rotation_matrix()
        ! Matrix multiplication of rotation_matrix and jac
        jac_sph2XYZ = jac   ! Store current Jacobian in jac2sphXYZ
        jac(:,:,:) = 0.0_r_def
        do i = 1, 3
          do j = 1, 3
            do k = 1, 3
              jac(:,i,j) = jac(:,i,j) + rotation_matrix(i,k) * jac_sph2XYZ(:,k,j)
            end do
          end do
        end do
      end if

    else
      ! Native coordinates for a limited area domain on the sphere
      radius = real(scaled_radius, kind=r_def)
      jac_sph2XYZ = jacobian_llr2XYZ(nlayers, chi_1_df, chi_2_df, chi_3_df+radius)

      ! Matrix multiplication of jac_sph2XYZ and jac_ref2sph
      jac(:,:,:) = 0.0_r_def
      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            jac(:,i,j) = jac(:,i,j) + jac_sph2XYZ(:,i,k) * jac_ref2sph(:,k,j)
          end do
        end do
      end do

      ! Apply rotation ---------------------------------------------------------
      to_rotate = get_to_rotate()
      if (to_rotate) then
        rotation_matrix = get_mesh_rotation_matrix()
        ! Matrix multiplication of rotation_matrix and jac
        jac_sph2XYZ = jac   ! Store current Jacobian in jac2sphXYZ
        jac(:,:,:) = 0.0_r_def
        do i = 1, 3
          do j = 1, 3
            do k = 1, 3
              jac(:,i,j) = jac(:,i,j) + rotation_matrix(i,k) * jac_sph2XYZ(:,k,j)
            end do
          end do
        end do
      end if

    end if

    ! Compute determinant ------------------------------------------------------
    dj(:) = jac(:,1,1) * (jac(:,2,2)*jac(:,3,3) - jac(:,2,3)*jac(:,3,2))       &
          - jac(:,1,2) * (jac(:,2,1)*jac(:,3,3) - jac(:,2,3)*jac(:,3,1))       &
          + jac(:,1,3) * (jac(:,2,1)*jac(:,3,2) - jac(:,2,2)*jac(:,3,1))

  end subroutine native_jacobian

  ! -------------------------------------------------------------------------- !
  ! Jacobian for transforming from equiangular cubed-sphere to Cartesian coords
  ! -------------------------------------------------------------------------- !
  !> @brief Compute the (lon,lat,r) -> (X,Y,Z) Jacobian for a column
  !> @param[in] nlayers    Number of layers in the mesh
  !> @param[in] alpha      Alpha coordinate in native system
  !> @param[in] beta       Beta coordinate in native system
  !> @param[in] radius     Radial coordinate
  !> @param[in] panel_id   Mesh panel ID value for the column
  !> @return    jac_out    3x3 matrix for the Jacobian of the transformation
  function jacobian_abr2XYZ(nlayers, alpha, beta, radius, panel_id) result(jac_out)
    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: alpha
    real(kind=r_def),    intent(in) :: beta
    real(kind=r_def),    intent(in) :: radius(nlayers)
    integer(kind=i_def), intent(in) :: panel_id

    real(kind=r_def) :: jac_abr2XYZ(nlayers,3,3)
    real(kind=r_def) :: jac_out(nlayers,3,3)
    real(kind=r_def) :: tan_alpha, tan_beta, panel_rho, inv_prho3

    integer(kind=i_def) :: i, j, k

    tan_alpha = tan(alpha)
    tan_beta = tan(beta)
    panel_rho = sqrt(1.0_r_def + tan_alpha**2 + tan_beta**2)
    inv_prho3 = 1.0_r_def / panel_rho**3

    ! First column, g_alpha
    jac_abr2XYZ(:,1,1) = -radius(:)*tan_alpha*(1.0_r_def + tan_alpha**2)*inv_prho3
    jac_abr2XYZ(:,2,1) = radius(:)*(1.0_r_def + tan_beta**2)*(1.0_r_def + tan_alpha**2)*inv_prho3
    jac_abr2XYZ(:,3,1) = -radius(:)*tan_alpha*tan_beta*(1.0_r_def + tan_alpha**2)*inv_prho3

    ! Second column, g_beta
    jac_abr2XYZ(:,1,2) = -radius(:)*tan_beta*(1.0_r_def + tan_beta**2)*inv_prho3
    jac_abr2XYZ(:,2,2) = -radius(:)*tan_alpha*tan_beta*(1.0_r_def + tan_beta**2)*inv_prho3
    jac_abr2XYZ(:,3,2) = radius(:)*(1.0_r_def + tan_alpha**2)*(1.0_r_def + tan_beta**2)*inv_prho3

    ! Third column, g_r
    jac_abr2XYZ(:,1,3) = 1.0_r_def/panel_rho
    jac_abr2XYZ(:,2,3) = tan_alpha/panel_rho
    jac_abr2XYZ(:,3,3) = tan_beta/panel_rho

    ! Rotate to the appropriate panel
    ! Matrix multiplication of panel rotation matrix with panel 1 Jacobian
    jac_out(:,:,:) = 0.0_r_def
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          jac_out(:,i,j) = jac_out(:,i,j) + PANEL_ROT_MATRIX(i,k,panel_id) * jac_abr2XYZ(:,k,j)
        end do
      end do
    end do

  end function jacobian_abr2XYZ

  ! -------------------------------------------------------------------------- !
  ! Jacobian for transforming from spherical polar to Cartesian coordinates
  ! -------------------------------------------------------------------------- !
  !> @brief Compute the (lon,lat,r) -> (X,Y,Z) Jacobian for a column
  !> @param[in] nlayers    Number of layers in the mesh
  !> @param[in] longitude  Longitudinal coordinate
  !> @param[in] latitude   Latitudinal coordinate
  !> @param[in] radius     Radial coordinate
  !> @return    jac_llr2XYZ  3x3 matrix for the Jacobian of the transformation
  function jacobian_llr2XYZ(nlayers, longitude, latitude, radius) &
                                                      result(jac_llr2XYZ)
    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: longitude
    real(kind=r_def),    intent(in) :: latitude
    real(kind=r_def),    intent(in) :: radius(nlayers)

    real(kind=r_def) :: jac_llr2XYZ(nlayers,3,3)
    real(kind=r_def) :: sin_lon, sin_lat, cos_lon, cos_lat

    sin_lat = sin(latitude)
    sin_lon = sin(longitude)
    cos_lat = cos(latitude)
    cos_lon = cos(longitude)

    ! First column, g_lon
    jac_llr2XYZ(:,1,1) = -radius(:)*sin_lon*cos_lat
    jac_llr2XYZ(:,2,1) = radius(:)*cos_lon*cos_lat
    jac_llr2XYZ(:,3,1) = 0.0_r_def

    ! Second column, g_lat
    jac_llr2XYZ(:,1,2) = -radius(:)*cos_lon*sin_lat
    jac_llr2XYZ(:,2,2) = -radius(:)*sin_lon*sin_lat
    jac_llr2XYZ(:,3,2) = radius(:)*cos_lat

    ! Third column, g_r
    jac_llr2XYZ(:,1,3) = cos_lon*cos_lat
    jac_llr2XYZ(:,2,3) = sin_lon*cos_lat
    jac_llr2XYZ(:,3,3) = sin_lat

  end function jacobian_llr2XYZ

  ! -------------------------------------------------------------------------- !
  ! Jacobian for transforming from Cartesian to spherical polar coordinates
  ! -------------------------------------------------------------------------- !
  !> @brief Compute the (X,Y,Z) -> (lon,lat,r) Jacobian for a column
  !> @param[in] nlayers    Number of layers in the mesh
  !> @param[in] longitude  Longitudinal coordinate
  !> @param[in] latitude   Latitudinal coordinate
  !> @param[in] radius     Radial coordinate
  !> @return    jac_XYZ2llr  3x3 matrix for the Jacobian of the transformation
  function jacobian_XYZ2llr(nlayers, longitude, latitude, radius) &
                                                      result(jac_XYZ2llr)
    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: longitude
    real(kind=r_def),    intent(in) :: latitude
    real(kind=r_def),    intent(in) :: radius(nlayers)

    real(kind=r_def)             :: jac_XYZ2llr(nlayers,3,3)
    real(kind=r_def)             :: sin_lon, sin_lat, cos_lon, cos_lat
    real(kind=r_def)             :: safe_cos_lat
    real(kind=r_def),  parameter :: tiny = 1.0e-15_r_def

    sin_lat = sin(latitude)
    sin_lon = sin(longitude)
    cos_lat = cos(latitude)
    cos_lon = cos(longitude)

    ! To avoid divide by zero errors at poles, add tiny number to cos(lat)
    safe_cos_lat = cos_lat + sign(tiny, cos_lat)

    jac_XYZ2llr(:,1,1) = -sin_lon / (radius(:) * safe_cos_lat)
    jac_XYZ2llr(:,1,2) = cos_lon / (radius(:) * safe_cos_lat)
    jac_XYZ2llr(:,1,3) = 0.0_r_single
    jac_XYZ2llr(:,2,1) = - cos_lon * sin_lat / radius(:)
    jac_XYZ2llr(:,2,2) = - sin_lon * sin_lat / radius(:)
    jac_XYZ2llr(:,2,3) = cos_lat / radius(:)
    jac_XYZ2llr(:,3,1) = cos_lon * cos_lat
    jac_XYZ2llr(:,3,2) = sin_lon * cos_lat
    jac_XYZ2llr(:,3,3) = sin_lat

  end function jacobian_XYZ2llr

  ! -------------------------------------------------------------------------- !
  ! Jacobian for Schmidt transform
  ! -------------------------------------------------------------------------- !
  !> @brief Compute the Jacobian for performing Schmidt transform for a column
  !> @param[in] nlayers    Number of layers in the mesh
  !> @param[in] longitude  Longitudinal coordinate in native (stretched) system
  !> @param[in] latitude   Latitudinal coordinate in native (stretched) system
  !> @param[in] radius     The radial coordinate
  !> @param[in] stretch    The stretching factor
  !> @return    jac_stretched  3x3 matrix for the Jacobian of the transformation
  function jacobian_stretched(nlayers, longitude, latitude, radius, stretch) result(jac_stretched)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: longitude, latitude
    real(kind=r_def),    intent(in) :: radius(nlayers), stretch

    real(kind=r_def) :: jac_llr2XYZ(nlayers,3,3)
    real(kind=r_def) :: jac_XYZ2llr(nlayers,3,3)

    real(kind=r_def), parameter :: one = 1.0_r_def

    real(kind=r_def) :: lat_stretched, psi
    real(kind=r_def) :: jac_stretched(nlayers,3,3)

    integer(kind=i_def) :: i, j, k

    ! Compute stretched variables
    lat_stretched = schmidt_transform_lat(latitude, stretch)
    psi = 2.0_r_def*stretch / (one + stretch**2 + (one - stretch**2)*sin(latitude))

    ! Get Jacobian for transformation from (X,Y,Z) to (lon,lat,r) coords
    ! Stretching Jacobian is:
    ! ( 1  0  0 )
    ! ( 0 psi 0 )
    ! ( 0  0  1 )
    ! So don't need to multiply out the whole matrix
    jac_XYZ2llr = jacobian_XYZ2llr(nlayers, longitude, latitude, radius)
    jac_XYZ2llr(:,2,1) = psi*jac_XYZ2llr(:,2,1)
    jac_XYZ2llr(:,2,2) = psi*jac_XYZ2llr(:,2,2)
    jac_XYZ2llr(:,2,3) = psi*jac_XYZ2llr(:,2,3)

    ! Get Jacobian for transformation from (lon,lat,r) to (X,Y,Z) coords on
    ! the stretched mesh
    jac_llr2XYZ = jacobian_llr2XYZ(nlayers, longitude, lat_stretched, radius)

    ! The resulting Jacobian is the product of the previous two Jacobians
    jac_stretched(:,:,:) = 0.0_r_def
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          jac_stretched(:,i,j) = jac_stretched(:,i,j) + jac_llr2XYZ(:,i,k) * jac_XYZ2llr(:,k,j)
        end do
      end do
    end do

  end function jacobian_stretched

end module sci_native_jacobian_mod
