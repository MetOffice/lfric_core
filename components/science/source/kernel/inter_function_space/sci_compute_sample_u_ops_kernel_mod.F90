!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes operators for sampling a W2 field from scalar components
!> @details This kernel builds the operators for converting three scalar
!!          components in (W3,W3,Wtheta) into a vector-valued field in W2, via
!!          sampling. The input components will take physical values, while the
!!          returned W2 field is in "computational" space.
!!          For spherical geometries, the components are assumed to be in
!!          spherical polar coordinates.
!!          This is only designed for the lowest-order function spaces.

module sci_compute_sample_u_ops_kernel_mod

  use argument_mod,            only : arg_type, func_type,         &
                                      GH_FIELD, GH_REAL,           &
                                      GH_OPERATOR,                 &
                                      GH_INC, GH_READ, GH_WRITE,   &
                                      ANY_DISCONTINUOUS_SPACE_3,   &
                                      GH_BASIS, GH_DIFF_BASIS,     &
                                      CELL_COLUMN, GH_EVALUATOR,   &
                                      reference_element_data_type, &
                                      normals_to_faces
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2broken, W3, Wtheta, Wchi
  use kernel_mod,              only : kernel_type
  use sci_chi_transform_mod,   only : chi2llr
  use sci_coordinate_jacobian_mod, only : coordinate_jacobian, &
                                          coordinate_jacobian_inverse
  use coord_transform_mod,     only : sphere2cart_vector
  use reference_element_mod,   only : W, S, N, E, T, B

  use finite_element_config_mod, only: coord_system
  use base_mesh_config_mod,      only: geometry, topology, &
                                       geometry_spherical, &
                                       geometry_planar
  use planet_config_mod,         only: scaled_radius

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_sample_u_ops_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2broken, W3),               &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2broken, W3),               &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2broken, WTHETA),           &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  Wchi),                       &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)   &
         /)
    type(func_type) :: meta_funcs(1) = (/                                      &
         func_type(Wchi, GH_BASIS, GH_DIFF_BASIS)                              &
         /)
    type(reference_element_data_type), dimension(1) ::                         &
    meta_reference_element =                                                   &
      (/ reference_element_data_type(normals_to_faces) /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: compute_sample_u_ops_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: compute_sample_u_ops_code

contains

!> @brief Computes operators for sampling three scalar wind components into W2
!> @param[in]     col                      Index of column
!> @param[in]     nlayers                  Number of layers in the mesh
!> @param[in]     ncells_3d_1              Number of cells in this partition of 3D mesh
!> @param[in,out] u_lon_op                 W3 -> W2 operator for longitudunal component
!> @param[in]     ncells_3d_1              Number of cells in this partition of 3D mesh
!> @param[in,out] u_lat_op                 W3 -> W2 operator for latitudinal component
!> @param[in]     ncells_3d_1              Number of cells in this partition of 3D mesh
!> @param[in,out] u_rad_op                 W3 -> W2 operator for radial component
!> @param[in]     chi1                     Coordinates in the first direction
!> @param[in]     chi2                     Coordinates in the second direction
!> @param[in]     chi3                     Coordinates in the third direction
!> @param[in]     panel_id                 A field giving the ID for mesh panels
!> @param[in]     ndf_w2b                  Number of DoFs per cell for broken W2
!> @param[in]     ndf_w3                   Number of DoFs per cell for W3
!> @param[in]     ndf_wt                   Number of DoFs per cell for Wtheta
!> @param[in]     ndf_chi                  Number of DoFs per cell for chi fields
!> @param[in]     undf_chi                 Number of unique DoFs per partition for chi fields
!> @param[in]     map_chi                  DoF map for the column's base cell for chi fields
!> @param[in]     chi_basis                Basis functions of the coordinate space
!!                                         evaluated at W2 nodal points
!> @param[in]     chi_diff_basis           Differential basis functions of the coordinate
!!                                         space evaluated at W2 nodal points
!> @param[in]     ndf_pid                  Number of DoFs per cell for panel ID field
!> @param[in]     undf_pid                 Number of unique DoFs per partition for panel ID field
!> @param[in]     map_pid                  DoF map for the column's base cell for panel ID field
!> @param[in]     nfaces                   Number of cell faces
!> @param[in]     face_normals             The normal vectors to each face
subroutine compute_sample_u_ops_code( col, nlayers,                   &
                                      ncell_3d_1, u_lon_op,           &
                                      ncell_3d_2, u_lat_op,           &
                                      ncell_3d_3, u_rad_op,           &
                                      chi1, chi2, chi3,               &
                                      panel_id,                       &
                                      ndf_w2b, ndf_w3, ndf_wt,        &
                                      ndf_chi, undf_chi, map_chi,     &
                                      chi_basis, chi_diff_basis,      &
                                      ndf_pid, undf_pid, map_pid,     &
                                      nfaces, face_normals            &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: col, nlayers, nfaces
  integer(kind=i_def), intent(in) :: ncell_3d_1, ncell_3d_2, ncell_3d_3
  integer(kind=i_def), intent(in) :: ndf_chi, ndf_w2b, ndf_w3, ndf_wt, ndf_pid
  integer(kind=i_def), intent(in) :: undf_pid, undf_chi

  ! DoF maps
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  ! Basis functions
  real(kind=r_def), intent(in), dimension(1,ndf_chi,ndf_w2b) :: chi_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,ndf_w2b) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(3,nfaces)          :: face_normals

  ! Fields
  real(kind=r_def), dimension(undf_pid), intent(in) :: panel_id
  real(kind=r_def), dimension(undf_chi), intent(in) :: chi1, chi2, chi3

  ! Operators
  real(kind=r_def), dimension(ncell_3d_1,ndf_w2b,ndf_w3), intent(inout) :: u_lon_op
  real(kind=r_def), dimension(ncell_3d_2,ndf_w2b,ndf_w3), intent(inout) :: u_lat_op
  real(kind=r_def), dimension(ncell_3d_3,ndf_w2b,ndf_wt), intent(inout) :: u_rad_op

  ! Internal variables
  integer(kind=i_def) :: df_w2, df_wt, df_chi, k, ipanel, cell_3d
  real(kind=r_def), dimension(3,3,ndf_w2b) :: jacobian, jac_inv
  real(kind=r_def), dimension(ndf_w2b)     :: dj
  real(kind=r_def), dimension(3)           :: llr, X_vector, Y_vector, Z_vector
  real(kind=r_def), dimension(3)           :: lon_vector_llr, lat_vector_llr
  real(kind=r_def), dimension(3)           :: rad_vector_llr
  real(kind=r_def), dimension(ndf_chi)     :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(ndf_w2b)     :: chi1_w2, chi2_w2, chi3_w2
  real(kind=r_def), dimension(ndf_w2b,3)   :: lon_vector_xyz, lat_vector_xyz
  real(kind=r_def), dimension(ndf_w2b,3)   :: rad_vector_xyz

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! For spherical geometry, need to rotate from (lon,lat,r) components
  ! For planar geometry, components should already be in (X,Y,Z) coordinates
  select case ( geometry )
  case ( geometry_planar )

    X_vector = (/ 1.0_r_def, 0.0_r_def, 0.0_r_def /)
    Y_vector = (/ 0.0_r_def, 1.0_r_def, 0.0_r_def /)
    Z_vector = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

    do k = 0, nlayers-1
      ! Get index of this cell
      cell_3d = k + 1 + (col-1)*nlayers
      u_lon_op(cell_3d,:,:) = 0.0_r_def
      u_lat_op(cell_3d,:,:) = 0.0_r_def
      u_rad_op(cell_3d,:,:) = 0.0_r_def

      ! Compute Jacobian for this cell
      do df_chi = 1, ndf_chi
        chi1_e(df_chi) = chi1(map_chi(df_chi) + k)
        chi2_e(df_chi) = chi2(map_chi(df_chi) + k)
        chi3_e(df_chi) = chi3(map_chi(df_chi) + k)
      end do

      call coordinate_jacobian(coord_system, geometry, topology, scaled_radius, &
                               ndf_chi, ndf_w2b, chi1_e, chi2_e, chi3_e,        &
                               ipanel, chi_basis, chi_diff_basis, jacobian, dj)
      call coordinate_jacobian_inverse(ndf_w2b, jacobian, dj, jac_inv)

      ! X and Y components contribute equally to all W2 DoFs
      do df_w2 = 1, ndf_w2b
        u_lon_op(cell_3d,df_w2,1) = dj(df_w2)* &
          dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), X_vector(:)))
        u_lat_op(cell_3d,df_w2,1) = dj(df_w2)* &
          dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), Y_vector(:)))
      end do

      ! Z component is in Wtheta -- assume each Wtheta DoF contributes equally
      ! to the horizontal W2 DoFs (so add an extra half factor here), but for
      ! vertical W2 DoFs only contribute once
      do df_w2 = 1, 4
        do df_wt = 1, ndf_wt
          u_rad_op(cell_3d,df_w2,df_wt) = 0.5_r_def*dj(df_w2)* &
            dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), Z_vector(:)))
        end do
      end do

      u_rad_op(cell_3d,B,1) = dj(B)* &
        dot_product(face_normals(:,B), matmul(jac_inv(:,:,B), Z_vector(:)))
      u_rad_op(cell_3d,T,2) = dj(T)* &
        dot_product(face_normals(:,T), matmul(jac_inv(:,:,T), Z_vector(:)))

    end do

  case ( geometry_spherical )

    lon_vector_llr = (/ 1.0_r_def, 0.0_r_def, 0.0_r_def /)
    lat_vector_llr = (/ 0.0_r_def, 1.0_r_def, 0.0_r_def /)
    rad_vector_llr = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

    do k = 0, nlayers-1
      ! Get index of this cell
      cell_3d = k + 1 + (col-1)*nlayers
      u_lon_op(cell_3d,:,:) = 0.0_r_def
      u_lat_op(cell_3d,:,:) = 0.0_r_def
      u_rad_op(cell_3d,:,:) = 0.0_r_def

      ! Compute Jacobian for this cell
      do df_chi = 1, ndf_chi
        chi1_e(df_chi) = chi1(map_chi(df_chi) + k)
        chi2_e(df_chi) = chi2(map_chi(df_chi) + k)
        chi3_e(df_chi) = chi3(map_chi(df_chi) + k)
      end do

      call coordinate_jacobian(coord_system, geometry, topology, scaled_radius, &
                               ndf_chi, ndf_w2b, chi1_e, chi2_e, chi3_e,        &
                               ipanel, chi_basis, chi_diff_basis, jacobian, dj)
      call coordinate_jacobian_inverse(ndf_w2b, jacobian, dj, jac_inv)

      ! Convert (lon,lat,r) vectors into (X,Y,Z) components
      chi1_w2(:) = 0.0_r_def
      chi2_w2(:) = 0.0_r_def
      chi3_w2(:) = 0.0_r_def

      do df_w2 = 1, ndf_w2b
        ! Get chi at this W2 DoF
        do df_chi = 1, ndf_chi
          chi1_w2(df_w2) = chi1_w2(df_w2) + chi_basis(1,df_chi,df_w2)*chi1_e(df_chi)
          chi2_w2(df_w2) = chi2_w2(df_w2) + chi_basis(1,df_chi,df_w2)*chi2_e(df_chi)
          chi3_w2(df_w2) = chi3_w2(df_w2) + chi_basis(1,df_chi,df_w2)*chi3_e(df_chi)
        end do

        ! Calculate (lon,lat,r) coordinates for W2 points in this cell
        call chi2llr(chi1_w2(df_w2), chi2_w2(df_w2), chi3_w2(df_w2), ipanel, &
                     llr(1), llr(2), llr(3))

        ! Rotate (lon,lat,r) unit vectors to (X,Y,Z) coordinates
        lon_vector_xyz(df_w2,:) = sphere2cart_vector(lon_vector_llr, llr)
        lat_vector_xyz(df_w2,:) = sphere2cart_vector(lat_vector_llr, llr)
        rad_vector_xyz(df_w2,:) = sphere2cart_vector(rad_vector_llr, llr)

      end do

      ! Lon and lat components contribute equally to all W2 DoFs
      do df_w2 = 1, ndf_w2b
        u_lon_op(cell_3d,df_w2,1) = dj(df_w2)* &
          dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), lon_vector_xyz(df_w2,:)))
        u_lat_op(cell_3d,df_w2,1) = dj(df_w2)* &
          dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), lat_vector_xyz(df_w2,:)))
      end do

      ! Radial component is in Wtheta -- assume each Wtheta DoF contributes equally
      ! to the horizontal W2 DoFs (so add an extra half factor here), but for
      ! vertical W2 DoFs only contribute once
      do df_w2 = 1, 4
        do df_wt = 1, ndf_wt
          u_rad_op(cell_3d,df_w2,df_wt) = 0.5_r_def*dj(df_w2)* &
            dot_product(face_normals(:,df_w2), matmul(jac_inv(:,:,df_w2), rad_vector_xyz(df_w2,:)))
        end do
      end do

      u_rad_op(cell_3d,B,1) = dj(B)* &
        dot_product(face_normals(:,B), matmul(jac_inv(:,:,B), rad_vector_xyz(B,:)))
      u_rad_op(cell_3d,T,2) = dj(T)* &
        dot_product(face_normals(:,T), matmul(jac_inv(:,:,T), rad_vector_xyz(T,:)))

    end do

  end select

  ! Enforce boundary condition at bottom and top
  cell_3d = 1 + (col-1)*nlayers
  u_lon_op(cell_3d,B,:) = 0.0_r_def
  u_lat_op(cell_3d,B,:) = 0.0_r_def
  u_rad_op(cell_3d,B,:) = 0.0_r_def

  cell_3d = nlayers + (col-1)*nlayers
  u_lon_op(cell_3d,T,:) = 0.0_r_def
  u_lat_op(cell_3d,T,:) = 0.0_r_def
  u_rad_op(cell_3d,T,:) = 0.0_r_def

end subroutine compute_sample_u_ops_code

end module sci_compute_sample_u_ops_kernel_mod
