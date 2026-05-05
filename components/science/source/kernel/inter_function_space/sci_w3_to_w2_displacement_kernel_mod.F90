!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the effective horizontal displacement corresponding to the
!!        error when averaging a W3 to W2 points.
!> @details Uses the coordinate fields to compute the displacement between a
!!          W2 point and the effective averaging point when averaging a scalar
!!          field from W3 to W2. Only intended to be used on the cubed-sphere.
!!          This kernel is only designed for lowest order finite elements.
module sci_w3_to_w2_displacement_kernel_mod

  use argument_mod,          only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_BASIS, GH_EVALUATOR,    &
                                    CELL_COLUMN, GH_SCALAR,    &
                                    GH_LOGICAL
  use fs_continuity_mod,     only : W3, W2H, Wchi
  use constants_mod,         only : r_def, i_def
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : E, W, N, S

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: w3_to_w2_displacement_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,   W2H),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3)                         &
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(Wchi, GH_BASIS)                                           &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
    integer :: gh_evaluator_targets(2) = (/ W2H, W3 /)
  contains
    procedure, nopass :: w3_to_w2_displacement_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: w3_to_w2_displacement_code

  contains

  !> @brief Calculates the effective horizontal displacement corresponding to
  !!        the error when averaging a W3 to W2 points
  !> @param[in]     nlayers       Number of layers
  !> @param[in,out] displacement  2D W2H field containing the displacements
  !!                              corresponding to the averaging error. This is
  !!                              dimensionless, being divided by the cell width
  !> @param[in]     chi_1         The first coordinate field
  !> @param[in]     chi_2         The second coordinate field
  !> @param[in]     chi_3         The third coordinate field
  !> @param[in]     panel_id      ID for panels of the underlying mesh
  !> @param[in]     dummy_w3      An unused dummy field in W3
  !> @param[in]     ndf_w2h       Number of DoFs for W2H per cell
  !> @param[in]     undf_w2h      Number of unique DoFs for W2H per partition
  !> @param[in]     map_w2h       The DoF map for bottom layer cells for W2H
  !> @param[in]     ndf_chi       Number of DoFs for Wchi per cell
  !> @param[in]     undf_chi      Number of unique DoFs for Wchi per partition
  !> @param[in]     map_chi       The DoF map for bottom layer cells for Wchi
  !> @param[in]     basis_chi_w2h Wchi basis functions evaluated at W2H points
  !> @param[in]     basis_chi_w3  Wchi basis functions evaluated at W3 points
  !> @param[in]     ndf_pid       Number of DoFs for panel id per cell
  !> @param[in]     undf_pid      Number of unique DoFs for panel id per partition
  !> @param[in]     map_pid       The DoF map for bottom layer cells for panel ID
  !> @param[in]     ndf_w3        Number of DoFs for W3 per cell
  !> @param[in]     undf_w3       Number of unique DoFs for W3 per partition
  !> @param[in]     map_w3        The DoF map for bottom layer cells for W3
  subroutine w3_to_w2_displacement_code( nlayers,       &
                                         displacement,  &
                                         chi_1,         &
                                         chi_2,         &
                                         chi_3,         &
                                         panel_id,      &
                                         dummy_w3,      &
                                         ndf_w2h,       &
                                         undf_w2h,      &
                                         map_w2h,       &
                                         ndf_chi,       &
                                         undf_chi,      &
                                         map_chi,       &
                                         basis_chi_w2h, &
                                         basis_chi_w3,  &
                                         ndf_pid,       &
                                         undf_pid,      &
                                         map_pid,       &
                                         ndf_w3,        &
                                         undf_w3,       &
                                         map_w3 )

    use sci_chi_transform_mod, only: chi2abr

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w2h, undf_w2h
    integer(kind=i_def), intent(in)    :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

    real(kind=r_def),    intent(inout) :: displacement(undf_w2h)
    real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
    real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
    real(kind=r_def),    intent(in)    :: chi_3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    real(kind=r_def),    intent(in)    :: dummy_w3(undf_w3)
    real(kind=r_def),    intent(in)    :: basis_chi_w2h(1,ndf_chi,ndf_w2h)
    real(kind=r_def),    intent(in)    :: basis_chi_w3(1,ndf_chi,ndf_w3)

    ! Vertical cell index
    integer(kind=i_def) :: df_w2h, df_w3, df_chi
    integer(kind=i_def) :: ipanel
    real(kind=r_def)    :: cell_width_opposite, cell_half_width_adjacent
    real(kind=r_def)    :: alpha_w3, beta_w3, dummy_r
    real(kind=r_def)    :: alpha_w2h(4), beta_w2h(4)
    real(kind=r_def)    :: chi1_at_dof, chi2_at_dof, chi3_at_dof
    real(kind=r_def)    :: e_alpha(3), e_beta(3)
    real(kind=r_def)    :: phi, varrho

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! The output field is 2D so we can ignore layers

    ! Get alpha and beta values at each DoF
    ! W3 points ----------------------------------------------------------------
    chi1_at_dof = 0.0_r_def
    chi2_at_dof = 0.0_r_def
    chi3_at_dof = 0.0_r_def
    ! Get chi at this point and then transform to alpha/beta coords
    df_w3 = 1
    do df_chi = 1, ndf_chi
      chi1_at_dof = chi1_at_dof + &
        basis_chi_w3(1,df_chi,df_w3) * chi_1(map_chi(df_chi))
      chi2_at_dof = chi2_at_dof + &
        basis_chi_w3(1,df_chi,df_w3) * chi_2(map_chi(df_chi))
      chi3_at_dof = chi3_at_dof + &
        basis_chi_w3(1,df_chi,df_w3) * chi_3(map_chi(df_chi))
    end do
    call chi2abr(chi1_at_dof, chi2_at_dof, chi3_at_dof, ipanel, &
                 alpha_w3, beta_w3, dummy_r)

    ! W2H points ---------------------------------------------------------------
    do df_w2h = 1, 4
      chi1_at_dof = 0.0_r_def
      chi2_at_dof = 0.0_r_def
      chi3_at_dof = 0.0_r_def
      ! Get chi at this point and then transform to alpha/beta coords
      df_w3 = 1
      do df_chi = 1, ndf_chi
        chi1_at_dof = chi1_at_dof + &
          basis_chi_w2h(1,df_chi,df_w2h) * chi_1(map_chi(df_chi))
        chi2_at_dof = chi2_at_dof + &
          basis_chi_w2h(1,df_chi,df_w2h) * chi_2(map_chi(df_chi))
        chi3_at_dof = chi3_at_dof + &
          basis_chi_w2h(1,df_chi,df_w2h) * chi_3(map_chi(df_chi))
      end do
      call chi2abr(chi1_at_dof, chi2_at_dof, chi3_at_dof, ipanel, &
                   alpha_w2h(df_w2h), beta_w2h(df_w2h), dummy_r)

    end do

    ! Compute angle between basis functions ------------------------------------
    varrho = sqrt(1.0_r_def + (tan(alpha_w3))**2.0_r_def + (tan(beta_w3))**2.0_r_def)
    e_alpha(1) = -tan(alpha_w3)*cos(beta_w3)/varrho
    e_alpha(2) = 1.0_r_def/cos(beta_w3)/varrho
    e_alpha(3) = -tan(alpha_w3)*sin(beta_w3)/varrho
    e_beta(1) = -tan(beta_w3)*cos(alpha_w3)/varrho
    e_beta(2) = -tan(beta_w3)*sin(alpha_w3)/varrho
    e_beta(3) = 1.0_r_def/cos(alpha_w3)/varrho
    phi = asin(dot_product(e_alpha, e_beta))

    ! Compute contribution to displacement for each face -----------------------
    do df_w2h = 1, 4
      ! Take alpha / beta depending on the face
      if (df_w2h == N .or. df_w2h == S) then
        cell_half_width_adjacent = beta_w3 - beta_w2h(df_w2h)
        cell_width_opposite = alpha_w2h(E) - alpha_w2h(W)
      else
        cell_half_width_adjacent = alpha_w3 - alpha_w2h(df_w2h)
        cell_width_opposite = beta_w2h(N) - beta_w2h(S)
      end if

      ! Half-factor for each side of the face -- could be rmultiplicity but
      ! as this is on the cubed-sphere we can just take it to be 0.5
      displacement(map_w2h(df_w2h)) = displacement(map_w2h(df_w2h)) + &
        0.5_r_def * cell_half_width_adjacent * sin(phi) / cell_width_opposite
    end do

  end subroutine w3_to_w2_displacement_code

end module sci_w3_to_w2_displacement_kernel_mod
