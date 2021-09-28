!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Computes the values of Det(J) at vertical W2 locations.
!> @details This kernel computes the Det(J) at vertical W2 locations. It uses
!>          either the cell above or the cell below for the calculation, based
!>          on the input parameter 'direction'. The Det(J) from a cell is given to
!>          the TOP dof if direction = 0 (below) and the BOTTOM dof if
!>          direction = 1 (above).
!>          This kernel is only designed for the lowest order finite-element spaces.
!>
module calc_directional_detj_at_w2_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_REAL, GH_INC, &
                                GH_SCALAR, GH_INTEGER,     &
                                GH_READ, ANY_SPACE_1,      &
                                GH_DIFF_BASIS, GH_BASIS,   &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                CELL_COLUMN, GH_EVALUATOR

  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: calc_directional_detj_at_w2_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD,   GH_REAL,    GH_INC,  W2),                        &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ, ANY_SPACE_1),               &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                             &
         /)
    type(func_type) :: meta_funcs(1) = (/                                      &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                       &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: calc_directional_detj_at_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_directional_detj_at_w2_code

contains

!> @brief Computes the values of Det(J) at vertical W2 locations using either
!!        the cell above or below.
!> @param[in]     nlayers        The number of layers
!> @param[in,out] detj_w2        The output field containing the Det(J) values at W2 locations
!> @param[in]     chi1           1st (spherical) coordinate field in Wchi
!> @param[in]     chi2           2nd (spherical) coordinate field in Wchi
!> @param[in]     chi3           3rd (spherical) coordinate field in Wchi
!> @param[in]     panel_id       Field giving the ID for mesh panels
!> @param[in]     direction      Parameter specifying cell above (1) or below (0) vertical W2 location
!> @param[in]     ndf_w2         The number of degrees of freedom per cell for the output field
!> @param[in]     undf_w2        The number of unique degrees of freedom for the output field
!> @param[in]     map_w2         Array holding the dofmap for the cell at the base
!!                               of the column for the output field
!> @param[in]     ndf_chi        The number of degrees of freedom per cell for the coordinate field
!> @param[in]     undf_chi       The number of unique degrees of freedom for the coordinate field
!> @param[in]     map_chi        Array holding the dofmap for the cell at the base
!!                               of the column for the coordinate field
!> @param[in]     basis_chi      Basis functions of the coordinate space evaluated at the nodal points
!> @param[in]     diff_basis_chi The differential basis functions of the coordinate
!!                               space evaluated at the nodal points
!> @param[in]     ndf_pid        Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid       Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid        Dofmap for the cell at the base of the column for panel_id
!>
subroutine calc_directional_detj_at_w2_code( nlayers,                    &
                                             detj_w2,                    &
                                             chi1, chi2, chi3,           &
                                             panel_id,                   &
                                             direction,                  &
                                             ndf_w2, undf_w2, map_w2,    &
                                             ndf_chi, undf_chi, map_chi, &
                                             basis_chi, diff_basis_chi,  &
                                             ndf_pid, undf_pid, map_pid  &
                                            )

  use coordinate_jacobian_mod, only: pointwise_coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in)    :: nlayers
  integer(kind=i_def),                            intent(in)    :: ndf_w2
  integer(kind=i_def),                            intent(in)    :: undf_w2
  integer(kind=i_def),                            intent(in)    :: ndf_chi
  integer(kind=i_def),                            intent(in)    :: undf_chi
  integer(kind=i_def),                            intent(in)    :: ndf_pid
  integer(kind=i_def),                            intent(in)    :: undf_pid
  real(kind=r_def), dimension(undf_w2),           intent(inout) :: detj_w2
  real(kind=r_def), dimension(undf_chi),          intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),          intent(in)    :: panel_id
  integer(kind=i_def),                            intent(in)    :: direction
  integer(kind=i_def), dimension(ndf_w2),         intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_chi),        intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),        intent(in)    :: map_pid
  real(kind=r_def), dimension(3,ndf_chi,ndf_w2),  intent(in)    :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_chi,ndf_w2),  intent(in)    :: basis_chi

  ! Internal variables
  integer(kind=i_def)                  :: df, cdf, k, k_start, k_end
  integer(kind=i_def)                  :: ipanel
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                     :: detj
  real(kind=r_def), dimension(3,3)     :: jacobian

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! If direction = 0 we use the cell below and set
  ! df = 6, k_start = 0, k_end = nlayers-2
  ! If direction = 1 we use the cell above and set
  ! df = 5, k_start = 1, k_end = nlayers-1

  df      = 6_i_def - direction
  k_start = direction
  k_end   = nlayers - 2_i_def + direction

  do k = k_start, k_end

    do cdf = 1,ndf_chi
      chi1_e(cdf) = chi1(map_chi(cdf) + k)
      chi2_e(cdf) = chi2(map_chi(cdf) + k)
      chi3_e(cdf) = chi3(map_chi(cdf) + k)
    end do

    call pointwise_coordinate_jacobian(ndf_chi, chi1_e, chi2_e, chi3_e, &
                                       ipanel, basis_chi(:,:,df),       &
                                       diff_basis_chi(:,:,df),          &
                                       jacobian, detj)

    detj_w2(map_w2(df)+k) = detj

  end do

end subroutine calc_directional_detj_at_w2_code

end module calc_directional_detj_at_w2_kernel_mod
