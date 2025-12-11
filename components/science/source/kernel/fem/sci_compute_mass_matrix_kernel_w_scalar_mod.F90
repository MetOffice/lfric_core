!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which
! you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Provides access to the members of the ws_kernel class.
!> Accessor functions for the ws_kernel class are defined in this module.
!>
!> @param RHS_w_scalar_code        Code to implement the RHS for a w0 or wtheta field
!> @param gaussian_quadrature      Contains result of gaussian quadrature
!>
module sci_compute_mass_matrix_kernel_w_scalar_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_OPERATOR, GH_FIELD,     &
                                     GH_LOGICAL, GH_SCALAR,     &
                                     GH_READ, GH_WRITE,         &
                                     GH_REAL, ANY_SPACE_2,      &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only: i_def, r_single, r_double, l_def
  use sci_coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,       only: W0, Wtheta, Wchi
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w_scalar_type
    private
    type(arg_type) :: meta_args(4) = (/                                          &
         arg_type(GH_OPERATOR, GH_REAL,    GH_WRITE, ANY_SPACE_2, ANY_SPACE_2),  &
         arg_type(GH_FIELD*3,  GH_REAL,    GH_READ,  Wchi),                      &
         arg_type(GH_FIELD,    GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,   GH_LOGICAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(ANY_SPACE_2, GH_BASIS),                                   &
         func_type(Wchi,        GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  end type compute_mass_matrix_kernel_w_scalar_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_mass_matrix_w_scalar_code

  ! Generic interface for real32 and real64 types
  interface compute_mass_matrix_w_scalar_code
     module procedure  compute_mass_matrix_w_scalar_code_r32,    &
                       compute_mass_matrix_w_scalar_code_r32r64, &
                       compute_mass_matrix_w_scalar_code_r64
  end interface

contains

  !> @brief This subroutine computes the mass matrix for scalar spaces w0 and wtheta
  !! @param[in] cell     The cell number
  !! @param[in] nlayers  The number of layers.
  !! @param[in] ncell_3d ncell*ndf
  !! @param[in,out] mm    The mass matrix data array
  !! @param[in] chi1     1st coordinate field in Wchi
  !! @param[in] chi2     2nd coordinate field in Wchi
  !! @param[in] chi3     3rd coordinate field in Wchi
  !! @param[in] panel_id Field giving the ID for mesh panels.
  !! @param[in] extended_mesh Compute on an extended mesh
  !! @param[in] ndf_w_scalar   The number of degrees of freedom per cell for w_scalar.
  !! @param[in] ndf_chi  The number of degrees of freedom per cell for chi.
  !! @param[in] basis_w_scalar 4-dim array holding SCALAR basis functions evaluated at
  !! quadrature points for w_scalar.
  !! @param[in] undf_chi Number of unique degrees of freedom for chi
  !! @param[in] map_chi  Array holding the dofmap for the cell at the base of the
  !! column for chi
  !! @param[in] basis_chi 4-dim array holding basis functions evaluated at quadrature
  !! points for chi
  !! @param[in] diff_basis_chi 4-dim array holding VECTOR differential basis functions
  !! evaluated at quadrature points for chi
  !! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  !! @param[in] nqp_h    Number of horizontal quadrature points
  !! @param[in] nqp_v    Number of vertical quadrature points
  !! @param[in] wqp_h    Quadrature weights horizontal
  !! @param[in] wqp_v    Quadrature weights vertical

  ! R_SINGLE PRECISION
  ! ==================
  subroutine compute_mass_matrix_w_scalar_code_r32(             &
                                 cell, nlayers, ncell_3d, mm,   &
                                 chi1, chi2, chi3, panel_id,    &
                                 extended_mesh,                 &
                                 ndf_w_scalar, basis_w_scalar,  &
                                 ndf_chi, undf_chi, map_chi,    &
                                 basis_chi, diff_basis_chi,     &
                                 ndf_pid, undf_pid, map_pid,    &
                                 nqp_h, nqp_v, wqp_h, wqp_v  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in) :: nlayers, ndf_w_scalar
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in) :: ncell_3d

    integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

    logical(kind=l_def), intent(in) :: extended_mesh

    real(kind=r_single), dimension(ncell_3d,ndf_w_scalar,ndf_w_scalar), &
                                          intent(inout) :: mm

    real(kind=r_single), dimension(1,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: basis_chi
    real(kind=r_single), dimension(3,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: diff_basis_chi
    real(kind=r_single), dimension(1,ndf_w_scalar,nqp_h,nqp_v), &
                                             intent(in) :: basis_w_scalar

    real(kind=r_single), dimension(undf_chi), intent(in) :: chi1
    real(kind=r_single), dimension(undf_chi), intent(in) :: chi2
    real(kind=r_single), dimension(undf_chi), intent(in) :: chi3
    real(kind=r_single), dimension(undf_pid), intent(in) :: panel_id

    real(kind=r_single), dimension(nqp_h),    intent(in) :: wqp_h
    real(kind=r_single), dimension(nqp_v),    intent(in) :: wqp_v

    ! Internal variables
    integer(kind=i_def)                        :: df, df2, k, ik
    integer(kind=i_def)                        :: ipanel, qp1, qp2
    real(kind=r_single), dimension(ndf_chi)         :: chi1_e
    real(kind=r_single), dimension(ndf_chi)         :: chi2_e
    real(kind=r_single), dimension(ndf_chi)         :: chi3_e
    real(kind=r_single)                             :: integrand
    real(kind=r_single), dimension(nqp_h,nqp_v)     :: dj
    real(kind=r_single), dimension(3,3,nqp_h,nqp_v) :: jac

    if ( extended_mesh ) then
      ipanel = int(panel_id(1), i_def)
    else
      ipanel = int(panel_id(map_pid(1)), i_def)
    end if

    ! loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell-1)*nlayers

      ! indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do

      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,          &
                               chi1_e, chi2_e, chi3_e, ipanel, &
                               basis_chi, diff_basis_chi,      &
                               jac, dj)

      do df2 = 1, ndf_w_scalar
        do df = df2, ndf_w_scalar ! mass matrix is symmetric
          mm(ik,df,df2) = 0.0_r_single
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              integrand = wqp_h(qp1) * wqp_v(qp2) * &
                   basis_w_scalar(1,df,qp1,qp2)   * &
                   basis_w_scalar(1,df2,qp1,qp2)  * &
                   dj(qp1,qp2)
              mm(ik,df,df2) = mm(ik,df,df2) + integrand
            end do
          end do
        end do
        do df = df2, 1, -1
          mm(ik,df,df2) = mm(ik,df2,df)
        end do
      end do

    end do ! end of k loop

  end subroutine compute_mass_matrix_w_scalar_code_r32

  ! MIXED PRECISION
  ! ===============
  subroutine compute_mass_matrix_w_scalar_code_r32r64(          &
                                 cell, nlayers, ncell_3d, mm,   &
                                 chi1, chi2, chi3, panel_id,    &
                                 extended_mesh,                 &
                                 ndf_w_scalar, basis_w_scalar,  &
                                 ndf_chi, undf_chi, map_chi,    &
                                 basis_chi, diff_basis_chi,     &
                                 ndf_pid, undf_pid, map_pid,    &
                                 nqp_h, nqp_v, wqp_h, wqp_v  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in) :: nlayers, ndf_w_scalar
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in) :: ncell_3d

    integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

    logical(kind=l_def), intent(in) :: extended_mesh

    real(kind=r_single), dimension(ncell_3d,ndf_w_scalar,ndf_w_scalar), &
                                          intent(inout) :: mm

    real(kind=r_double), dimension(1,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: basis_chi
    real(kind=r_double), dimension(3,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: diff_basis_chi
    real(kind=r_double), dimension(1,ndf_w_scalar,nqp_h,nqp_v), &
                                             intent(in) :: basis_w_scalar

    real(kind=r_double), dimension(undf_chi), intent(in) :: chi1, chi2, chi3
    real(kind=r_double), dimension(undf_pid), intent(in) :: panel_id

    real(kind=r_double), dimension(nqp_h),    intent(in) :: wqp_h
    real(kind=r_double), dimension(nqp_v),    intent(in) :: wqp_v

    ! Internal variables
    integer(kind=i_def)                        :: df, df2, k, ik
    integer(kind=i_def)                        :: ipanel, qp1, qp2
    real(kind=r_double), dimension(ndf_chi)         :: chi1_e
    real(kind=r_double), dimension(ndf_chi)         :: chi2_e
    real(kind=r_double), dimension(ndf_chi)         :: chi3_e
    real(kind=r_double)                             :: integrand
    real(kind=r_double), dimension(nqp_h,nqp_v)     :: dj
    real(kind=r_double), dimension(3,3,nqp_h,nqp_v) :: jac

    if ( extended_mesh ) then
      ipanel = int(panel_id(1), i_def)
    else
      ipanel = int(panel_id(map_pid(1)), i_def)
    end if

    ! loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell-1)*nlayers

      ! indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do

      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,          &
                               chi1_e, chi2_e, chi3_e, ipanel, &
                               basis_chi, diff_basis_chi,      &
                               jac, dj)

      do df2 = 1, ndf_w_scalar
        do df = df2, ndf_w_scalar ! mass matrix is symmetric
          mm(ik,df,df2) = 0.0_r_single
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              integrand = wqp_h(qp1) * wqp_v(qp2) * &
                   basis_w_scalar(1,df,qp1,qp2)   * &
                   basis_w_scalar(1,df2,qp1,qp2)  * &
                   dj(qp1,qp2)
              mm(ik,df,df2) = mm(ik,df,df2) + real(integrand, kind=r_single)
            end do
          end do
        end do
        do df = df2, 1, -1
          mm(ik,df,df2) = mm(ik,df2,df)
        end do
      end do

    end do ! end of k loop

  end subroutine compute_mass_matrix_w_scalar_code_r32r64


  ! R_DOUBLE PRECISION
  ! ==================
  subroutine compute_mass_matrix_w_scalar_code_r64(             &
                                 cell, nlayers, ncell_3d, mm,   &
                                 chi1, chi2, chi3, panel_id,    &
                                 extended_mesh,                 &
                                 ndf_w_scalar, basis_w_scalar,  &
                                 ndf_chi, undf_chi, map_chi,    &
                                 basis_chi, diff_basis_chi,     &
                                 ndf_pid, undf_pid, map_pid,    &
                                 nqp_h, nqp_v, wqp_h, wqp_v  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in) :: nlayers, ndf_w_scalar
    integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
    integer(kind=i_def), intent(in) :: ncell_3d

    integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

    logical(kind=l_def), intent(in) :: extended_mesh

    real(kind=r_double), dimension(ncell_3d,ndf_w_scalar,ndf_w_scalar), &
                                          intent(inout) :: mm

    real(kind=r_double), dimension(1,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: basis_chi
    real(kind=r_double), dimension(3,ndf_chi,nqp_h,nqp_v), &
                                             intent(in) :: diff_basis_chi
    real(kind=r_double), dimension(1,ndf_w_scalar,nqp_h,nqp_v), &
                                             intent(in) :: basis_w_scalar

    real(kind=r_double), dimension(undf_chi), intent(in) :: chi1
    real(kind=r_double), dimension(undf_chi), intent(in) :: chi2
    real(kind=r_double), dimension(undf_chi), intent(in) :: chi3
    real(kind=r_double), dimension(undf_pid), intent(in) :: panel_id

    real(kind=r_double), dimension(nqp_h),    intent(in) :: wqp_h
    real(kind=r_double), dimension(nqp_v),    intent(in) :: wqp_v

    ! Internal variables
    integer(kind=i_def)                        :: df, df2, k, ik
    integer(kind=i_def)                        :: ipanel, qp1, qp2
    real(kind=r_double), dimension(ndf_chi)         :: chi1_e
    real(kind=r_double), dimension(ndf_chi)         :: chi2_e
    real(kind=r_double), dimension(ndf_chi)         :: chi3_e
    real(kind=r_double)                             :: integrand
    real(kind=r_double), dimension(nqp_h,nqp_v)     :: dj
    real(kind=r_double), dimension(3,3,nqp_h,nqp_v) :: jac

    if ( extended_mesh ) then
      ipanel = int(panel_id(1), i_def)
    else
      ipanel = int(panel_id(map_pid(1)), i_def)
    end if

    ! loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell-1)*nlayers

      ! indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do

      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,          &
                               chi1_e, chi2_e, chi3_e, ipanel, &
                               basis_chi, diff_basis_chi,      &
                               jac, dj)

      do df2 = 1, ndf_w_scalar
        do df = df2, ndf_w_scalar ! mass matrix is symmetric
          mm(ik,df,df2) = 0.0_r_double
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              integrand = wqp_h(qp1) * wqp_v(qp2) * &
                   basis_w_scalar(1,df,qp1,qp2)   * &
                   basis_w_scalar(1,df2,qp1,qp2)  * &
                   dj(qp1,qp2)
               mm(ik,df,df2) = mm(ik,df,df2) + integrand
            end do
          end do
        end do
        do df = df2, 1, -1
          mm(ik,df,df2) = mm(ik,df2,df)
        end do
      end do

    end do ! end of k loop

  end subroutine compute_mass_matrix_w_scalar_code_r64

end module sci_compute_mass_matrix_kernel_w_scalar_mod
