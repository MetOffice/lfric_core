!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Returns latitude and longitude fields
!>
module sci_compute_latlon_kernel_mod

  use argument_mod,         only: arg_type, func_type,       &
                                  GH_FIELD, GH_REAL,         &
                                  GH_WRITE, GH_READ,         &
                                  ANY_SPACE_1, &
                                  ANY_DISCONTINUOUS_SPACE_3, &
                                  ANY_SPACE_9, GH_BASIS,     &
                                  CELL_COLUMN, GH_EVALUATOR
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use sci_chi_transform_mod, only: chi2llr

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: compute_latlon_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1), &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS)                                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: compute_latlon_code
  end type


  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_latlon_code

contains

!> @brief Calculates the latitude and longitude fields from the x, y and z components
!> @details Will only work at lowest order for now
!> @param[in]     nlayers   The number of layers (always 1)
!> @param[in,out] latitude  Latitude field data
!> @param[in,out] longitude Longitude field data
!> @param[in]     chi_1     First component of the coordinate field
!> @param[in]     chi_2     Second component of the coordinate field
!> @param[in]     chi_3     Third component of the coordinate field
!> @param[in]     panel_id  A field giving the ID for mesh panels
!> @param[in]     ndf_x     Number of degrees of freedom per cell for height
!> @param[in]     undf_x    Number of unique degrees of freedom for height
!> @param[in]     map_x     Dofmap for the cell at the base of the column for height
!> @param[in]     ndf_chi   The number of degrees of freedom per cell for chi
!> @param[in]     undf_chi  The number of unique degrees of freedom for chi
!> @param[in]     map_chi   Dofmap for the cell at the base of the column for chi
!> @param[in]     basis_chi Basis functions evaluated at nodal points for height
!> @param[in]     ndf_pid   Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid  Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid   Dofmap for the cell at the base of the column for panel_id
subroutine compute_latlon_code(nlayers,                         &
                               latitude, longitude,             &
                               chi_1, chi_2, chi_3,             &
                               panel_id,                        &
                               ndf_x, undf_x, map_x,            &
                               ndf_chi, undf_chi, map_chi,      &
                               basis_chi,                       &
                               ndf_pid, undf_pid, map_pid       &
                               )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_x, undf_x
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid

  real(kind=r_def), dimension(undf_x), intent(inout) :: latitude, longitude
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)  :: panel_id

  integer(kind=i_def), dimension(ndf_x), intent(in)          :: map_x
  integer(kind=i_def), dimension(ndf_chi), intent(in)        :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in)        :: map_pid
  real(kind=r_def), dimension(1, ndf_chi, ndf_x), intent(in) :: basis_chi

  ! Internal variables
  integer(kind=i_def) :: df_chi, df_x, k, ipanel
  real(kind=r_def)    :: coords(3), lat, lon, radius

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df_x = 1, ndf_x
      coords(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        coords(1) = coords(1) + chi_1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        coords(2) = coords(2) + chi_2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        coords(3) = coords(3) + chi_3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do
      call chi2llr(coords(1), coords(2), coords(3), ipanel, lon, lat, radius)
      latitude(map_x(df_x) + k) = lat
      longitude(map_x(df_x) + k) = lon
    end do
  end do

end subroutine compute_latlon_code

end module sci_compute_latlon_kernel_mod
