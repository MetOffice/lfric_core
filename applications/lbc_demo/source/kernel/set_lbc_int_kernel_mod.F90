!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module set_lbc_int_kernel_mod

  use argument_mod,  only: arg_type, func_type,             &
                           GH_FIELD, GH_SCALAR,             &
                           GH_REAL, GH_INTEGER, GH_LOGICAL, &
                           GH_READ, GH_WRITE,               &
                           ANY_SPACE_8, ANY_SPACE_9,        &
                           GH_BASIS, CELL_COLUMN, GH_EVALUATOR
  use constants_mod, only: r_def, i_def, l_def, radians_to_degrees
  use kernel_mod,    only: kernel_type

  use base_mesh_config_mod, only: geometry_spherical

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata
  !> needed by the Psy layer
  type, public, extends(kernel_type) :: set_lbc_int_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                           &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),               & ! ndata
         arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ),               & ! ndata_first
         arg_type(GH_FIELD,   GH_INTEGER, GH_WRITE, ANY_SPACE_8), & ! lbc
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),               & ! geometry
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9)  & ! chi_1,chi_2,chi_3
         /)
    type(func_type) :: meta_funcs(1) = (/ &
         func_type(ANY_SPACE_9, GH_BASIS) &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: set_lbc_int_code
  end type set_lbc_int_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_lbc_int_code

contains

!!----------------------------------------------------------------------------
!! @brief         Assigns an real value to dofs based on their
!!                x-lat location.
!! @description   This kernel is to initialise in memory an lbc real field
!!                to check that the correct orientation of LBC field data is
!!                to a LAM field. DoFs are assign reals based on their x-lat
!!                location with repect to (0,0)
!!
!!                              4    |    1
!!                               8   |   5
!!                                   |
!!                            -----(0,0)-----
!!                                   |
!!                               7   |   6
!!                              3    |    2
!!
!! @param[in]      nlayers     Number of layers
!! @param[in]      ndata       Number of data values at each DoF location
!! @param[in]      ndata_first True if data at DoFs are contiguous in memory
!! @param[in,out]  lbc         LBC integer field data
!! @param[in]      geometry    Geometry enumeration value
!! @param[in]      chi_1       1st physical coordinate field
!! @param[in]      chi_2       2nd physical coordinate field
!! @param[in]      chi_3       3rd physical coordinate field
!! @param[in]      ndf         Number of degrees of freedom per cell
!! @param[in]      undf        Number of unique degrees of freedom
!! @param[in]      map         Dofmap for the cell at the base of the column
!! @param[in]      ndf_chi     Number of degrees of freedom per cell for chi
!! @param[in]      undf_chi    Number of unique degrees of freedom for chi
!! @param[in]      map_chi     Dofmap for the cell at the base of the column for chi
!! @param[in]      chi_basis   Chi basis functions evaluated at w3 nodes
subroutine set_lbc_int_code( nlayers, ndata, ndata_first, lbc, &
                             geometry, chi_1, chi_2, chi_3,    &
                             ndf, undf, map,                   &
                             ndf_chi, undf_chi, map_chi,       &
                             chi_basis )

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndata
  integer(i_def), intent(in) :: geometry
  integer(i_def), intent(in) :: ndf
  integer(i_def), intent(in) :: undf
  integer(i_def), intent(in) :: ndf_chi
  integer(i_def), intent(in) :: undf_chi

  logical(l_def), intent(in) :: ndata_first

  integer(i_def), dimension(ndf),     intent(in)    :: map
  integer(i_def), dimension(ndf_chi), intent(in)    :: map_chi
  integer(i_def), dimension(undf),    intent(inout) :: lbc

  real(r_def), dimension(undf_chi),   intent(in)    :: chi_1
  real(r_def), dimension(undf_chi),   intent(in)    :: chi_2
  real(r_def), dimension(undf_chi),   intent(in)    :: chi_3

  real(r_def), dimension(1, ndf_chi, ndf), intent(in) :: chi_basis

  ! Internal variables
  integer(i_def) :: df, dfc, k
  integer(i_def) :: offset, idata
  real(r_def)    :: coord(3)
  real(r_def)    :: x, y
  real(r_def)    :: radius

  real(r_def), parameter :: threshold_value = 10.0_r_def

  integer(i_def) :: lbc_value

  real(r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e

  ! Note: (x,y) does not imply cartesian coords.

  do k=0, nlayers-1

    do dfc=1, ndf_chi
      chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
      chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
      chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
    end do

    do df=1, ndf

      x = 0.0_r_def
      y = 0.0_r_def

      coord(:) = 0.0_r_def
      do dfc=1, ndf_chi
        coord(1) = coord(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
        coord(2) = coord(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
        coord(3) = coord(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
      end do

      if (geometry == geometry_spherical) then
        ! Convert coords to degrees
        x = coord(1)*radians_to_degrees
        y = coord(2)*radians_to_degrees
      else
        x = coord(1)
        y = coord(2)
      end if

      radius = SQRT(x**2.0_r_def + y**2.0_r_def)

      ! Should give (clockwise):
      !   1, 2, 3, 4 in outer quadrants
      !   5, 6, 7, 8 in inner quadrants
      ! around (0,0) starting in +x,+y quadrant
      if (x > 0.0_r_def) then

        if (y > 0.0_r_def) then

          ! (+ve x, +ve y)
          if ( radius >= threshold_value ) then
            lbc_value = 1_i_def
          else
            lbc_value = 5_i_def
          end if

        else

          ! (+ve x, -ve y)
          if ( radius >= threshold_value ) then
            lbc_value = 2_i_def
          else
            lbc_value = 6_i_def
          end if

        end if ! test on y

      else

        if (y > 0.0_r_def) then

          ! (-ve x, +ve y)
          if ( radius >= threshold_value ) then
            lbc_value = 4_i_def
          else
            lbc_value = 8_i_def
          end if

        else

          ! (-ve x, -ve y)
          if ( radius >= threshold_value ) then
            lbc_value = 3_i_def
          else
            lbc_value = 7_i_def
          end if

        end if ! test on y

      end if ! test on x

      ! Allow for multidata fields
      do idata=0, ndata-1
        if (ndata_first) then
          offset = k*ndata + idata
        else
          offset = k + idata*nlayers
        end if
        lbc(map(df) + offset) = lbc_value
      end do

    end do ! df
  end do ! k

end subroutine set_lbc_int_code

end module set_lbc_int_kernel_mod
