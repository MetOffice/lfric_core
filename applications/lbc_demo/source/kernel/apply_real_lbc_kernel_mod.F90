!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Updates LAM field with data from an LBC field.
!>
!> @details LAM/LBC fields are given defined as being on two different
!>          function spaces as they are on differing meshes. Although the
!>          Function Space Definition (W0,W1...) and element order will
!>          need to be the same.
!>
!>          The LBC mesh is a mesh object which is a sub-domain of the full
!>          LAM mesh. For every cell in the LBC mesh there is a corresponding
!>          geo-located cell in the LAM mesh.
module apply_real_lbc_kernel_mod

  use argument_mod, only: arg_type,    &
                    GH_FIELD,          &
                    GH_REAL,           &
                    GH_READ, GH_WRITE, &
                    ANY_SPACE_1,       &
                    ANY_SPACE_2,       &
                    CELL_COLUMN


  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: apply_real_lbc_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                        &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_SPACE_1),  & ! lam_field
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_2)   & ! lbc_field
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: apply_real_lbc_kernel_code
  end type apply_real_lbc_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: apply_real_lbc_kernel_code

contains

!> @brief Updates LAM dof data with LBC dof data
!! @param[in]     nlayers      Number of layers
!! @param[in,out] lam_field    The lam field data to be updated
!! @param[in]     lbc_field    The lbc field data to be applied
!! @param[in]     undf_lam     Number of unique degrees of freedom for LAM
!! @param[in]     map_lam      LAM dofmap for cell at base of column
!! @param[in]     ndf          Number of degrees of freedom per cell
!! @param[in]     undf_lbc     Number of unique degrees of freedom for LBC
!! @param[in]     map_lbc      LBC dofmap for cell at base of column
subroutine apply_real_lbc_kernel_code( nlayers,     &
                                       lam_field,   &
                                       lbc_field,   &
                                       undf_lam,    &
                                       map_lam,     &
                                       ndf,         &
                                       ndata,       &
                                       ndata_first, &
                                       undf_lbc,    &
                                       map_lbc )

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf
  integer(i_def), intent(in) :: ndata

  logical, intent(in) :: ndata_first

  integer(i_def), intent(in) :: undf_lam
  integer(i_def), intent(in) :: undf_lbc

  ! Fields
  real(r_def), intent(inout) :: lam_field( undf_lam )
  real(r_def), intent(in)    :: lbc_field( undf_lbc )

  ! Maps
  integer(i_def), intent(in) :: map_lam( ndf )
  integer(i_def), intent(in) :: map_lbc( ndf )

  ! Internal variables
  integer(i_def) :: ilayer, idof, idata, offset

  integer(i_def) :: lbc_index, lam_index

  if (ndata_first) then

    do idof=1, ndf
      do idata=0, ndata-1
        do ilayer=0, nlayers-1

          offset    = ilayer + idata*nlayers

          lbc_index = map_lbc(idof) + offset
          lam_index = map_lam(idof) + offset

          lam_field(lam_index) = lbc_field(lbc_index)

        end do
      end do
    end do

  else

    do idof=1, ndf
      do ilayer=0, nlayers-1
        do idata=0, ndata-1

          offset = idata + ilayer*ndata

          lbc_index = map_lbc(idof) + offset
          lam_index = map_lam(idof) + offset

          lam_field(lam_index) = lbc_field(lbc_index)

        end do
      end do
    end do

  end if

end subroutine apply_real_lbc_kernel_code

end module apply_real_lbc_kernel_mod
