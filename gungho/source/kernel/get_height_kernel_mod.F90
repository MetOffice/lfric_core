!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Returns a height field (r or z) from the chi array

!> @details Returns a height field (r or z) from the chi array
module get_height_kernel_mod
  
use kernel_mod,               only: kernel_type
use argument_mod,             only: arg_type, func_type,                 &
                                    GH_FIELD, GH_WRITE, GH_READ, GH_INC, &
                                    W0, W2, W3, ANY_SPACE_9,             &
                                    GH_BASIS,                            &
                                    CELLS, GH_EVALUATOR, EVALUATOR

use constants_mod,            only: r_def
use base_mesh_config_mod,     only: geometry,                            &
                                    base_mesh_geometry_spherical
use planet_config_mod,        only: scaled_radius

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy 
!> layer
type, public, extends(kernel_type) :: get_height_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                   &
       arg_type(GH_FIELD,   GH_WRITE,  W3),             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)      &
       /)
  type(func_type) :: meta_funcs(1) = (/                 &
       func_type(ANY_SPACE_9, GH_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = EVALUATOR
contains
  procedure, nopass :: get_height_code
end type


!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface get_height_kernel_type
   module procedure get_height_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public get_height_code
contains

type(get_height_kernel_type) function get_height_kernel_constructor() result(self)
  return
end function get_height_kernel_constructor

!> @brief Returns a height field (r or z) from the chi array
!>        Will only work at lowest order for now
!! @param[in] nlayers Number of layers
!! @param[inout] height The height field
!! @param[in] chi_1 X component of the coordinate 
!! @param[in] chi_2 Y component of the coordinate 
!! @param[in] chi_3 Z component of the coordinate
!! @param[in] ndf_x Number of degrees of freedom per cell for height
!! @param[in] undf_x Number of unique degrees of freedom for height
!! @param[in] map_x Dofmap for the cell at the base of the column for height
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] basis_chi Basis functions evaluated at nodal points for height
subroutine get_height_code(nlayers,                         &
                           height,                          &
                           chi_1, chi_2, chi_3,             &
                           ndf_x, undf_x, map_x,            &
                           ndf_chi, undf_chi, map_chi,      &
                           basis_chi                        &
                           )
  implicit none

  !Arguments
  integer, intent(in) :: nlayers

  integer, intent(in)                                :: ndf_x, undf_x
  integer, intent(in)                                :: ndf_chi, undf_chi
  real(kind=r_def), dimension(undf_x), intent(inout) :: height
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi_1, chi_2, chi_3

  integer, dimension(ndf_x), intent(in)                       :: map_x
  integer, dimension(ndf_chi), intent(in)                     :: map_chi
  real(kind=r_def), dimension(1, ndf_chi, ndf_x), intent(in)  :: basis_chi

  !Internal variables
  integer          :: df_chi, df_x, k

  real(kind=r_def) :: xyz(3), r
  
  do k = 0, nlayers-1
    do df_x = 1, ndf_x
      xyz(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        xyz(1) = xyz(1) + chi_1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi_2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi_3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do
      if (geometry == base_mesh_geometry_spherical) then
        r = sqrt(xyz(1)**2 + xyz(2)**2 + xyz(3)**2)
        ! NB This will result in the height above 
        ! the spherical representation of the planet
        ! by not necessarily the height above the bottom
        ! of the mesh
        ! This should be reviewed with ticket #562
        r = r - scaled_radius 
      else
        r = xyz(3)
      end if
      height( map_x(df_x) + k ) = r
    end do
  end do

end subroutine get_height_code

end module get_height_kernel_mod
