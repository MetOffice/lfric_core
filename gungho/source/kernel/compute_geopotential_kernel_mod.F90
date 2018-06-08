!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the geopotential field.
!>
!> Computes the geopotential field Phi = g*r or g*z for Cartesian domains.
!>
module compute_geopotential_kernel_mod

  use argument_mod,           only : arg_type, func_type,         &
                                     GH_FIELD, GH_READ, GH_INC,   &
                                     ANY_SPACE_9, GH_BASIS,       &
                                     CELLS, GH_EVALUATOR
  use base_mesh_config_mod,   only : geometry, &
                                     base_mesh_geometry_spherical
  use constants_mod,          only : r_def, i_def
  use coord_transform_mod,    only : xyz2llr
  use formulation_config_mod, only : shallow
  use fs_continuity_mod,      only : W0
  use kernel_mod,             only : kernel_type
  use planet_config_mod,      only : gravity, scaled_radius

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: compute_geopotential_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/            &
        arg_type(GH_FIELD,   GH_INC,  W0),         &
        arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9) &
        /)
    type(func_type) :: meta_funcs(1) = (/ &
        func_type(ANY_SPACE_9, GH_BASIS)  &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: compute_geopotential_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface compute_geopotential_kernel_type
    module procedure compute_geopotential_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_geopotential_code

contains

type(compute_geopotential_kernel_type) function compute_geopotential_kernel_constructor() result(self)
  implicit none
  return
end function compute_geopotential_kernel_constructor

!! @param[in] nlayers Number of layers
!! @param[inout] phi Geopotential array
!! @param[in] chi_1 Physical x coordinates
!! @param[in] chi_2 Physical y coordinates
!! @param[in] chi_3 Physical z coordinates
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number of unique degrees of freedom for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Chi basis functions evaluated at w0 nodes
subroutine compute_geopotential_code(nlayers, phi,               &
                                     chi_1, chi_2, chi_3,        &
                                     ndf_w0, undf_w0, map_w0,    &
                                     ndf_chi, undf_chi, map_chi, &
                                     chi_basis)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                           :: nlayers
  integer(kind=i_def), intent(in)                           :: ndf_w0
  integer(kind=i_def), intent(in)                           :: undf_w0
  integer(kind=i_def), intent(in)                           :: ndf_chi
  integer(kind=i_def), intent(in)                           :: undf_chi
  integer(kind=i_def), dimension(ndf_w0), intent(in)        :: map_w0
  integer(kind=i_def), dimension(ndf_chi), intent(in)       :: map_chi
  real(kind=r_def), dimension(undf_w0), intent(inout)       :: phi
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_1
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_2
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi_3
  real(kind=r_def), dimension(1,ndf_chi,ndf_w0), intent(in) :: chi_basis

  ! Internal variables
  integer(kind=i_def) :: df, dfc, k
  real(kind=r_def)    :: x(3)
  real(kind=r_def)    :: lat, lon, r, s

  real(kind=r_def), dimension(ndf_chi)             :: chi_1_e, chi_2_e, chi_3_e


  if ( geometry == base_mesh_geometry_spherical ) then
    if ( shallow ) then
      s = 1.0_r_def
    else
      s = 0.0_r_def
    end if
    do k = 0, nlayers-1
        do dfc = 1, ndf_chi
            chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
            chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
            chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
        end do

        do df = 1, ndf_w0
            x(:) = 0.0_r_def
            do dfc = 1, ndf_chi
                x(1) = x(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
                x(2) = x(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
                x(3) = x(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
            end do
            call xyz2llr(x(1), x(2), x(3), lon, lat, r)

            phi(map_w0(df) + k) =  gravity*(s*r - (1.0_r_def-s)*scaled_radius**2/r)

        end do
    end do

  else

    do k = 0, nlayers-1
        do dfc = 1, ndf_chi
            chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
        end do

        do df = 1, ndf_w0
            x(:) = 0.0_r_def
            do dfc = 1, ndf_chi
                x(3) = x(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
            end do

            phi(map_w0(df) + k) =  gravity*x(3)

        end do
    end do

  end if

end subroutine compute_geopotential_code

end module compute_geopotential_kernel_mod
