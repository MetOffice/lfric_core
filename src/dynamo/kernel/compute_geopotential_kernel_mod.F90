!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the geopotential field

!> @details Computes the geopotential field Phi = g*r or g*z for Cartesian
!!         domains

module compute_geopotential_kernel_mod

use argument_mod,         only : arg_type, func_type,                      &
                                 GH_FIELD, GH_READ, GH_WRITE,              &
                                 W0, ANY_SPACE_9, GH_BASIS,                &
                                 CELLS
use base_mesh_config_mod, only : geometry, &
                                 base_mesh_geometry_spherical
use constants_mod,        only : r_def
use coord_transform_mod,  only : xyz2llr
use kernel_mod,           only : kernel_type
use planet_config_mod,    only : gravity

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: compute_geopotential_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W0),                             &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                               &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: compute_geopotential_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_geopotential_kernel_type
   module procedure compute_geopotential_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_geopotential_code
contains

type(compute_geopotential_kernel_type) function compute_geopotential_kernel_constructor() result(self)
  return
end function compute_geopotential_kernel_constructor

!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] undf_w0 The number of unique degrees of freedom for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[inout] phi The actual data
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi The dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis 3-dim array chi basis functions evaluated at w0 nodes
!! @param[in] chi_1 Real array, the physical x coordinates
!! @param[in] chi_2 Real array, the physical y coordinates
!! @param[in] chi_3 Real array, the physical z coordinates
subroutine compute_geopotential_code(nlayers,phi, &
                                     chi_1,chi_2,chi_3, &
                                     ndf_w0,undf_w0,map_w0, &
                                     ndf_chi, undf_chi, map_chi, &
                                     chi_basis)

  !Arguments
  integer, intent(in) :: nlayers, ndf_w0, undf_w0, ndf_chi, undf_chi
  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(undf_w0), intent(inout)   :: phi
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_2
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_3
  real(kind=r_def), dimension(1,ndf_chi,ndf_w0),  intent(in)    :: chi_basis

  !Internal variables
  integer          :: df, dfc, k
  real(kind=r_def) :: x(3)
  real(kind=r_def) :: lat, lon, r

  real(kind=r_def), dimension(ndf_chi)             :: chi_1_e, chi_2_e, chi_3_e

  if ( geometry == base_mesh_geometry_spherical ) then
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

            phi(map_w0(df) + k) =  gravity*r

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
