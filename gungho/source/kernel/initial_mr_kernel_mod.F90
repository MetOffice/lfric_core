!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial mr field

!> @detail The kernel computes initial mixing ratio fields for mr in the same
!>         space as that of theta

module initial_mr_kernel_mod

    use argument_mod,                  only: arg_type,  &
        GH_FIELD, GH_WRITE, GH_READ, ANY_SPACE_9, GH_BASIS,  &
        CELLS
    use fs_continuity_mod, only : W3, Wtheta
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use planet_config_mod, only : p_zero, Rd, kappa

    !physics routines
    use physics_common_mod, only: qsaturation
    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_mr_kernel_type
        private
        type(arg_type) :: meta_args(4) = (/              &
            arg_type(GH_FIELD, GH_READ, WTHETA),         &
            arg_type(GH_FIELD, GH_READ, W3),             &
            arg_type(GH_FIELD*5, GH_WRITE, WTHETA),      &
            arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)   &
            /)
        integer :: iterates_over = CELLS

    contains
        procedure, nopass :: initial_mr_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! overload the default structure constructor for function space
    interface initial_mr_kernel_type
        module procedure initial_mr_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public initial_mr_code
contains

    type(initial_mr_kernel_type) function initial_mr_kernel_constructor() result(self)
        return
    end function initial_mr_kernel_constructor

    !> @brief The subroutine which is called directly by the Psy layer
    !! @param[in] nlayers Integer the number of layers
    !! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta The number of total degrees of freedom for wtheta
    !! @param[in] map_wtheta Integer array holding the dofmap for the cell at the base of the column
    !! @param[in] theta Real array the data: theta
    !! @param[in] rho Real array the data: dry rho
    !! @param[inout] mr_v Real array the data: vapour
    !! @param[inout] mr_c Real array the data: cloud
    !! @param[inout] mr_r Real array the data: rain
    !! @param[inout] mr_nc Real array the data: cloud number
    !! @param[inout] mr_nr Real array the data: rain number
    !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
    !! @param[in] undf_chi Number of total degrees of freedom for chi
    !! @param[in] map_chi Dofmap for the cell at the base of the column
    !! @param[in] chi_1 X component of the chi coordinate field
    !! @param[in] chi_2 Y component of the chi coordinate field
    !! @param[in] chi_3 Z component of the chi coordinate field
    subroutine initial_mr_code(nlayers, ndf_wtheta, undf_wtheta, map_wtheta,     &
                               theta, rho, ndf_w3, undf_w3, map_w3,       &
                               mr_v, mr_c, mr_r, mr_nc, mr_nr,                   &
                               ndf_chi, undf_chi, map_chi, chi_1, chi_2, chi_3)

        implicit none

        !Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, ndf_chi, undf_wtheta, undf_chi
        integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
        integer(kind=i_def), dimension(ndf_wtheta), intent(in)  :: map_wtheta
        integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3
        integer(kind=i_def), dimension(ndf_chi), intent(in)     :: map_chi
        real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_v, mr_c, mr_r, mr_nc, mr_nr
        real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
        real(kind=r_def), dimension(undf_w3), intent(in)        :: rho
        real(kind=r_def), dimension(undf_chi), intent(in)       :: chi_1, chi_2, chi_3

        !Internal variables
        integer(kind=i_def)                 :: k, df, kp1

        real(kind=r_def)                    :: theta_at_dof, rho_at_dof, pressure_at_dof, &
                                               exner_at_dof, temperature_at_dof
        ! compute the pointwise mr profile
        do k = 0, nlayers-1
          kp1=min(k+1,nlayers-1)
          do df = 1, ndf_wtheta
            theta_at_dof=theta(map_wtheta(df) + k)
            rho_at_dof=0.5*(rho(map_w3(1) + k) + rho(map_w3(1) + kp1))
            pressure_at_dof = p_zero * &
               (rho_at_dof*Rd/p_zero*theta_at_dof)**(1.0_r_def/(1.0_r_def-kappa))
            exner_at_dof    = (pressure_at_dof / p_zero ) ** kappa
            temperature_at_dof=theta_at_dof*exner_at_dof
            mr_v(map_wtheta(df) + k) = 0.99_r_def *  &
               qsaturation(temperature_at_dof, 0.01_r_def*pressure_at_dof)
            mr_c(map_wtheta(df) + k) = 0.0_r_def
            mr_r(map_wtheta(df) + k) = 0.0_r_def
            mr_nc(map_wtheta(df) + k) = 0.0_r_def
            mr_nr(map_wtheta(df) + k) = 0.0_r_def
          end do
          
        end do

    end subroutine initial_mr_code

end module initial_mr_kernel_mod
