!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes the initial theta field

!> @details The kernel computes initial theta perturbation field for theta in the space
!>          of horizontally discontinuous, vertically continuous polynomials

module initial_theta_kernel_mod

    use argument_mod, only: arg_type, func_type,        &
        GH_FIELD, GH_WRITE, GH_READ,                    &
        ANY_SPACE_9, ANY_SPACE_1, GH_BASIS,             &
        GH_DIFF_BASIS,                                  &
        CELLS, GH_EVALUATOR, EVALUATOR
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use idealised_config_mod,          only: test

    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_theta_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                               &
            arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1),                  &
            arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                    &
            /)
        type(func_type) :: meta_funcs(1) = (/                             &
             func_type(ANY_SPACE_9, GH_BASIS)                             &
             /)
        integer :: iterates_over = CELLS
        integer :: gh_shape = GH_EVALUATOR
        ! gh_shape replaces evaluator_shape and will be removed by #1066
        integer :: evaluator_shape = EVALUATOR
    contains
        procedure, nopass :: initial_theta_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! Overload the default structure constructor for function space
    interface initial_theta_kernel_type
        module procedure initial_theta_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public initial_theta_code
contains

    type(initial_theta_kernel_type) function initial_theta_kernel_constructor() result(self)
        return
    end function initial_theta_kernel_constructor

    !> @brief Computes the initial theta field
    !! @param[in] nlayers Number of layers
    !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta Number of total degrees of freedom for wtheta
    !! @param[in] map_wtheta Dofmap for the cell at the base of the column
    !! @param[inout] theta Potential temperature
    !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
    !! @param[in] undf_chi Number of total degrees of freedom for chi
    !! @param[in] map_chi Dofmap for the cell at the base of the column
    !! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
    !! @param[inout] chi_1 X component of the chi coordinate field
    !! @param[inout] chi_2 Y component of the chi coordinate field
    !! @param[inout] chi_3 Z component of the chi coordinate field
    subroutine initial_theta_code(nlayers, ndf_wtheta, undf_wtheta, map_wtheta, theta, &
                                  ndf_chi, undf_chi, map_chi, chi_basis, chi_1, chi_2, chi_3)

        use analytic_temperature_profiles_mod, only : analytic_temperature

        implicit none

        !Arguments
        integer(kind=i_def),                        intent(in) :: nlayers, ndf_wtheta, ndf_chi, undf_wtheta, undf_chi
        integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
        integer(kind=i_def), dimension(ndf_chi),    intent(in) :: map_chi

        real(kind=r_def),    dimension(undf_wtheta),           intent(out) :: theta
        real(kind=r_def),    dimension(undf_chi),              intent(in)  :: chi_1, chi_2, chi_3
        real(kind=r_def),    dimension(1,ndf_chi,ndf_wtheta),  intent(in)  :: chi_basis

        !Internal variables
        integer(kind=i_def)                    :: df, dfc, k
        real(kind=r_def),   dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
        real(kind=r_def)                       :: x(3), x_surf(3)

        ! compute the pointwise theta profile

        do k = 0, nlayers-1
          do dfc = 1, ndf_chi
            chi_1_e(dfc) = chi_1( map_chi(dfc) + k)
            chi_2_e(dfc) = chi_2( map_chi(dfc) + k)
            chi_3_e(dfc) = chi_3( map_chi(dfc) + k)
          end do

          do df = 1, ndf_wtheta
            x(:) = 0.0_r_def
            do dfc = 1, ndf_chi
              x(1) = x(1) + chi_1_e(dfc)*chi_basis(1,dfc,df)
              x(2) = x(2) + chi_2_e(dfc)*chi_basis(1,dfc,df)
              x(3) = x(3) + chi_3_e(dfc)*chi_basis(1,dfc,df)
            end do
            if(df == 1 .and. k .eq. 0 ) x_surf = x

            theta(map_wtheta(df) + k) = analytic_temperature(x, test, x_surf)

          end do
        end do

    end subroutine initial_theta_code

end module initial_theta_kernel_mod
