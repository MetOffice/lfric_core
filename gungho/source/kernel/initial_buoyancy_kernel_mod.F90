!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial theta field

!> @detail The kernel computes initial theta perturbation field for theta in the space
!>         of horizontally discontinuous, vertically continuous polynomials

module initial_buoyancy_kernel_mod

    use argument_mod, only: arg_type, func_type,        &
        GH_FIELD, GH_WRITE, GH_READ,                    &
        ANY_SPACE_9, ANY_SPACE_1, GH_BASIS,             &
        GH_DIFF_BASIS,                                  &
        CELLS, GH_EVALUATOR, EVALUATOR
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type

    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_buoyancy_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                               &
            arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1),                  &
            arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                    &
            /)
        type(func_type) :: meta_funcs(1) = (/                             &
            func_type(ANY_SPACE_9, GH_BASIS)                              &
            /)
        integer :: iterates_over = CELLS
        integer :: gh_shape = GH_EVALUATOR
        ! gh_shape replaces evaluator_shape and will be removed by #1066
        integer :: evaluator_shape = EVALUATOR
    contains
        procedure, nopass :: initial_buoyancy_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! Overload the default structure constructor for function space
    interface initial_buoyancy_kernel_type
        module procedure initial_buoyancy_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public initial_buoyancy_code
contains

    type(initial_buoyancy_kernel_type) function initial_buoyancy_kernel_constructor() result(self)
        return
    end function initial_buoyancy_kernel_constructor

      !> @brief The subroutine which is called directly by the Psy layer
      !! @param[in] nlayers Integer the number of layers
      !! @param[in] ndf_wt The number of degrees of freedom per cell for wt
      !! @param[in] udf_wt The number of total degrees of freedom for wt
      !! @param[in] map_wt Integer array holding the dofmap for the cell at the base of the column
      !! @param[inout] theta Real array the data
      !! @param[in] wt_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
      !! @param[in] ndf_w0 The number of degrees of freedom per cell
      !! @param[in] ndf_w0 The total number of degrees of freedom
      !! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column
      !! @param[in] w0_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
      !! @param[inout] chi_1 Real array, the x component of the w0 coordinate field
      !! @param[inout] chi_2 Real array, the y component of the w0 coordinate field
      !! @param[inout] chi_3 Real array, the z component of the w0 coordinate field

      ! In Psyclone

    subroutine initial_buoyancy_code(nlayers, &
                                     buoyancy, &
                                     chi_1, chi_2, chi_3, &
                                     ndf_wt, undf_wt, map_wt, &
                                     ndf_w0, undf_w0, map_w0, w0_basis)

        use analytic_buoyancy_profiles_mod, only : analytic_buoyancy

        implicit none

        !Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wt, ndf_w0, undf_wt, undf_w0
        integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
        integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0
        real(kind=r_def), dimension(undf_wt),          intent(inout) :: buoyancy
        real(kind=r_def), dimension(undf_w0),          intent(in)    :: chi_1, chi_2, chi_3
        real(kind=r_def), dimension(1,ndf_w0,ndf_wt),  intent(in)    :: w0_basis

        !Internal variables
        integer(kind=i_def)                 :: df, df0, k
        real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
        real(kind=r_def)                    :: x(3)

        ! compute the pointwise buoyancy profile
        do k = 0, nlayers-1
          do df0 = 1, ndf_w0
            chi_1_e(df0) = chi_1( map_w0(df0) + k)
            chi_2_e(df0) = chi_2( map_w0(df0) + k)
            chi_3_e(df0) = chi_3( map_w0(df0) + k)
          end do

          do df = 1, ndf_wt
            x(:) = 0.0_r_def
            do df0 = 1, ndf_w0
              x(1) = x(1) + chi_1_e(df0)*w0_basis(1,df0,df)
              x(2) = x(2) + chi_2_e(df0)*w0_basis(1,df0,df)
              x(3) = x(3) + chi_3_e(df0)*w0_basis(1,df0,df)
            end do

            buoyancy(map_wt(df) + k) = analytic_buoyancy(x)
          end do
        end do

    end subroutine initial_buoyancy_code

end module initial_buoyancy_kernel_mod
