!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the boundary integral part of the projection operator from the velocity space to
!>        the potential temperature space weighted by the potential temperature gradient

!> @detail Kernel which computes the boundary integral part of the projection operator from the velocity space to
!>         the potential temperature space weighted by the potential temperature gradient
!>         Compute the projection operator \f<[<\gamma, flux(\theta*v)\cdot n>\f]
!>         where v is in W2 and gamma is in the potential temperature space
module weighted_proj_theta2_bd_kernel_mod
    use kernel_mod,              only : kernel_type
    use argument_mod,            only : arg_type, func_type,                       &
                                        GH_OPERATOR, GH_FIELD, GH_READ, GH_INC,    &
                                        W2, Wtheta, GH_BASIS,                      &
                                        CELLS, GH_QUADRATURE_XYoZ, QUADRATURE_XYoZ
    use constants_mod,           only : r_def, i_def
    use cross_product_mod,       only : cross_product
    use planet_config_mod,       only : cp
    use reference_element_mod,   only : nfaces_h, normal_to_face, out_face_normal



    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: weighted_proj_theta2_bd_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                               &
            arg_type(GH_OPERATOR, GH_INC, Wtheta, W2),                    &
            arg_type(GH_FIELD,   GH_READ, Wtheta)                         &
            /)
        type(func_type) :: meta_funcs(2) = (/                             &
            func_type(Wtheta, GH_BASIS),                                  &
            func_type(W2, GH_BASIS)                                       &
            /)
        integer :: iterates_over = CELLS
        integer :: gh_shape = GH_QUADRATURE_XYoZ
        ! gh_shape replaces evaluator_shape and will be removed by #1066
        integer :: evaluator_shape = QUADRATURE_XYoZ
    contains
        procedure, nopass ::weighted_proj_theta2_bd_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! overload the default structure constructor for function space
    interface weighted_proj_theta2_bd_kernel_type
        module procedure weighted_proj_theta2_bd_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public weighted_proj_theta2_bd_code
contains

    type(weighted_proj_theta2_bd_kernel_type) function weighted_proj_theta2_bd_kernel_constructor() result(self)
        return
    end function weighted_proj_theta2_bd_kernel_constructor

    !> @brief The subroutine which is called directly by the Psy layer
    !! @param[in] cell Cell number
    !! @param[in] nlayers Integer the number of layers
    !! @param[in] ncell_3d ncell*ndf
    !! @param[inout] projection Projection operator to compute
    !! @param[in] theta Potential temperature
    !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
    !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta Number of unique of degrees of freedom for wtheta
    !! @param[in] map_wtheta Dofmap for wtheta
    !! @param[inout] r_theta_bd Real array the data
    !! @param[in] nqp_h_1d Number of quadrature points in a single horizontal direction
    !! @param[in] nqp_v Integer, number of quadrature points in the vertical
    !! @param[in] wqp_v Real array. Quadrature weights vertical
    !! @param[in] w2_basis_face Real 5-dim array holding w2 basis functions evaluated at gaussian quadrature points on horizontal faces
    !! @param[in] wtheta_basis_face Real 5-dim array holding wtheta basis functions evaluated at gaussian quadrature points on horizontal faces
    !! @param[in] adjacent_face Vector containing information on neighbouring face index for the current cell

    subroutine weighted_proj_theta2_bd_code(cell, nlayers, ncell_3d,             &
                                            projection,                          &
                                            theta,                               &
                                            ndf_w2,                              &
                                            ndf_wtheta, undf_wtheta,             &
                                            stencil_wtheta_map,                  &
                                            stencil_wtheta_size,                 &
                                            nqp_h_1d, nqp_v, wqp_v,              &
                                            w2_basis_face,                       &
                                            wtheta_basis_face, adjacent_face )

        implicit none
        ! Arguments
        integer(kind=i_def), intent(in) :: cell, nlayers, ncell_3d, nqp_h_1d, nqp_v
        integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta, ndf_w2
        integer(kind=i_def), intent(in) :: stencil_wtheta_size

        integer(kind=i_def), dimension(ndf_wtheta, stencil_wtheta_size), intent(in) :: stencil_wtheta_map
        integer(kind=i_def), dimension(nfaces_h),                        intent(in) :: adjacent_face

        real(kind=r_def), dimension(ndf_wtheta,ndf_w2,ncell_3d),    intent(inout) :: projection
        real(kind=r_def), dimension(undf_wtheta),                   intent(in)    :: theta
        real(kind=r_def), dimension(3,ndf_w2,    nqp_h_1d,nqp_v,4), intent(in)    :: w2_basis_face
        real(kind=r_def), dimension(1,ndf_wtheta,nqp_h_1d,nqp_v,4), intent(in)    :: wtheta_basis_face
        real(kind=r_def), dimension(nqp_v),                         intent(in)    :: wqp_v

        !Internal variables
        integer(kind=i_def) :: df, k, ik, face, face_next, dft, df2
        integer(kind=i_def) :: qp1, qp2

        real(kind=r_def), dimension(ndf_wtheta) :: theta_e, theta_next_e

        real(kind=r_def) :: theta_at_fquad, theta_next_at_fquad, v_dot_n
        real(kind=r_def) :: flux_term, theta_av

        logical, parameter :: upwind = .false.


        ! Assumes same number of horizontal qp in x and y
        do k = 0, nlayers-1
          ik = k + 1 + (cell-1)*nlayers
          do face = 1, nfaces_h
            ! Storing opposite face number on neighbouring cell
            face_next = adjacent_face(face)

            ! Computing theta in adjacent cells
            do df = 1, ndf_wtheta
              theta_e(df)      = theta(stencil_wtheta_map(df, 1) + k )
              theta_next_e(df) = theta(stencil_wtheta_map(df, face+1) + k )
            end do

            do qp2 = 1, nqp_v
              do qp1 = 1, nqp_h_1d

                theta_at_fquad = 0.0_r_def
                theta_next_at_fquad = 0.0_r_def
                do df = 1, ndf_wtheta
                  theta_at_fquad       = theta_at_fquad      + theta_e(df)     *wtheta_basis_face(1,df,qp1,qp2,face)
                  theta_next_at_fquad  = theta_next_at_fquad + theta_next_e(df)*wtheta_basis_face(1,df,qp1,qp2,face_next)
                end do
                theta_av = 0.5_r_def * (theta_at_fquad + theta_next_at_fquad)

                do df2 = 1,ndf_w2
                  v_dot_n = dot_product(w2_basis_face(:,df2,qp1,qp2,face),out_face_normal(:, face))
                  flux_term = wqp_v(qp1)*wqp_v(qp2) * theta_av * v_dot_n
                  if (upwind) then                      
                    flux_term = flux_term + 0.5_r_def * abs(v_dot_n) * &
                              (theta_at_fquad - theta_next_at_fquad)       
                  end if
                  do dft = 1,ndf_wtheta
                    projection(dft,df2,ik) = projection(dft,df2,ik) &
                                           + wtheta_basis_face(1,dft,qp1,qp2,face) * flux_term

                  end do ! dft
                end do ! df2
              end do ! qp1
           end do ! qp2
         end do ! faces
      end do ! layers
    end subroutine weighted_proj_theta2_bd_code


end module weighted_proj_theta2_bd_kernel_mod
