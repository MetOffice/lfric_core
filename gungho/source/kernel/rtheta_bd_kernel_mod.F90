!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the boundary integral part of rhs of the thermodynamic equation for the nonlinear equations


!> @detail The kernel computes the boundary integral on rhs of the thermodynamic equation for the nonlinear equations
!>         This consists of
!>         rtheta_bd = theta * gamma * u * normal
module rtheta_bd_kernel_mod
    use kernel_mod,              only : kernel_type
    use argument_mod,            only : arg_type, func_type,                 &
                                        GH_FIELD, GH_READ, GH_INC,           &
                                        W2, Wtheta, GH_BASIS, GH_DIFF_BASIS, &
                                        CELLS, GH_QUADRATURE_XYoZ,           &
                                        QUADRATURE_XYoZ
    use constants_mod,           only : r_def, i_def
    use cross_product_mod,       only : cross_product
    use planet_config_mod,       only : cp
    use reference_element_mod,   only : nfaces_h, normal_to_face, out_face_normal



    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: rtheta_bd_kernel_type
        private
        type(arg_type) :: meta_args(3) = (/                               &
            arg_type(GH_FIELD,   GH_INC,  Wtheta),                        &
            arg_type(GH_FIELD,   GH_READ, Wtheta),                        &
            arg_type(GH_FIELD,   GH_READ, W2)                             &
            /)
        type(func_type) :: meta_funcs(2) = (/                             &
            func_type(W2, GH_BASIS),                                      &
            func_type(Wtheta, GH_BASIS)                                   &
            /)
        integer :: iterates_over = CELLS
        integer :: gh_shape = GH_QUADRATURE_XYoZ
        ! gh_shape replaces evaluator_shape and will be removed by #1066
        integer :: evaluator_shape = QUADRATURE_XYoZ
    contains
        procedure, nopass ::rtheta_bd_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! overload the default structure constructor for function space
    interface rtheta_bd_kernel_type
        module procedure rtheta_bd_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public rtheta_bd_code
contains

    type(rtheta_bd_kernel_type) function rtheta_bd_kernel_constructor() result(self)
        return
    end function rtheta_bd_kernel_constructor

    !> @brief The subroutine which is called directly by the Psy layer
    !! @param[in] nlayers Integer the number of layers
    !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
    !! @param[in] undf_w2 Number of unique of degrees of freedom  for w2
    !! @param[in] stencil_w2_map W2 dofmaps for the stencil
    !! @param[in] stencil_w2_size Size of the W2 stencil (number of cells)
    !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta Number of unique of degrees of freedom for wtheta
    !! @param[in] stencil_wtheta_map W2 dofmaps for the stencil
    !! @param[in] stencil_wtheta_size Size of the W2 stencil (number of cells)
    !! @param[inout] r_theta_bd Real array the data
    !! @param[in] theta Real array. potential temperature
    !! @param[in] u wind field
    !! @param[in] nqp_v Integer, number of quadrature points in the vertical
    !! @param[in] nqp_h_1d Integer, number of quadrature points in a single horizontal direction
    !! @param[in] wqp_h Real array. Quadrature weights horizontal
    !! @param[in] wqp_v Real array. Quadrature weights vertical
    !! @param[in] w2_basis_face Real 5-dim array holding w2 basis functions evaluated at gaussian quadrature points on horizontal faces
    !! @param[in] wtheta_basis_face Real 5-dim array holding wtheta basis functions evaluated at gaussian quadrature points on horizontal faces
    !! @param[in] adjacent_face Vector containing information on neighbouring face index for the current cell

    subroutine rtheta_bd_code(nlayers,                                                    &
                              ndf_w2, undf_w2,                                            &
                              stencil_w2_map,                                             &
                              stencil_w2_size,                                            &
                              ndf_wtheta, undf_wtheta,                                    &
                              stencil_wtheta_map,                                         &
                              stencil_wtheta_size,                                        &
                              r_theta_bd,                                                 &
                              theta, u,                                                   &
                              nqp_v, nqp_h_1d, wqp_v, w2_basis_face,                      &
                              wtheta_basis_face, adjacent_face )


        ! Arguments
        integer(kind=i_def), intent(in) :: nlayers, nqp_h_1d, nqp_v
        integer(kind=i_def), intent(in) :: ndf_w2
        integer(kind=i_def), intent(in) :: undf_w2
        integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta

        integer(kind=i_def), intent(in) :: stencil_w2_size
        integer(kind=i_def), intent(in) :: stencil_wtheta_size

        integer(kind=i_def), dimension(ndf_w2, stencil_w2_size), intent(in)          :: stencil_w2_map
        integer(kind=i_def), dimension(ndf_wtheta, stencil_wtheta_size), intent(in)  :: stencil_wtheta_map

        real(kind=r_def), dimension(4,3,ndf_w2,nqp_h_1d,nqp_v), intent(in)     :: w2_basis_face
        real(kind=r_def), dimension(4,1,ndf_wtheta,nqp_h_1d,nqp_v), intent(in) :: wtheta_basis_face

        integer(kind=i_def), dimension(nfaces_h), intent(in) :: adjacent_face

        real(kind=r_def), dimension(undf_wtheta), intent(inout) :: r_theta_bd
        real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
        real(kind=r_def), dimension(undf_w2), intent(in)        :: u

        real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

        !Internal variables
        integer(kind=i_def)               :: df, k, face, face_next
        integer(kind=i_def)               :: qp1, qp2

        real(kind=r_def), dimension(ndf_wtheta)      :: theta_e, theta_next_e
        real(kind=r_def), dimension(ndf_wtheta)      :: rtheta_bd_e
        real(kind=r_def), dimension(ndf_w2)          :: u_e, u_next_e

        real(kind=r_def) :: u_at_uquad(3), u_next_at_uquad(3), face_next_inward_normal(3)
        real(kind=r_def) :: theta_at_uquad, theta_next_at_uquad
        real(kind=r_def) :: bdary_term, gamma_wtheta, sign_face_next_outward, flux_term

        logical          :: upwind = .false.

        ! Assumes same number of horizontal qp in x and y

        do k = 0, nlayers-1

            do df = 1, ndf_wtheta
                rtheta_bd_e(df) = 0.0_r_def
            end do
            do face = 1, nfaces_h

              ! Storing opposite face number on neighbouring cell

              face_next = adjacent_face(face)

              sign_face_next_outward = (-1.0_r_def)**(int(floor(real(mod(face_next, 4))/2.0) + 1.0_r_def))
              face_next_inward_normal(:) = -sign_face_next_outward * normal_to_face(face_next, :)



              ! Computing theta and f in adjacent cells

              do df = 1, ndf_w2
                u_e(df)      = u( stencil_w2_map(df, 1) + k )
                u_next_e(df) = u( stencil_w2_map(df, face+1) + k )
              end do

              do df = 1, ndf_wtheta
                theta_e(df)      = theta( stencil_wtheta_map(df, 1) + k )
                theta_next_e(df) = theta( stencil_wtheta_map(df, face+1) + k )
              end do

              ! compute the boundary RHS integrated over one cell
              do qp2 = 1, nqp_v
                do qp1 = 1, nqp_h_1d
                  theta_at_uquad = 0.0_r_def
                  theta_next_at_uquad = 0.0_r_def

                  do df = 1, ndf_wtheta
                    theta_at_uquad       = theta_at_uquad + theta_e(df)*wtheta_basis_face(face,1,df,qp1,qp2)
                    theta_next_at_uquad  = theta_next_at_uquad + theta_next_e(df)*wtheta_basis_face(face_next,1,df,qp1,qp2)
                  end do

                  u_at_uquad(:) = 0.0_r_def
                  u_next_at_uquad(:) = 0.0_r_def

                  do df = 1, ndf_w2
                    u_at_uquad(:)       = u_at_uquad(:)      + u_e(df)     *w2_basis_face(face,     :,df,qp1,qp2)
                    u_next_at_uquad(:)  = u_next_at_uquad(:) + u_next_e(df)*w2_basis_face(face_next,:,df,qp1,qp2)
                  end do

                  flux_term = 0.5_r_def * (theta_next_at_uquad * dot_product(u_next_at_uquad, face_next_inward_normal) + &
                                           theta_at_uquad      * dot_product(u_at_uquad,      out_face_normal(:, face)))

                  if (upwind) then
                    flux_term = flux_term + 0.5_r_def * abs(dot_product(u_at_uquad, out_face_normal(:, face))) * &
                                            (dot_product(theta_at_uquad      * out_face_normal(:, face), out_face_normal(:, face)) - &
                                             dot_product(theta_next_at_uquad * face_next_inward_normal , face_next_inward_normal))
                  end if

                  do df = 1, ndf_wtheta
                    gamma_wtheta  = wtheta_basis_face(face,1,df,qp1,qp2)

                    bdary_term = gamma_wtheta * flux_term
                    rtheta_bd_e(df) = rtheta_bd_e(df) +  wqp_v(qp1)*wqp_v(qp2) * bdary_term
                  end do

                end do ! qp1
              end do ! qp2

          end do ! faces

          do df = 1, ndf_wtheta
            r_theta_bd( stencil_wtheta_map(df, 1) + k ) =  r_theta_bd( stencil_wtheta_map(df, 1) + k ) + rtheta_bd_e(df)
          end do

        end do ! layers

    end subroutine rtheta_bd_code


end module rtheta_bd_kernel_mod
