!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes boundary integral part of rhs of the momentum equation for the nonlinear equations,
!>         written in the vector invariant form

!> @details The kernel computes the boundary integral on rhs of the momentum equation for the nonlinear equations,
!>         written in the vector invariant form
!>         This consists of
!>         ru_bd = -cp*theta*v*normal_vector*average(pi)
!>
!>         where average(pi) needs to be considered as both rho and theta are discontinuous in the horizontal direction.
module ru_bd_kernel_mod
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  use kernel_mod,               only : kernel_type
  use argument_mod,             only : arg_type, func_type,                 &
    GH_FIELD, GH_READ, GH_INC,                                              &
    W2, W3, Wtheta, GH_BASIS,                                               &
    GH_DIFF_BASIS, CELLS,                                                   &
    GH_QUADRATURE_XYoZ, QUADRATURE_XYoZ
  use constants_mod,            only : r_def, i_def
  use cross_product_mod,        only : cross_product
  use planet_config_mod,        only : cp
  use reference_element_mod,    only : nfaces_h, out_face_normal


  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: ru_bd_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                               &
      arg_type(GH_FIELD,   GH_INC,  W2),                              &
      arg_type(GH_FIELD,   GH_READ, W3),                              &
      arg_type(GH_FIELD,   GH_READ, Wtheta)                           &
      /)
    type(func_type) :: meta_funcs(3) = (/                             &
      func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
      func_type(W3, GH_BASIS),                                        &
      func_type(Wtheta, GH_BASIS)                                     &
      /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
    ! gh_shape replaces evaluator_shape and will be removed by #1066
    integer :: evaluator_shape = QUADRATURE_XYoZ
  contains
    procedure, nopass ::ru_bd_code
  end type

  !-------------------------------------------------------------------------------
  ! Constructors
  !-------------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface ru_bd_kernel_type
    module procedure ru_bd_kernel_constructor
  end interface

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public ru_bd_code
contains

  type(ru_bd_kernel_type) function ru_bd_kernel_constructor() result(self)
    return
  end function ru_bd_kernel_constructor

  !> @brief Compute the boundary integral terms in the pressure gradient
  !! @param[in] nlayers Number of layers
  !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
  !! @param[in] undf_w2 Number unique of degrees of freedom  for w2
  !! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
  !! @param[inout] r_u_bd Right hand side of the momentum equation
  !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
  !! @param[in] undf_w3 Number unique of degrees of freedom  for w3
  !! @param[in] stencil_w3_map W3 dofmaps for the stencil
  !! @param[in] stencil_w3_size Size of the W3 stencil (number of cells)
  !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
  !! @param[in] undf_wtheta Number unique of degrees of freedom  for wtheta
  !! @param[in] stencil_wtheta_map W2 dofmaps for the stencil
  !! @param[in] stencil_wtheta_size Size of the W2 stencil (number of cells)
  !! @param[in] rho Density
  !! @param[in] theta Potential temperature
  !! @param[in] nqp_v Number of quadrature points in the vertical
  !! @param[in] nqp_h_1d Number of quadrature points in a single horizontal direction
  !! @param[in] wqp_v Vertical quadrature weights
  !! @param[in] w2_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] w3_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] wtheta_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] adjacent_face Vector containing information on neighbouring face index for the current cell

  subroutine ru_bd_code(nlayers,                                                    &
                        ndf_w2, undf_w2,                                            &
                        map_w2,                                                     &
                        ndf_w3, undf_w3,                                            &
                        stencil_w3_map,                                             &
                        stencil_w3_size,                                            &
                        ndf_wtheta, undf_wtheta,                                    &
                        stencil_wtheta_map,                                         &
                        stencil_wtheta_size,                                        &
                        r_u_bd,                                                     &
                        rho, theta,                                                 &
                        nqp_v, nqp_h_1d, wqp_v, w2_basis_face, w3_basis_face,       &
                        wtheta_basis_face, adjacent_face)


    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, nqp_v, nqp_h_1d
    integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2, undf_w3
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta
    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

    integer(kind=i_def), intent(in) :: stencil_w3_size
    integer(kind=i_def), dimension(ndf_w3, stencil_w3_size), intent(in)  :: stencil_w3_map

    integer(kind=i_def), intent(in) :: stencil_wtheta_size
    integer(kind=i_def), dimension(ndf_wtheta, stencil_wtheta_size), intent(in)  :: stencil_wtheta_map

    real(kind=r_def), dimension(4,3,ndf_w2,nqp_h_1d,nqp_v), intent(in)  :: w2_basis_face
    real(kind=r_def), dimension(4,1,ndf_w3,nqp_h_1d,nqp_v), intent(in)  :: w3_basis_face
    real(kind=r_def), dimension(4,1,ndf_wtheta,nqp_h_1d,nqp_v), intent(in) :: wtheta_basis_face

    integer(kind=i_def), dimension(nfaces_h), intent(in) :: adjacent_face

    real(kind=r_def), dimension(undf_w2), intent(inout)     :: r_u_bd
    real(kind=r_def), dimension(undf_w3), intent(in)        :: rho
    real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta

    real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

    !Internal variables
    integer(kind=i_def)              :: df, k, face, face_next
    integer(kind=i_def)              :: qp1, qp2

    real(kind=r_def), dimension(ndf_w3)          :: rho_e, rho_next_e
    real(kind=r_def), dimension(ndf_wtheta)      :: theta_e, theta_next_e
    real(kind=r_def), dimension(ndf_w2)          :: ru_bd_e

    real(kind=r_def) :: v(3)
    real(kind=r_def) :: exner_at_fquad, exner_next_at_fquad, rho_at_fquad
    real(kind=r_def) :: theta_at_fquad, theta_next_at_fquad, rho_next_at_fquad, bdary_term
    real(kind=r_def) :: av_pi_at_fquad

    do k = 0, nlayers-1

      do df = 1, ndf_w2
          ru_bd_e(df) = 0.0_r_def
      end do
      do face = 1, nfaces_h

        bdary_term     = 0.0_r_def
        av_pi_at_fquad = 0.0_r_def

        ! Storing opposite face number on neighbouring cell

        face_next = adjacent_face(face)

        ! Computing rho and theta in local and adjacent cell

        do df = 1, ndf_w3
          rho_e(df)      = rho( stencil_w3_map(df, 1) + k )
          rho_next_e(df) = rho( stencil_w3_map(df, face+1) + k )
        end do

        do df = 1, ndf_wtheta
          theta_e(df)      = theta( stencil_wtheta_map(df, 1) + k )
          theta_next_e(df) = theta( stencil_wtheta_map(df, face+1) + k )
        end do

        ! Compute the boundary RHS integrated over one horizontal face
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h_1d
            rho_at_fquad = 0.0_r_def
            rho_next_at_fquad = 0.0_r_def

            do df = 1, ndf_w3
              rho_at_fquad  = rho_at_fquad + rho_e(df)*w3_basis_face(face,1,df,qp1,qp2)
              rho_next_at_fquad  = rho_next_at_fquad + rho_next_e(df)*w3_basis_face(face_next,1,df,qp1,qp2)
            end do

            theta_at_fquad = 0.0_r_def
            theta_next_at_fquad = 0.0_r_def

            do df = 1, ndf_wtheta
              theta_at_fquad   = theta_at_fquad + theta_e(df)*wtheta_basis_face(face,1,df,qp1,qp2)
              theta_next_at_fquad  = theta_next_at_fquad + theta_next_e(df)*wtheta_basis_face(face_next,1,df,qp1,qp2)
            end do

            exner_at_fquad = calc_exner_pointwise(rho_at_fquad, theta_at_fquad)
            exner_next_at_fquad = calc_exner_pointwise(rho_next_at_fquad, theta_next_at_fquad)

            av_pi_at_fquad = .5 * (exner_at_fquad + exner_next_at_fquad)

            do df = 1, ndf_w2
              v  = w2_basis_face(face,:,df,qp1,qp2)

              bdary_term = - cp * dot_product(v, out_face_normal(:, face)) *  theta_at_fquad * av_pi_at_fquad
              ru_bd_e(df) = ru_bd_e(df) + wqp_v(qp1)*wqp_v(qp2) * bdary_term
            end do

          end do ! qp1
        end do ! qp2
      end do ! faces

      do df = 1, ndf_w2
        r_u_bd( map_w2(df) + k ) =  r_u_bd( map_w2(df) + k ) + ru_bd_e(df)
      end do

    end do ! layers

  end subroutine ru_bd_code

end module ru_bd_kernel_mod
