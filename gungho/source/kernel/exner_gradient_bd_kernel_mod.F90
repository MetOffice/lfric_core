!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes boundary integral of the pressure gradient term for the nonlinear equations,
!>         written in the vector invariant form

!> @details The kernel computes the boundary integral of the pressure gradient term for the nonlinear equations,
!>         written in the vector invariant form
!>         This consists of
!>         ru_bd = -cp*theta*v*normal_vector*average(pi)
!>
!>         where average(pi) needs to be considered as  exner is discontinuous in the horizontal direction.
module exner_gradient_bd_kernel_mod
  use kernel_mod,              only : kernel_type
  use argument_mod,            only : arg_type, func_type, mesh_data_type, &
                                      GH_FIELD, GH_READ, GH_INC,           &
                                      W2, W3, Wtheta, GH_BASIS,            &
                                      GH_DIFF_BASIS, CELLS,                &
                                      GH_QUADRATURE_XYoZ,                  &
                                      adjacent_face,                       &
                                      reference_element_out_face_normal
  use constants_mod,           only : r_def, i_def
  use cross_product_mod,       only : cross_product
  use planet_config_mod,       only : cp

  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: exner_gradient_bd_kernel_type
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
        type(mesh_data_type) :: meta_init(2) = (/                         &
            mesh_data_type( adjacent_face ),                              &
            mesh_data_type( reference_element_out_face_normal )           &
          /)
  contains
    procedure, nopass ::exner_gradient_bd_code
  end type

  !-------------------------------------------------------------------------------
  ! Constructors
  !-------------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface exner_gradient_bd_kernel_type
    module procedure exner_gradient_bd_kernel_constructor
  end interface

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public exner_gradient_bd_code
contains

  type(exner_gradient_bd_kernel_type) function exner_gradient_bd_kernel_constructor() result(self)
    return
  end function exner_gradient_bd_kernel_constructor

  !> @brief Compute the boundary integral terms in the pressure gradient
  !>
  !! @param[in] nlayers Number of layers
  !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
  !! @param[in] undf_w2 Number unique of degrees of freedom  for w2
  !! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
  !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
  !! @param[in] undf_w3 Number unique of degrees of freedom  for w3
  !! @param[in] stencil_w3_map W3 dofmaps for the stencil
  !! @param[in] stencil_w3_size Size of the W3 stencil (number of cells)
  !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
  !! @param[in] undf_wtheta Number unique of degrees of freedom  for wtheta
  !! @param[in] map_wtheta Integer array holding the dofmap for the cell at the base of the column for wtheta
  !! @param[inout] r_u_bd Right hand side of the momentum equation
  !! @param[in] exner The Exner pressure
  !! @param[in] theta Potential temperature
  !! @param[in] nqp_v Number of quadrature points in the vertical
  !! @param[in] nqp_h_1d Number of quadrature points in a single horizontal direction
  !! @param[in] wqp_v Vertical quadrature weights
  !! @param[in] w2_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] w3_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] wtheta_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] adjacent_face Vector containing information on neighbouring face index for the current cell
  !! @param[in] out_face_normal Vector of normals to the reference element
  !!                            "out face".
  !!
  subroutine exner_gradient_bd_code( nlayers,                      &
                                     ndf_w2, undf_w2,              &
                                     map_w2,                       &
                                     ndf_w3, undf_w3,              &
                                     stencil_w3_map,               &
                                     stencil_w3_size,              &
                                     ndf_wtheta, undf_wtheta,      &
                                     map_wtheta,                   &
                                     r_u_bd,                       &
                                     exner, theta,                 &
                                     nqp_v, nqp_h_1d, wqp_v,       &
                                     w2_basis_face, w3_basis_face, &
                                     wtheta_basis_face,            &
                                     adjacent_face, out_face_normal )

    implicit none

    ! Arguments
    integer(i_def), intent(in)    :: nlayers
    integer(i_def), intent(in)    :: ndf_w2
    integer(i_def), intent(in)    :: undf_w2
    integer(i_def), intent(in)    :: map_w2(ndf_w2)
    integer(i_def), intent(in)    :: ndf_w3
    integer(i_def), intent(in)    :: undf_w3

    ! These two arguments must be defined out of declaration order so as to
    ! avoid type conflicts.
    integer(i_def), intent(in)    :: stencil_w3_size
    integer(i_def), intent(in)    :: stencil_w3_map(ndf_w3, stencil_w3_size)

    integer(i_def), intent(in)    :: ndf_wtheta
    integer(i_def), intent(in)    :: undf_wtheta
    integer(i_def), intent(in)    :: map_wtheta(ndf_wtheta)
    real(r_def),    intent(inout) :: r_u_bd(undf_w2)
    real(r_def),    intent(in)    :: exner(undf_w3)
    real(r_def),    intent(in)    :: theta(undf_wtheta)
    integer(i_def), intent(in)    :: nqp_v
    integer(i_def), intent(in)    :: nqp_h_1d
    real(r_def),    intent(in)    :: wqp_v(nqp_v)
    real(r_def),    intent(in)    :: w2_basis_face(4, 3, &
                                                   ndf_w2, nqp_h_1d, nqp_v)
    real(r_def),    intent(in)    :: w3_basis_face(4, 1, &
                                                   ndf_w3, nqp_h_1d, nqp_v)
    real(r_def),    intent(in)    :: wtheta_basis_face(4, 1, &
                                                       ndf_wtheta, &
                                                       nqp_h_1d, nqp_v)
    integer(i_def), intent(in)    :: adjacent_face(:)
    real(r_def),    intent(in)    :: out_face_normal(:,:)


    !Internal variables
    integer(kind=i_def)              :: df, k, face, face_next
    integer(kind=i_def)              :: qp1, qp2

    real(kind=r_def), dimension(ndf_w3)          :: exner_e, exner_next_e
    real(kind=r_def), dimension(ndf_wtheta)      :: theta_e
    real(kind=r_def), dimension(ndf_w2)          :: ru_bd_e

    real(kind=r_def) :: v(3)
    real(kind=r_def) :: exner_at_fquad, exner_next_at_fquad
    real(kind=r_def) :: theta_at_fquad, bdary_term, av_pi_at_fquad

    do k = 0, nlayers-1

      do df = 1, ndf_w2
        ru_bd_e(df) = 0.0_r_def
      end do
      do face = 1, size( adjacent_face, 1 )

        bdary_term     = 0.0_r_def
        av_pi_at_fquad = 0.0_r_def

        ! Storing opposite face number on neighbouring cell

        face_next = adjacent_face(face)

        ! Computing exner and theta in local and adjacent cell

        do df = 1, ndf_w3
          exner_e(df)      = exner( stencil_w3_map(df, 1) + k )
          exner_next_e(df) = exner( stencil_w3_map(df, face+1) + k )
        end do

        do df = 1, ndf_wtheta
          theta_e(df)      = theta( map_wtheta(df) + k )
        end do

        ! Compute the boundary RHS integrated over one horizontal face
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h_1d
            exner_at_fquad = 0.0_r_def
            exner_next_at_fquad = 0.0_r_def

            do df = 1, ndf_w3
              exner_at_fquad  = exner_at_fquad + exner_e(df)*w3_basis_face(face,1,df,qp1,qp2)
              exner_next_at_fquad  = exner_next_at_fquad + exner_next_e(df)*w3_basis_face(face_next,1,df,qp1,qp2)
            end do

            theta_at_fquad = 0.0_r_def

            do df = 1, ndf_wtheta
              theta_at_fquad   = theta_at_fquad + theta_e(df)*wtheta_basis_face(face,1,df,qp1,qp2)
            end do

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

  end subroutine exner_gradient_bd_code

end module exner_gradient_bd_kernel_mod
