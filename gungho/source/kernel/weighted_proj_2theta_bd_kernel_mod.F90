!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the boundary part of the projection operator from the
!>        potential temperature space to the velocity space weighted by the
!>        pressure gradient.
!>
!> Compute the boundary projection operator \f[<v.n,{\Pi}*\gamma>\f] where v
!> is in W2, gamma is in the potential temperature space and exner is computed
!> pointwise from the equation of state.
!>
module weighted_proj_2theta_bd_kernel_mod

  use argument_mod,      only: arg_type, func_type, mesh_data_type, &
                               GH_OPERATOR, GH_FIELD, GH_REAL,      &
                               GH_READ, GH_READWRITE,               &
                               GH_BASIS, GH_DIFF_BASIS,             &
                               CELLS, GH_QUADRATURE_XYoZ,           &
                               adjacent_face,                       &
                               reference_element_out_face_normal
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W2, W3, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: weighted_proj_2theta_bd_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                  &
        arg_type(GH_OPERATOR, GH_READWRITE, W2, Wtheta), &
        arg_type(GH_FIELD,    GH_READ,      W3),         &
        arg_type(GH_REAL,     GH_READ)                   &
        /)
    type(func_type) :: meta_funcs(3) = (/ &
        func_type(W2,     GH_BASIS),      &
        func_type(Wtheta, GH_BASIS),      &
        func_type(W3,     GH_BASIS)       &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
      type(mesh_data_type) :: meta_init(2) = (/               &
          mesh_data_type( adjacent_face ),                    &
          mesh_data_type( reference_element_out_face_normal ) &
        /)
  contains
    procedure, nopass :: weighted_proj_2theta_bd_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface weighted_proj_2theta_bd_kernel_type
    module procedure weighted_proj_2theta_bd_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public weighted_proj_2theta_bd_code

contains

  type(weighted_proj_2theta_bd_kernel_type) &
  function weighted_proj_2theta_bd_constructor() result(self)
    implicit none
    return
  end function weighted_proj_2theta_bd_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Computes the weigthed projection from Wtheta to W2.
  !>
  !! @param[in] cell Cell number
  !! @param[in] nlayers Number of layers.
  !! @param[in] ncell_3d ncell*ndf
  !! @param[inout] projection Projection operator to compute
  !! @param[in] exner Exner presssure
  !! @param[in] scalar Real to scale matrix by
  !! @param[in] ndf_w2 Number of degrees of freedom per cell.
  !! @param[in] ndf_wtheta Number of degrees of freedom per cell.
  !! @param[in] ndf_w3 Number of degrees of freedom per cell.
  !! @param[in] undf_w3 Total number of degrees.
  !! @param[in] stencil_w3_map W3 dofmaps for the stencil
  !! @param[in] stencil_w3_size Size of the W3 stencil (number of cells)
  !! @param[in] nqp_h_1d  Number of quadrature points in a single horizontal &
  !!                      direction.
  !! @param[in] nqp_v Number of vertical quadrature points
  !! @param[in] wqp_v Vertical quadrature weights
  !! @param[in] w2_basis_face  Basis functions evaluated at gaussian &
  !!                           quadrature points on horizontal faces.
  !! @param[in] w3_basis_face  Basis functions evaluated at gaussian &
  !!                           quadrature points on horizontal faces.
  !! @param[in] wtheta_basis_face  Basis functions evaluated at gaussian &
  !!                               quadrature points on horizontal faces.
  !! @param[in] adjacent_face  Vector containing information on neighbouring &
  !!                           face index for the current cell.
  !! @param[in] out_face_normal  Vectors on "out faces" of reference element.
  !!
  subroutine weighted_proj_2theta_bd_code( cell, nlayers, ncell_3d, &
                                           projection,              &
                                           exner,                   &
                                           scalar,                  &
                                           ndf_w2,                  &
                                           ndf_wtheta,              &
                                           ndf_w3, undf_w3,         &
                                           stencil_w3_map,          &
                                           stencil_w3_size,         &
                                           nqp_h_1d, nqp_v, wqp_v,  &
                                           w2_basis_face,           &
                                           w3_basis_face,           &
                                           wtheta_basis_face,       &
                                           adjacent_face, out_face_normal )

    use calc_exner_pointwise_mod, only: calc_exner_pointwise

    implicit none

    ! Arguments
    integer(kind=i_def),                     intent(in) :: cell, nqp_h_1d, nqp_v
    integer(kind=i_def),                     intent(in) :: nlayers
    integer(kind=i_def),                     intent(in) :: ncell_3d
    integer(kind=i_def),                     intent(in) :: undf_w3, ndf_w3, ndf_w2, ndf_wtheta

    integer(kind=i_def), intent(in) :: stencil_w3_size
    integer(kind=i_def), dimension(ndf_w3, stencil_w3_size), intent(in)  :: stencil_w3_map

    real(kind=r_def), dimension(3,ndf_w2,nqp_h_1d,nqp_v,4), intent(in)     :: w2_basis_face
    real(kind=r_def), dimension(1,ndf_w3,nqp_h_1d,nqp_v,4), intent(in)     :: w3_basis_face
    real(kind=r_def), dimension(1,ndf_wtheta,nqp_h_1d,nqp_v,4), intent(in) :: wtheta_basis_face

    real(kind=r_def), dimension(ndf_w2,ndf_wtheta,ncell_3d), intent(inout) :: projection
    real(kind=r_def), dimension(undf_w3),                    intent(in)    :: exner
    real(kind=r_def),                                        intent(in)    :: scalar
    real(kind=r_def), dimension(nqp_v),                      intent(in)    :: wqp_v

    integer(kind=i_def), intent(in) :: adjacent_face(:)
    real(kind=r_def),    intent(in) :: out_face_normal(:,:)

    ! Internal variables
    integer(kind=i_def)                      :: df, df0, df2, k, ik, face, face_next
    integer(kind=i_def)                      :: qp1, qp2
    real(kind=r_def), dimension(ndf_w3)      :: exner_e, exner_next_e

    real(kind=r_def)                         :: v(3), normal(3), integrand
    real(kind=r_def)                         :: exner_at_fquad, &
                                                exner_next_at_fquad, &
                                                exner_av

    do k = 0, nlayers - 1
      ik = k + 1 + (cell-1)*nlayers
      do face = 1, size( adjacent_face, 1 )

        face_next = adjacent_face(face)

        do df = 1,ndf_w3
          exner_e(df)      = exner(stencil_w3_map(df, 1) + k)
          exner_next_e(df) = exner(stencil_w3_map(df, face+1) + k)
        end do

        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h_1d
            exner_at_fquad      = 0.0_r_def
            exner_next_at_fquad = 0.0_r_def
            do df = 1, ndf_w3
              exner_at_fquad      = exner_at_fquad + exner_e(df)*w3_basis_face(1,df,qp1,qp2,face)
              exner_next_at_fquad = exner_next_at_fquad + exner_next_e(df)*w3_basis_face(1,df,qp1,qp2,face_next)
            end do
            exner_av = 0.5_r_def*(exner_at_fquad + exner_next_at_fquad)

            do df0 = 1, ndf_wtheta
              normal =  out_face_normal(:,face)*wtheta_basis_face(1,df0,qp1,qp2,face)
              do df2 = 1, ndf_w2
                v  = w2_basis_face(:,df2,qp1,qp2,face)

                integrand = wqp_v(qp1)*wqp_v(qp2)*exner_av*dot_product(v, normal)
                projection(df2,df0,ik) = projection(df2,df0,ik) - scalar*integrand
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine weighted_proj_2theta_bd_code

end module weighted_proj_2theta_bd_kernel_mod
