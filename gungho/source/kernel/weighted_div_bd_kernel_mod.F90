!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes boundary part of the weighted divergence operator

!> @details The kernel computes the boundary part of the weighted divergence operator
!>         This consists of \f[<\sigma,\theta*\mathbf{v}\cdot\mathbf{n}> \f]
!>         where sigma is the W3 test function, v is the W2 trial function,
!>         theta is the potential temperature, and \mathbf{n} is the outward pointing normal.
!>         Each face provides two contributions to a pressure point, one from
!>         the right side, using the potential temperature on the right side of
!>         the face and one from the left side, using the potential temperature on the left 
!>         side of the face

module weighted_div_bd_kernel_mod
  use kernel_mod,              only : kernel_type
  use argument_mod,            only : arg_type, func_type, mesh_data_type, &
                                      GH_OPERATOR, GH_FIELD, GH_REAL,      &
                                      GH_READ, GH_INC,                     &
                                      W2, W3, Wtheta, GH_BASIS,            &
                                      CELLS, GH_QUADRATURE_XYoZ,           &
                                      adjacent_face,                       &
                               reference_element_number_horizontal_faces,  &
                                      reference_element_out_face_normal
  use constants_mod,           only : r_def, i_def


  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: weighted_div_bd_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                               &
       arg_type(GH_OPERATOR, GH_INC,  W2, W3),                        &
       arg_type(GH_FIELD,    GH_READ, Wtheta),                        &
       arg_type(GH_REAL,     GH_READ)                                 &
      /)
    type(func_type) :: meta_funcs(3) = (/                             &
       func_type(W2, GH_BASIS),                                       &
       func_type(W3, GH_BASIS),                                       &
       func_type(Wtheta, GH_BASIS)                                    &
      /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
    type(mesh_data_type) :: meta_init(3) = (/                        &
        mesh_data_type( adjacent_face ),                             &
        mesh_data_type( reference_element_number_horizontal_faces ), &
        mesh_data_type( reference_element_out_face_normal )          &
      /)
  contains
    procedure, nopass :: weighted_div_bd_code
  end type

  !-------------------------------------------------------------------------------
  ! Constructors
  !-------------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface weighted_div_bd_kernel_type
    module procedure weighted_div_bd_kernel_constructor
  end interface

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public weighted_div_bd_code

contains

  type(weighted_div_bd_kernel_type) function weighted_div_bd_kernel_constructor() result(self)
    return
  end function weighted_div_bd_kernel_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Computes the boundary terms in the weighted divergence for the
  !!        Helmholtz lhs
  !!
  !! @param[in] cell Cell number
  !! @param[in] nlayers Number of layers
  !! @param[in] ncell_3d ncell*ndf
  !! @param[in] div Local stencil of the div operator
  !! @param[in] theta Potential temperature
  !! @param[in] scalar Real to scale matrix by
  !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
  !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
  !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
  !! @param[in] undf_wtheta Number unique of degrees of freedom  for wtheta
  !! @param[in] stencil_wtheta_map W2 dofmaps for the stencil
  !! @param[in] stencil_wtheta_size Size of the W2 stencil (number of cells)
  !! @param[in] nqp_v Number of quadrature points in the vertical
  !! @param[in] nqp_h_1d  Number of quadrature points in a single horizontal
  !!                      direction.
  !! @param[in] wqp_v Vertical quadrature weights
  !! @param[in] w2_basis_face  Basis functions evaluated at gaussian &
  !!                           quadrature points on horizontal faces.
  !! @param[in] w3_basis_face  Basis functions evaluated at gaussian &
  !!                           quadrature points on horizontal faces.
  !! @param[in] wtheta_basis_face  Basis functions evaluated at gaussian &
  !!                               quadrature points on horizontal faces.
  !! @param[in] adjacent_face  Vector containing information on neighbouring &
  !!                           face index for the current cell.
  !! @param[in] number_horizontal_faces  Number of reference element faces
  !!                                     bisected by a horizontal plane.
  !! @param[in] out_face_normal  Vectors normal to the faces of the &
  !!                             reference element.
  !!
  !!
  subroutine weighted_div_bd_code( cell, nlayers, ncell_3d, &
                                   div, theta, scalar,      &
                                   ndf_w2, ndf_w3,          &
                                   ndf_wtheta, undf_wtheta, &
                                   stencil_wtheta_map,      &
                                   stencil_wtheta_size,     &
                                   nqp_v, nqp_h_1d, wqp_v,  &
                                   w2_basis_face,           &
                                   w3_basis_face,           &
                                   wtheta_basis_face,       &
                                   adjacent_face,           &
                                   number_horizontal_faces, &
                                   out_face_normal)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, nlayers, nqp_v, nqp_h_1d, ncell_3d
    integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in) :: stencil_wtheta_size

    integer(kind=i_def), dimension(ndf_wtheta,stencil_wtheta_size),  intent(in) :: stencil_wtheta_map

    real(kind=r_def), dimension(3,ndf_w2,nqp_h_1d,nqp_v,4), intent(in)     :: w2_basis_face
    real(kind=r_def), dimension(1,ndf_w3,nqp_h_1d,nqp_v,4), intent(in)     :: w3_basis_face
    real(kind=r_def), dimension(1,ndf_wtheta,nqp_h_1d,nqp_v,4), intent(in) :: wtheta_basis_face

    real(kind=r_def), dimension(ndf_w2,ndf_w3,ncell_3d), intent(inout) :: div
    real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
    real(kind=r_def),                         intent(in)    :: scalar

    real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

    integer(i_def), intent(in) :: adjacent_face(number_horizontal_faces)
    integer(i_def), intent(in) :: number_horizontal_faces
    real(r_def),    intent(in) :: out_face_normal(:,:)

    ! Internal variables
    integer(kind=i_def)              :: df, df2, df3, k, ik, face, face_next
    integer(kind=i_def)              :: qp1, qp2

    real(kind=r_def), dimension(ndf_wtheta) :: theta_e, theta_next_e

    real(kind=r_def) :: v_dot_n, integrand
    real(kind=r_def) :: theta_at_fquad, theta_next_at_fquad
    real(kind=r_def) :: this_bd_term, next_bd_term


    do k = 0, nlayers - 1
      ik = k + 1 + (cell-1)*nlayers
      do face = 1, number_horizontal_faces

        face_next = adjacent_face(face)

        do df = 1,ndf_wtheta
          theta_e(df)      = theta(stencil_wtheta_map(df, 1)      + k)
          theta_next_e(df) = theta(stencil_wtheta_map(df, face+1) + k)
        end do
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h_1d
            theta_at_fquad      = 0.0_r_def 
            theta_next_at_fquad = 0.0_r_def               
            do df = 1, ndf_wtheta
              theta_at_fquad       = theta_at_fquad      + theta_e(df)     *wtheta_basis_face(1,df,qp1,qp2,face)
              theta_next_at_fquad  = theta_next_at_fquad + theta_next_e(df)*wtheta_basis_face(1,df,qp1,qp2,face_next)
            end do                
            do df3 = 1, ndf_w3
              do df2 = 1, ndf_w2
                v_dot_n  = dot_product(w2_basis_face(:,df2,qp1,qp2,face),out_face_normal(:,face))
                this_bd_term =  v_dot_n*theta_at_fquad
                next_bd_term = -v_dot_n*theta_next_at_fquad
                integrand = wqp_v(qp1)*wqp_v(qp2)*w3_basis_face(1,df3,qp1,qp2,face) &
                             * 0.5_r_def*(this_bd_term + next_bd_term)
                div(df2,df3,ik) = div(df2,df3,ik) - scalar*integrand
              end do ! df2
            end do ! df3
          end do ! qp1
        end do ! qp2
      end do ! faces
    end do ! layers

  end subroutine weighted_div_bd_code

 end module weighted_div_bd_kernel_mod
