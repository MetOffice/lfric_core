!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the coefficients for reconstructing a
!>        a  1D horizontal upwind polynomial representation of a tracer field on the
!>        faces of a cell
!> @details Compute the coefficients of the advective update of a tracer field using a high order
!>          1D polynomial fit to the integrated tracer values over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A symmetric polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + x^2
!>          This is exactly fitted over all cells in the stencil
!>          The methodology is inspired by that of Thuburn et.al GMD 2014 for
!>          2D reconstructions
!>          This method is only valid for lowest order elements
module poly1d_advective_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type, mesh_data_type,  &
                              GH_FIELD, GH_INTEGER,                 &
                              GH_WRITE, GH_READ,                    &
                              ANY_SPACE_1,                          &
                              GH_BASIS, CELLS, GH_QUADRATURE_XYoZ,  &
                              GH_QUADRATURE_edge, STENCIL, CROSS
use constants_mod,     only : r_def, i_def, l_def
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly1d_advective_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, Wtheta),                         &
       arg_type(GH_FIELD,   GH_READ,  Wtheta,      STENCIL(CROSS)),    &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_1, STENCIL(CROSS)),    &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape(2)  = (/ GH_QUADRATURE_XYoZ, GH_QUADRATURE_edge /)
contains
  procedure, nopass :: poly1d_advective_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_advective_coeffs_code
contains

!>@brief Compute the coefficients needed for a 1D horizontal reconstruction
!>       of a tracer field on horizontal faces
!>@param[in] nlayers Number of vertical layers
!>@param[out] coeff Array of fields to store the coefficients for the polynomial
!!                  reconstruction
!>@param[in]  mdwt Mass matrix diagonal for the Wtheta space, this is used to give
!!                 the cell volume
!>@param[in] chi1 1st component of the physical coordinate field
!>@param[in] chi2 2nd component of the physical coordinate field
!>@param[in] chi3 3rd component of the physical coordinate field
!>@param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!>@param[in] undf_wt Total number of degrees of freedom for Wtheta
!>@param[in] stencil_size_wt Number of cells in the Wtheta stencil
!>@param[in] smap_wt Stencil dofmap of the Wtheta stencil
!>@param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!>@param[in] undf_wx Total number of degrees of freedom for the coordinate space
!>@param[in] stencil_size_wx Number of cells in the coordinate space stencil
!>@param[in] smap_wx Stencil dofmap of the coordinate space stencil
!>@param[in] basis_wx Basis function of the coordinate space evaluated on
!!                    quadrature points. The vertical aspect must be on GLL points
!>@param[in] edge_basis_wx Basis function of the coordinate space evaluated on
!!                         quadrature points on the horizontal edges
!>@param[in] order Polynomial order for flux computations
!>@param[in] nfaces_h Number of horizontal neighbours
!>@param[in] nqp_h Number of horizontal quadrature points
!>@param[in] nqp_v Number of vertical quadrature points
!>@param[in] wqp_h Weights of horizontal quadrature points
!>@param[in] wqp_v Weights of vertical quadrature points
!>@param[in] n_edges Number of edges in the quadrature rule
!>@param[in] nqp_e Number of edge quadrature points
!>@param[in] wqp_e Weights of edge quadrature points
subroutine poly1d_advective_coeffs_code(nlayers,                   &
                                        coeff,                     &
                                        mdwt,                      &
                                        chi1, chi2, chi3,          &
                                        ndf_wt,                    &
                                        undf_wt,                   &
                                        stencil_size_wt,           &
                                        smap_wt,                   &
                                        ndf_wx,                    &
                                        undf_wx,                   &
                                        stencil_size_wx,           &
                                        smap_wx,                   &
                                        basis_wx,                  &
                                        edge_basis_wx,             &
                                        order,                     &
                                        nfaces_h,                  &
                                        nqp_h, nqp_v, wqp_h, wqp_v,&
                                        n_edges, nqp_e, wqp_e )

  use matrix_invert_mod,         only: matrix_invert
  use cross_product_mod,         only: cross_product
  use base_mesh_config_mod,      only: geometry, &
                                       geometry_spherical
  use poly_helper_functions_mod, only: local_distance_1d

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: order
  integer(kind=i_def), intent(in) :: nfaces_h
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, &
                                     ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: nqp_v, nqp_h, nqp_e, n_edges
  integer(kind=i_def), intent(in) :: stencil_size_wt, stencil_size_wx

  integer(kind=i_def), dimension(ndf_wt,stencil_size_wt), intent(in) :: smap_wt
  integer(kind=i_def), dimension(ndf_wx,stencil_size_wx), intent(in) :: smap_wx

  real(kind=r_def), dimension(undf_wt),                    intent(in)  :: mdwt
  real(kind=r_def), dimension(undf_wx),                    intent(in)  :: chi1, chi2, chi3
  real(kind=r_def), dimension(order+1, nfaces_h, undf_wt), intent(out) :: coeff

  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v),   intent(in) :: basis_wx
  real(kind=r_def), dimension(1,ndf_wx,nqp_e,n_edges), intent(in) :: edge_basis_wx

  real(kind=r_def), dimension(nqp_h),         intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v),         intent(in) ::  wqp_v
  real(kind=r_def), dimension(nqp_e,n_edges), intent(in) ::  wqp_e

  ! Local variables
  logical(kind=l_def) :: spherical
  integer(kind=i_def) :: ispherical
  integer(kind=i_def) :: k, ijk, df, qf0, stencil, nmonomial, qp, &
                         m, edge, stencil_depth, depth, face_mod
  integer(kind=i_def),           dimension(order+1,nfaces_h) :: map1d
  integer(kind=i_def),           dimension(0:nlayers)        :: kx, vert_face, qv
  real(kind=r_def)                                           :: xx, fn
  real(kind=r_def),              dimension(3)                :: x0, x1, xq, xn1
  real(kind=r_def), allocatable, dimension(:,:)              :: int_monomial, inv_int_monomial
  real(kind=r_def),              dimension(order+1)          :: beta, delta, monomial
  real(kind=r_def),              dimension(order+1)          :: area

  if ( geometry == geometry_spherical ) then
    spherical = .true.
    ispherical = 1_i_def
  else
    spherical = .false.
    ispherical = 0_i_def
  end if

  ! Compute 1d map from the cross stencil
  ! i.e for order = 2 the stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |
  ! so map1d is
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! First cell is always the centre cell
  stencil_depth = order/2
  map1d(1,:) = 1
  do edge = 1,nfaces_h
    depth=1
    face_mod = mod(edge+1,2) * stencil_depth
    do stencil = 2,stencil_depth+1
      map1d(stencil+depth-1, edge) = stencil + face_mod
      map1d(stencil+depth, edge) = stencil + order + face_mod
      depth=depth+1
    end do
  end do

  ! Avoid compile warning for unused variables
  fn = wqp_v(1)
  fn = wqp_e(1,1)

  ! Number of monomials to use (all polynomials up to total degree of order)
  nmonomial = (order + 1)

  ! Index of quadrature point in the centre of the vertical face
  ! (this is only true if the number of quadrature points is odd
  qf0 = (nqp_h+1)/2

  ! Step 1: Build integrals of monomials over all cells in advection stencils
  ! Initialize to zero
  allocate( int_monomial(nmonomial, nmonomial),  &
            inv_int_monomial(nmonomial, nmonomial) )

  vert_face = 0
  vert_face(nlayers) = 4
  ! k for the coordinate field: ensure it never equals nlayers
  kx = (/(k, k=0, nlayers)/)
  kx(nlayers) = nlayers-1
  ! Quadrature point for the coordinate field
  qv = 1
  qv(nlayers) = nqp_v

  ! Loop over all layers: goes to nalyers to pick up top point
  layer_loop: do k = 0, nlayers
    ! Position vector of bottom of this cell unless very last point in
    ! which case use the top of the cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = smap_wx(df, 1) + kx(k)
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qf0,qv(k))
    end do

    ! Compute the coefficients of each cell in the stencil for
    ! each edge when this is the upwind cell
    edge_loop: do edge = 1,nfaces_h
      int_monomial = 0.0_r_def

      ! Find direction of first neighbour to establish axes of
      ! Local coordinate system
      ! Position vector of neighbour cell centr
      x1 = 0.0_r_def
      do df = 1, ndf_wx
        ijk = smap_wx(df, edge+1) + kx(k)
        x1(:) = x1(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qf0,qv(k))
      end do

      ! Unit normal to plane containing points 0 and 1
      ! cross_product is zero if x0(3) = x1(3) = 0, which occurs for the first level in
      ! cartesian domains so use the direct formulation
      xn1 = ispherical*cross_product(x0,x1) &
          + (1_i_def-ispherical)*(/ x0(2) - x1(2), x0(1) - x1(1), x0(1)*x1(2) - x0(2)*x1(1) /)

      xn1 = xn1/sqrt(xn1(1)**2 + xn1(2)**2 + xn1(3)**2)

      ! Loop over all cells in the stencil
      stencil_loop: do stencil = 1, order+1
        area(stencil) = mdwt(smap_wt( 1, map1d(stencil,edge)) + k)
        if ( k > 0 .and. k < nlayers ) &
            area(stencil) =  area(stencil) &
            + mdwt(smap_wt( 2, map1d(stencil,edge)) + k-1)
        ! Integrate monomials over this cell
        quadrature_loop: do qp = 1, nqp_h
          ! First: Compute physical coordinate of each quadrature point
          xq = 0.0_r_def
          do df = 1, ndf_wx
            ijk = smap_wx(df, map1d(stencil,edge)) + kx(k)
            xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qp,qv(k))
          end do
          xq(3) = ispherical*xq(3) + (1_i_def-ispherical)*x0(3)
          ! Second: Compute the local coordinate of each quadrature point from the
          !         physical coordinate
          xx = local_distance_1d(x0, xq, xn1, spherical)
          ! Third: Compute each needed monomial in terms of the local coordinate
          !        on each quadrature point
          ! Loop over monomials
          do m = 1, nmonomial
            fn = xx**(m-1)
            int_monomial(stencil,m) = int_monomial(stencil,m) &
                                    + wqp_h(qp)*fn*area(stencil)
          end do
        end do quadrature_loop
      end do stencil_loop

      ! Manipulate the integrals of monomials
      call matrix_invert(int_monomial, inv_int_monomial, nmonomial)

      ! Initialise polynomial coeffficients to zero
      coeff(:,edge,smap_wt(1,1)+k) = 0.0_r_def
      ! Loop over quadrature points on this edge
      edge_quadrature_loop: do qp = 1,nqp_e

        ! Obtain physical coordinates of gauss points on this edge
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = smap_wx(df, 1) + kx(k)
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*edge_basis_wx(1,df,qp,edge+vert_face(k))
        end do

        ! Obtain local coordinates of gauss points on this edge
        xx = local_distance_1d(x0, xq, xn1, spherical)

        ! Evaluate polynomial fit
        ! Loop over monomials
        do stencil = 1, order+1
          monomial(stencil) = xx**(stencil-1)
        end do
        do stencil = 1, order+1
          delta(:) = 0.0_r_def
          delta(stencil) = 1.0_r_def
          beta = matmul(inv_int_monomial,delta)
          coeff(stencil,edge,smap_wt(1,1)+k) = dot_product(monomial,beta)*area(stencil)
        end do
      end do edge_quadrature_loop
    end do edge_loop
  end do layer_loop
  deallocate( int_monomial, inv_int_monomial )

end subroutine poly1d_advective_coeffs_code

end module poly1d_advective_coeffs_kernel_mod
