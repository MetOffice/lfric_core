!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Projects the components of a vector field into a scalar space.
!>
module sci_gp_vector_rhs_kernel_mod

  use argument_mod,              only : arg_type, func_type,       &
                                        GH_FIELD, GH_REAL, GH_INC, &
                                        GH_READ, ANY_SPACE_1,      &
                                        ANY_SPACE_2, ANY_SPACE_9,  &
                                        ANY_DISCONTINUOUS_SPACE_3, &
                                        GH_BASIS, GH_DIFF_BASIS,   &
                                        CELL_COLUMN, GH_QUADRATURE_XYoZ
  use base_mesh_config_mod,      only : geometry,                   &
                                        geometry_spherical
  use sci_chi_transform_mod,     only : chi2xyz
  use constants_mod,             only : r_def, i_def
  use sci_coordinate_jacobian_mod, only : coordinate_jacobian,        &
                                          coordinate_jacobian_inverse
  use coord_transform_mod,       only : cart2sphere_vector
  use fs_continuity_mod,         only : W0, W2
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: gp_vector_rhs_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
         arg_type(GH_FIELD*3, GH_REAL, GH_INC,  ANY_SPACE_1),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_2),               &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2)                         &
         /)
    type(func_type) :: meta_funcs(3) = (/                                   &
         func_type(ANY_SPACE_1, GH_BASIS),                                  &
         func_type(ANY_SPACE_2, GH_BASIS),                                  &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: gp_vector_rhs_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: gp_vector_rhs_code

contains

!> @brief Computes the right-hand-side of the Galerkin projection for a vector
!> field into a scalar space by decomposing the vector into orthogonal
!> components in cartesian or (lon,lat,r) coordinates.
!> @details Computes rhs_i = int (gamma * f_i dx) for a vector field f which  is
!>          decomposed into orthogonal components and a separate right hand side
!>          field is computed for each component, this allows a vector field to
!>          be projected into three separate scalar fields suitable for further
!>          manipulation
!! @param[in] nlayers   Number of layers
!! @param[in,out] rhs1  Field to compute
!! @param[in,out] rhs2  Field to compute
!! @param[in,out] rhs3  Field to compute
!! @param[in] field     Field to be projected
!! @param[in] chi_1     1st coordinate field
!! @param[in] chi_2     2nd coordinate field
!! @param[in] chi_3     3rd coordinate field
!! @param[in] panel_id  Field giving the ID for mesh panels.
!! @param[in] w2_field  W2_field needed to get function space components
!! @param[in] ndf       Number of degrees of freedom per cell
!! @param[in] undf      Number of degrees of freedom
!! @param[in] map       Dofmap for the cell at the base of the column
!! @param[in] basis     Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_f     Number of degrees of freedom per cell for the field to be projected
!! @param[in] undf_f    Number of degrees of freedom for the field to be projected
!! @param[in] map_f     Dofmap for the cell at the base of the column
!! @param[in] f_basis   Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi   Number of degrees of freedom per cell for chi
!! @param[in] undf_chi  Number of unique degrees of freedom for chi
!! @param[in] map_chi   Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions for Wchi evaluated at
!!                      gaussian quadrature points
!! @param[in] chi_diff_basis Differential of the Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] ndf_w2  Number of degrees of freedom per cell for the w2 field
!! @param[in] undf_w2 Total number of degrees of freedom for the w2 field
!! @param[in] map_w2  Dofmap for the cell at the base of the column
!! @param[in] nqp_h   Number of quadrature points in the horizontal
!! @param[in] nqp_v   Number of quadrature points in the vertical
!! @param[in] wqp_h   Horizontal quadrature weights
!! @param[in] wqp_v   Vertical quadrature weights
subroutine gp_vector_rhs_code(nlayers,                           &
                              rhs1, rhs2, rhs3, field,           &
                              chi_1, chi_2, chi_3,               &
                              panel_id, w2_field,                &
                              ndf, undf, map, basis,             &
                              ndf_f, undf_f, map_f, f_basis,     &
                              ndf_chi, undf_chi,                 &
                              map_chi,                           &
                              chi_basis, chi_diff_basis,         &
                              ndf_pid, undf_pid, map_pid,        &
                              ndf_w2, undf_w2, map_w2,           &
                              nqp_h, nqp_v, wqp_h, wqp_v         &
                             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_f, ndf_w2, ndf_pid
  integer(kind=i_def), intent(in) :: undf, undf_f, undf_w2, undf_pid
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf),         intent(in) :: map
  integer(kind=i_def), dimension(ndf_f),       intent(in) :: map_f
  integer(kind=i_def), dimension(ndf_chi),     intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w2),      intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid


  real(kind=r_def), intent(in), dimension(1,ndf,    nqp_h,nqp_v) :: basis
  real(kind=r_def), intent(in), dimension(3,ndf_f,  nqp_h,nqp_v) :: f_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),      intent(inout) :: rhs1, rhs2, rhs3
  real(kind=r_def), dimension(undf_chi),     intent(in) :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_f),       intent(in) :: field
  real(kind=r_def), dimension(undf_w2),      intent(in) :: w2_field
  real(kind=r_def), dimension(undf_pid),     intent(in) :: panel_id
  real(kind=r_def), dimension(nqp_h),        intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v),        intent(in) ::  wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian, jacobian_inv
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: u_at_quad, x_at_quad, u_physical, coords
  real(kind=r_def)                             :: integrand
  logical                                      :: hdiv
  integer(kind=i_def)                          :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Check if this is hdiv (W2) field or a hcurl (W1) field
  if ( ndf_f == ndf_w2 ) then
    hdiv = .true.
  else
    hdiv = .false.
  end if

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k )
      chi_2_cell(df) = chi_2( map_chi(df) + k )
      chi_3_cell(df) = chi_3( map_chi(df) + k )
    end do
    call coordinate_jacobian(ndf_chi,        &
                             nqp_h,          &
                             nqp_v,          &
                             chi_1_cell,     &
                             chi_2_cell,     &
                             chi_3_cell,     &
                             ipanel,         &
                             chi_basis,      &
                             chi_diff_basis, &
                             jacobian,       &
                             dj              )

    if ( .not. hdiv) call coordinate_jacobian_inverse(nqp_h, nqp_v, &
                                                      jacobian, dj, &
                                                      jacobian_inv)

    do df = 1, ndf
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          ! Compute vector in computational space
          u_at_quad(:) = 0.0_r_def
          do df2 = 1,ndf_f
            u_at_quad(:) = u_at_quad(:) &
                         + f_basis(:,df2,qp1,qp2)*field(map_f(df2) + k)
          end do
          if ( hdiv ) then
            ! For W2 space
            u_at_quad(:) = matmul(jacobian(:,:,qp1,qp2),u_at_quad(:))/dj(qp1,qp2)
          else
            ! For W1 space
            u_at_quad(:) = matmul(transpose(jacobian_inv(:,:,qp1,qp2)),u_at_quad(:))
          end if
          ! Compute physical coordinate of quadrature point
          if ( geometry == geometry_spherical ) then
            coords(:) = 0.0_r_def
            do df2 = 1,ndf_chi
              coords(1) = coords(1) + chi_1_cell(df2)*chi_basis(1,df2,qp1,qp2)
              coords(2) = coords(2) + chi_2_cell(df2)*chi_basis(1,df2,qp1,qp2)
              coords(3) = coords(3) + chi_3_cell(df2)*chi_basis(1,df2,qp1,qp2)
            end do

            ! Obtain (X,Y,Z) coordinates for converting components of u
            call chi2xyz(coords(1), coords(2), coords(3), ipanel, &
                         x_at_quad(1), x_at_quad(2), x_at_quad(3))

            u_physical(:) = cart2sphere_vector(x_at_quad, u_at_quad)
          else
            u_physical(:) = u_at_quad(:)
          end if
          integrand = wqp_h(qp1)*wqp_v(qp2)*basis(1,df,qp1,qp2)*dj(qp1,qp2)
          rhs1(map(df) + k) = rhs1(map(df) + k) + integrand * u_physical(1)
          rhs2(map(df) + k) = rhs2(map(df) + k) + integrand * u_physical(2)
          rhs3(map(df) + k) = rhs3(map(df) + k) + integrand * u_physical(3)
        end do
      end do
    end do
  end do

end subroutine gp_vector_rhs_code

end module sci_gp_vector_rhs_kernel_mod
