!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the operators for rhs for the mapping of the W2 physical
!> wind field from the finite difference W3 lat, W3 lon and Wtheta
!> up computational wind components
!>
!> The kernel computes the linear operators to evaluate
!>  a very crude approximation to the rhs of the equation u = u0
!> where u0 is the physical wind field. These are used to project
!> the computational wind field onto the physical wind field
!> using Galerkin projection.
!>
module sci_compute_map_u_operators_kernel_mod

  use argument_mod,            only : arg_type, func_type,       &
                                      GH_FIELD, GH_REAL,         &
                                      GH_OPERATOR,               &
                                      GH_INC, GH_READ, GH_WRITE, &
                                      ANY_SPACE_9,               &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2, W3, Wtheta
  use kernel_mod,              only : kernel_type
  use log_mod,                 only : log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_map_u_operators_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W3),                 &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W3),                 &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, WTHETA),             &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(4) = (/                                  &
         func_type(W2,          GH_BASIS),                                 &
         func_type(W3,          GH_BASIS),                                 &
         func_type(WTHETA,      GH_BASIS),                                 &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_map_u_operators_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: compute_map_u_operators_code

contains

!> @brief Compute the operators to evaluate
!>        the right hand side in order to to map the wind field.
!! @param[in] cell     Cell number
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d_1 ncell*ndf
!! @param[in,out] u_lon_op Operator to map u_lon from W3 to W2
!! @param[in] ncell_3d_2 ncell*ndf
!! @param[in,out] u_lat_op Operator to map u_lat from W3 to W2
!! @param[in] ncell_3d_3 ncell*ndf
!! @param[in,out] u_up_op Operator to map u_up from WTHETA to W2
!! @param[in] chi_sph_1 1st coordinate in spherical Wchi
!! @param[in] chi_sph_2 2nd coordinate in spherical Wchi
!! @param[in] chi_sph_3 3rd coordinate in spherical Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] basis_w2 W2 basis functions evaluated at quadrature points
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] basis_w3 W3 basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_wt Number of degrees of freedom per cell for WTHETA
!! @param[in] basis_wt WTHETA basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi_sph Number of degrees of freedom per cell for spherical chi
!! @param[in] undf_chi_sph Number of unique degrees of freedom for spherical chi
!! @param[in] map_chi_sph Dofmap for the cell at the base of the column for spherical chi
!! @param[in] chi_sph_basis Basis functions for spherical Wchi evaluated at
!!                          gaussian quadrature points
!! @param[in] chi_sph_diff_basis Differential of the spherical Wchi basis functions
!!                               evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_map_u_operators_code(cell, nlayers, ncell_3d_1, &
                        u_lon_op, ncell_3d_2, u_lat_op,            &
                        ncell_3d_3, u_up_op,                       &
                        chi_sph_1, chi_sph_2, chi_sph_3, panel_id, &
                        ndf_w2, basis_w2,                          &
                        ndf_w3, basis_w3,                          &
                        ndf_wt, basis_wt,                          &
                        ndf_chi_sph, undf_chi_sph, map_chi_sph,    &
                        chi_sph_basis, chi_sph_diff_basis,         &
                        ndf_pid, undf_pid, map_pid,                &
                        nqp_h, nqp_v, wqp_h, wqp_v                 &
                        )

  use sci_chi_transform_mod,       only : chi2llr
  use sci_coordinate_jacobian_mod, only : coordinate_jacobian
  use coord_transform_mod,         only : sphere2cart_vector

  use finite_element_config_mod, only: coord_system
  use base_mesh_config_mod,      only: geometry, topology, &
                                       geometry_spherical, &
                                       geometry_planar
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nlayers
  integer(kind=i_def), intent(in) :: ncell_3d_1, ncell_3d_2, ncell_3d_3
  integer(kind=i_def), intent(in) :: ndf_chi_sph, ndf_w2, ndf_w3, ndf_wt, ndf_pid
  integer(kind=i_def), intent(in) :: undf_pid
  integer(kind=i_def), intent(in) :: undf_chi_sph
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf_chi_sph), intent(in) :: map_chi_sph
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(3,ndf_w2, nqp_h,nqp_v)     :: basis_w2
  real(kind=r_def), intent(in), dimension(1,ndf_w3, nqp_h,nqp_v)     :: basis_w3
  real(kind=r_def), intent(in), dimension(1,ndf_wt, nqp_h,nqp_v)     :: basis_wt
  real(kind=r_def), intent(in), dimension(1,ndf_chi_sph,nqp_h,nqp_v) :: chi_sph_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi_sph,nqp_h,nqp_v) :: chi_sph_diff_basis

  real(kind=r_def), dimension(ncell_3d_1,ndf_w2,ndf_w3), intent(inout) :: u_lon_op
  real(kind=r_def), dimension(ncell_3d_2,ndf_w2,ndf_w3), intent(inout) :: u_lat_op
  real(kind=r_def), dimension(ncell_3d_3,ndf_w2,ndf_wt), intent(inout) :: u_up_op
  real(kind=r_def), dimension(undf_pid),     intent(in) :: panel_id
  real(kind=r_def), dimension(undf_chi_sph), intent(in) :: chi_sph_1, chi_sph_2, chi_sph_3

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, df3, dft
  integer(kind=i_def)                          :: k, ik, qp1, qp2, ipanel
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi_sph)     :: chi_sph_1_cell, chi_sph_2_cell, chi_sph_3_cell
  real(kind=r_def), dimension(3)               :: u_lon_phys_basis, u_lon_sph_basis, &
                                                  u_lat_phys_basis, u_lat_sph_basis, &
                                                  u_up_phys_basis, u_up_sph_basis,   &
                                                  coords, llr
  real(kind=r_def)                             :: integrand

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    ik = k + 1 + (cell-1)*nlayers
    u_lon_op(ik,:,:) = 0.0_r_def
    u_lat_op(ik,:,:) = 0.0_r_def
    u_up_op(ik,:,:) = 0.0_r_def
    do df = 1, ndf_chi_sph
      chi_sph_1_cell(df) = chi_sph_1( map_chi_sph(df) + k )
      chi_sph_2_cell(df) = chi_sph_2( map_chi_sph(df) + k )
      chi_sph_3_cell(df) = chi_sph_3( map_chi_sph(df) + k )
    end do

    call coordinate_jacobian(coord_system,       &
                             geometry,           &
                             topology,           &
                             scaled_radius,      &
                             ndf_chi_sph,        &
                             nqp_h,              &
                             nqp_v,              &
                             chi_sph_1_cell,     &
                             chi_sph_2_cell,     &
                             chi_sph_3_cell,     &
                             ipanel,             &
                             chi_sph_basis,      &
                             chi_sph_diff_basis, &
                             jacobian,           &
                             dj)

      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h

          ! Compute analytical vector wind in physical space
          if ( geometry == geometry_spherical ) then
          ! Need position vector for obtaining (X,Y,Z) components of the physical u
          ! basis functions
            coords(:) = 0.0_r_def
            do df = 1, ndf_chi_sph
              coords(1) = coords(1) + chi_sph_1_cell(df)*chi_sph_basis(1,df,qp1,qp2)
              coords(2) = coords(2) + chi_sph_2_cell(df)*chi_sph_basis(1,df,qp1,qp2)
              coords(3) = coords(3) + chi_sph_3_cell(df)*chi_sph_basis(1,df,qp1,qp2)
            end do

            llr(:) = 0.0_r_def

            call chi2llr(coords(1), coords(2), coords(3), &
                          ipanel, llr(1), llr(2), llr(3))
          end if


          do df3 = 1, ndf_w3
            do df2 = 1, ndf_w2

              ! Compute analytical vector wind basis functions in physical space
              if ( geometry == geometry_spherical ) then
                u_lon_sph_basis(:) = 0.0_r_def
                u_lat_sph_basis(:) = 0.0_r_def
                u_lon_sph_basis(1) = basis_w3(1,df3,qp1,qp2)
                u_lat_sph_basis(2) = basis_w3(1,df3,qp1,qp2)
                u_lon_phys_basis   = sphere2cart_vector(u_lon_sph_basis,llr)
                u_lat_phys_basis   = sphere2cart_vector(u_lat_sph_basis,llr)

              else if ( geometry == geometry_planar ) then
                u_lon_phys_basis(:) = 0.0_r_def
                u_lat_phys_basis(:) = 0.0_r_def
                u_lon_phys_basis(1) = basis_w3(1,df3,qp1,qp2)
                u_lat_phys_basis(2) = basis_w3(1,df3,qp1,qp2)
              else

                call log_event('compute_map_u_operators_kernel is not implemented ' //    &
                                'with your geometry',       &
                                LOG_LEVEL_ERROR)

              end if

              integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                             basis_w2(:,df2,qp1,qp2)),u_lon_phys_basis)
              u_lon_op(ik,df2,df3) = u_lon_op(ik,df2,df3) &
                                    + wqp_h(qp1)*wqp_v(qp2)*integrand
              integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                             basis_w2(:,df2,qp1,qp2)),u_lat_phys_basis)
              u_lat_op(ik,df2,df3) = u_lat_op(ik,df2,df3)&
                                    + wqp_h(qp1)*wqp_v(qp2)*integrand
           end do
         end do

         do dft = 1, ndf_wt
           do df2 = 1, ndf_w2
             ! Compute analytical vector wind basis functions in physical space
             if ( geometry == geometry_spherical ) then

              u_up_sph_basis(:) = 0.0_r_def
              u_up_sph_basis(3) = basis_wt(1,dft,qp1,qp2)
              u_up_phys_basis   = sphere2cart_vector(u_up_sph_basis,llr)

              else if ( geometry == geometry_planar ) then
                u_up_phys_basis(:) = 0.0_r_def
                u_up_phys_basis(3) = basis_wt(1,dft,qp1,qp2)
              else

                call log_event('compute_map_u_operators_kernel is not implemented ' //    &
                               'with your geometry',       &
                               LOG_LEVEL_ERROR)

              end if

              integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                             basis_w2(:,df2,qp1,qp2)),u_up_phys_basis)
              u_up_op(ik,df2,dft) = u_up_op(ik,df2,dft) &
                                    + wqp_h(qp1)*wqp_v(qp2)*integrand
          end do
        end do
      end do
    end do

  end do


end subroutine compute_map_u_operators_code

end module sci_compute_map_u_operators_kernel_mod
