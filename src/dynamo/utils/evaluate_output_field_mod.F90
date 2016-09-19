!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
module evaluate_output_field_mod

use base_mesh_config_mod,    only: geometry, &
                                   base_mesh_geometry_spherical
use constants_mod,           only: r_def, i_def
use coordinate_jacobian_mod, only: coordinate_jacobian, &
                                   coordinate_jacobian_inverse
use coord_transform_mod,     only: cartesian_distance, llr2xyz
use field_mod,               only: field_type, field_proxy_type
use log_mod,                 only: log_event, LOG_LEVEL_ERROR
use mesh_mod,                only: mesh_type
use planet_config_mod,       only: scaled_radius

implicit none

contains
!>@brief  Evaluates a field at a point (x_in) that horizontally lies within a 
!>        given cell
!>@details Subroutine that evaluates a field at a vertical column of given
!>         points.
!>         The column of points is taken to lie within a given column of grid
!>         cells. A newton method is used to compute the exact point within each
!>         cell where a point lies and then the field is evaluated at this point
!>         using the known basis functions.
!>@deprecated This is a temporary solution until a better output routine is
!>implemented as which point this routine will be reviewed to see if it will be
!>needed elsewhere in the model
!>@param[in]  field     Field object to evaluate
!>@param[in]  chi       3D coordinate field
!>@param[in]  x_in      Point to evaluate the field at
!>@param[in]  cell      Horizontal cell that x_in lies within
!>@param[in]  nz        Number of vertical points to evaluate
!>@param[out] field_out Array containing field evaluates at x_in
subroutine evaluate_output_field( field, chi, x_in, cell, nz, field_out )

  type(field_type), intent(in)  :: field, chi(3) 
  integer,          intent(in)  :: cell, nz
  real(kind=r_def), intent(in)  :: x_in(3,nz)
  real(kind=r_def), intent(out) :: field_out(nz)

  type(field_proxy_type)        :: chi_proxy(3), field_proxy

  integer, pointer              :: map(:)   => null(), &
                                   map_f(:) => null()
  integer                       :: iter, ndf, ndf_f, k, df, &
                                   dir, nlayers, dfk
  integer                       :: kk, alloc_error
  integer, allocatable          :: out_layer(:)
  integer,          parameter   :: NEWTON_ITERS = 4
  real(kind=r_def)              :: jac(3,3), jac_inv(3,3), dj(1,1), g_func(3), &
                                   gamma(1), x_loc(3), x_out(3), offset 
  real(kind=r_def)              :: domain_top, eta_in
  real(kind=r_def), allocatable :: eta(:)  ! eta(0:nlayers)
  real(kind=r_def), allocatable :: chi_cell(:,:), dgamma(:,:) 
  type(mesh_type),  pointer  :: mesh => null()

  offset = 0.0_r_def

  if ( geometry == base_mesh_geometry_spherical ) offset = scaled_radius

  chi_proxy(1) = chi(1)%get_proxy()
  chi_proxy(2) = chi(2)%get_proxy()
  chi_proxy(3) = chi(3)%get_proxy()
  field_proxy  = field%get_proxy()

  map => chi_proxy(1)%vspace%get_cell_dofmap(cell)
  ndf = chi_proxy(1)%vspace%get_ndf() 
  map_f => field_proxy%vspace%get_cell_dofmap(cell)
  ndf_f = field_proxy%vspace%get_ndf() 
  nlayers = field_proxy%vspace%get_nlayers()
  mesh => field_proxy%vspace%get_mesh()
  domain_top = mesh%get_domain_top()

  allocate ( chi_cell(ndf,3), dgamma(3,ndf), out_layer(nz) )
  allocate ( eta(0:nlayers), STAT = alloc_error )
  if ( alloc_error /= 0 ) then
    call log_event( " evaluate_output_field: Unable to allocate "// &
                    "local array eta(0:nlayers) ", LOG_LEVEL_ERROR )
  end if
  call mesh%get_eta(eta)

! Compute layer each output point is in using non-dimensional 
! eta_in = x_in(3,k)/domain_top (rather than dz like previously used)
  do k = 1,nz
    eta_in = x_in(3,k)/domain_top
    ! Initialise position layer to next-to-top layer of mesh
    out_layer(k) = nlayers-1
    do kk = 0,nlayers-1
      if ( eta_in >= eta(kk) .and. eta_in < eta(kk+1) ) &
           out_layer(k) = min(kk,nlayers-1)
    end do
  end do

! Find the horizontal coordinates (x_out in [0,1]^2) corresponding to each 
! input point using a newton method with  a fixed number of iterations
  do df = 1,ndf 
    chi_cell(df,1) = chi_proxy(1)%data( map(df) )
    chi_cell(df,2) = chi_proxy(2)%data( map(df) )
    chi_cell(df,3) = chi_proxy(3)%data( map(df) )
  end do
! First guess of out point in reference element
  x_out(:) = (/ 0.5_r_def, 0.5_r_def, 0.0_r_def /)
  x_loc(:) = x_in(:,1)
  if ( geometry == base_mesh_geometry_spherical ) call llr2xyz( x_in(1,1), &
                                   x_in(2,1), &
                                   x_in(3,1) + offset, &
                                   x_loc(1),  &
                                   x_loc(2),  &
                                   x_loc(3))

! Find location in computational space of point to evaluate field at
  do iter = 1,NEWTON_ITERS
    do df = 1,ndf
      dgamma(:,df) = chi_proxy(1)%vspace%evaluate_diff_basis(df, x_out)
    end do

    call coordinate_jacobian( ndf, &
                              1,   &
                              1,   &
                              chi_cell(:,1), &
                              chi_cell(:,2), &
                              chi_cell(:,3), &
                              dgamma, &
                              jac, &
                              dj)
    call coordinate_jacobian_inverse(1, 1, jac, dj, jac_inv)
! Compute g(xi^n) - [x,y,z]
    g_func = - x_loc(:)
    do df = 1, ndf
      gamma(:) = chi_proxy(1)%vspace%evaluate_basis(df, x_out)
      do dir = 1,3
        g_func(dir) = g_func(dir) + chi_cell(df,dir)*gamma(1)
      end do
    end do
    x_out(:) = x_out(:) - matmul(jac_inv, g_func)
    x_out(3) = 0.0_r_def ! fixed to surface
  end do
! Evaluate field at xi
  do k = 1,nz
  ! Find vertical xi point between 0 and 1 
    kk = out_layer(k)
    eta_in = x_in(3,k)/domain_top
    x_out(3) = (eta_in - eta(kk))/(eta(kk+1) - eta(kk))
    field_out(k) = 0.0_r_def
    do df = 1, ndf_f
      gamma(:) = field_proxy%vspace%evaluate_basis(df,x_out)
      dfk = map_f(df) + out_layer(k) 
      field_out(k) = field_out(k) + gamma(1)*field_proxy%data(dfk)
    end do
  end do

  deallocate ( chi_cell, dgamma, out_layer, eta )

  end subroutine evaluate_output_field


end module evaluate_output_field_mod


