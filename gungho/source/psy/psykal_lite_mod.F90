!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_mod

  use field_mod,             only : field_type, field_proxy_type 
  use scalar_mod,            only : scalar_type
  use operator_mod,          only : operator_type, operator_proxy_type
  use quadrature_mod,        only : quadrature_type
  use constants_mod,         only : r_def, i_def, cache_block
  use mesh_mod,              only : mesh_type
  use function_space_mod,    only : BASIS, DIFF_BASIS

  ! The following modules are not currently implemented as part of the main
  ! code but they do have unit tests so need to be declared here so they are
  ! built ready for unit testing.

  use quadrature_xoyoz_mod
  use quadrature_xyz_mod
  use quadrature_xyoz_mod, only : quadrature_xyoz_type, &
                                  quadrature_xyoz_proxy_type

  implicit none
  public

contains

    !------------------------------------------------------------------------------------------
    !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
    !> (not quadrature) will be removed and the functionality will be implemented via
    !> PSY using kernal meta data. 
    !> PSyclone support is required - documented in #942.
    SUBROUTINE invoke_compute_geopotential_kernel_type(geopotential, chi)
      USE compute_geopotential_kernel_mod, ONLY: compute_geopotential_code
      use mesh_mod, only: mesh_type
      TYPE(field_type), intent(inout)      :: geopotential, chi(3)

      INTEGER, pointer :: map_w0(:) => null()
      INTEGER, pointer :: map_chi(:) => null()
      REAL(kind=r_def), pointer :: nodes_w0(:,:) => null()
      INTEGER cell
      REAL(KIND=r_def), allocatable :: basis_chi(:,:,:)
      INTEGER ndf_w0, undf_w0, ndf_chi, undf_chi, dim_chi
      INTEGER :: df_w0, df_chi
      type(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) geopotential_proxy, chi_proxy(3)
      !
      ! Initialise field proxies
      !
      geopotential_proxy = geopotential%get_proxy()
      chi_proxy(1) = chi(1)%get_proxy()
      chi_proxy(2) = chi(2)%get_proxy()
      chi_proxy(3) = chi(3)%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = geopotential_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => geopotential%get_mesh()
      !
      ! Initialise sizes and allocate any basis arrays for w0
      !
      ndf_w0 = geopotential_proxy%vspace%get_ndf()
      undf_w0 = geopotential_proxy%vspace%get_undf()
      nodes_w0 => geopotential_proxy%vspace%get_nodes()

      !
      ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
      undf_chi  = chi_proxy(1)%vspace%get_undf( )
      dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

      ! Evaluate the basis function
      allocate( basis_chi(dim_chi, ndf_chi, ndf_w0) )
      do df_w0 = 1, ndf_w0
        do df_chi = 1, ndf_chi
          basis_chi(:,df_chi,df_w0) = chi_proxy(1)%vspace%call_function(BASIS,df_chi,nodes_w0(:,df_w0))
        end do
      end do

      DO cell=1, geopotential_proxy%vspace%get_ncell()

        map_w0 => geopotential_proxy%vspace%get_cell_dofmap(cell)
        map_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)

        CALL compute_geopotential_code(nlayers, geopotential_proxy%data, chi_proxy(1)%data, chi_proxy(2)%data, chi_proxy(3)%data, ndf_w0, undf_w0, map_w0, ndf_chi, undf_chi, map_chi, basis_chi)
      END DO
    END SUBROUTINE invoke_compute_geopotential_kernel_type

  !------------------------------------------------------------------------------
  !> invoke_initial_theta_kernel: invoke the potential temperature initialization for a generic space
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data.
  !> PSyclone support is required - documented in #942
  subroutine invoke_initial_theta_kernel( theta, chi )

    use initial_theta_kernel_mod, only : initial_theta_code
    use mesh_mod,                 only : mesh_type
    implicit none

    type( field_type ), intent( inout )  :: theta
    type( field_type ), intent( in )     :: chi(3)

    integer          :: cell
    integer          :: ndf_wtheta, undf_wtheta, &
                        ndf_chi, undf_chi, dim_chi, &
                        df_wtheta, df_chi
    integer, pointer :: map_wtheta(:) => null()
    integer, pointer :: map_chi(:)    => null()
    real(kind=r_def), pointer :: nodes_wtheta(:,:) => null()

    type( field_proxy_type ) :: theta_proxy
    type( field_proxy_type ) :: chi_proxy(3)

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    type(mesh_type), pointer :: mesh => null()

    theta_proxy  = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ndf_wtheta   = theta_proxy%vspace%get_ndf( )
    undf_wtheta  = theta_proxy%vspace%get_undf( )
    nodes_wtheta => theta_proxy%vspace%get_nodes()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi  = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate( basis_chi(dim_chi, ndf_chi, ndf_wtheta) )
    do df_wtheta = 1, ndf_wtheta
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_wtheta) = &
                  chi_proxy(1)%vspace%call_function(BASIS,df_chi,nodes_wtheta(:,df_wtheta))
      end do
    end do


    if (chi_proxy(1)%is_dirty(depth=1)) call chi_proxy(1)%halo_exchange(depth=1)
    if (chi_proxy(2)%is_dirty(depth=1)) call chi_proxy(2)%halo_exchange(depth=1)
    if (chi_proxy(3)%is_dirty(depth=1)) call chi_proxy(3)%halo_exchange(depth=1)

    mesh => theta%get_mesh()
    do cell = 1, mesh%get_last_halo_cell(1)

      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )
      map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      call initial_theta_code(       &
        theta_proxy%vspace%get_nlayers(),   &
        ndf_wtheta,                         &
        undf_wtheta,                        &
        map_wtheta,                         &
        theta_proxy%data,                   &
        ndf_chi,                            &
        undf_chi,                           &
        map_chi,                            &
        basis_chi,                          &
        chi_proxy(1)%data,                  &
        chi_proxy(2)%data,                  &
        chi_proxy(3)%data                   &
        )
    end do

    call theta_proxy%set_dirty()
  end subroutine invoke_initial_theta_kernel

  !------------------------------------------------------------------------------
  !> invoke_initial_rho_sample_kernel: invoke the density initialization for a generic space
  !> Computation of nodal basis function for coordinates chi not currently supported PSyClone,
  !> will be introduced in modification of quadrature strategy, see ticket #723.

  subroutine invoke_initial_rho_sample_kernel( rho, chi, time )

    use initial_rho_sample_kernel_mod, only : initial_rho_sample_code
    use mesh_mod,                      only : mesh_type
    implicit none

    type( field_type ), intent( inout )  :: rho
    type( field_type ), intent( in )     :: chi(3)
    real(kind=r_def),   intent( in )     :: time

    integer          :: cell
    integer          :: ndf_w3, undf_w3, &
                        ndf_chi, undf_chi, dim_chi
    integer, pointer :: map_w3(:) => null()
    integer, pointer :: map_chi(:)    => null()

    type( field_proxy_type ) :: rho_proxy
    type( field_proxy_type ) :: chi_proxy(3)

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    type(mesh_type), pointer :: mesh => null()

    integer          :: df_w3, df_chi
    real(kind=r_def), pointer :: nodes_w3(:,:) => null()

    rho_proxy    = rho%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ndf_w3       = rho_proxy%vspace%get_ndf( )
    undf_w3      = rho_proxy%vspace%get_undf( )

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi  = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

    allocate( basis_chi(dim_chi, ndf_chi, ndf_w3) )

    nodes_w3 => rho_proxy%vspace%get_nodes()
    do df_w3 = 1, ndf_w3
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_w3) = &
                  chi_proxy(1)%vspace%call_function(BASIS,df_chi,nodes_w3(:,df_w3))
      end do
    end do

    if (chi_proxy(1)%is_dirty(depth=1)) call chi_proxy(1)%halo_exchange(depth=1)
    if (chi_proxy(2)%is_dirty(depth=1)) call chi_proxy(2)%halo_exchange(depth=1)
    if (chi_proxy(3)%is_dirty(depth=1)) call chi_proxy(3)%halo_exchange(depth=1)

    mesh => rho%get_mesh()
    do cell = 1, mesh%get_last_halo_cell(1)

      map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
      map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      call initial_rho_sample_code(         &
        rho_proxy%vspace%get_nlayers(),     &
        ndf_w3,                             &
        undf_w3,                            &
        map_w3,                             &
        rho_proxy%data,                     &
        ndf_chi,                            &
        undf_chi,                           &
        map_chi,                            &
        basis_chi,                          &
        chi_proxy(1)%data,                  &
        chi_proxy(2)%data,                  &
        chi_proxy(3)%data,                  &
        time                                &
        )
    end do

    call rho_proxy%set_dirty()
  end subroutine invoke_initial_rho_sample_kernel


  !------------------------------------------------------------------------------
  !> invoke_initial_mr_kernel: invoke the moisture initialization
  subroutine invoke_initial_mr_kernel( theta, rho_in_wth, mr, chi )

    use initial_mr_kernel_mod, only : initial_mr_code
    use mr_indices_mod,        only : imr_v, imr_c, imr_r, imr_nc, imr_nr, nummr
    use mesh_mod,              only : mesh_type
    implicit none

    type( field_type ), intent( inout ) :: mr(nummr)
    type( field_type ), intent( in ) :: theta, rho_in_wth
    type( field_type ), intent( in ) :: chi(3)

    integer          :: cell
    integer          :: ndf_wtheta, undf_wtheta, &
                        ndf_chi, undf_chi
    integer, pointer :: map_wtheta(:) => null()
    integer, pointer :: map_chi(:)    => null()

    type( field_proxy_type ) :: mr_proxy(nummr)
    type( field_proxy_type ) :: theta_proxy, rho_proxy
    type( field_proxy_type ) :: chi_proxy(3)

    type(mesh_type), pointer :: mesh => null()

    integer :: imr

    do imr = 1,nummr
      mr_proxy(imr)  = mr(imr)%get_proxy()
    end do
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    theta_proxy = theta%get_proxy()
    rho_proxy = rho_in_wth%get_proxy()

    ndf_wtheta  = theta_proxy%vspace%get_ndf( )
    undf_wtheta = theta_proxy%vspace%get_undf( )

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf( )

    if (chi_proxy(1)%is_dirty(depth=1)) call chi_proxy(1)%halo_exchange(depth=1)
    if (chi_proxy(2)%is_dirty(depth=1)) call chi_proxy(2)%halo_exchange(depth=1)
    if (chi_proxy(3)%is_dirty(depth=1)) call chi_proxy(3)%halo_exchange(depth=1)
    do imr = 1,nummr
      if (mr_proxy(imr)%is_dirty(depth=1))  call mr_proxy(imr)%halo_exchange(depth=1)
    end do

    mesh => theta%get_mesh()
    do cell = 1, mesh%get_last_halo_cell(1)

      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )
      map_chi    => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      !!todo: nummr must be 5 for this to work, since we currently cant
      ! pass variable sized vectors into a kernel
      call initial_mr_code(       &
        theta_proxy%vspace%get_nlayers(),   &
        ndf_wtheta,                         &
        undf_wtheta,                        &
        map_wtheta,                         &
        theta_proxy%data,                   &
        rho_proxy%data,                     &
        mr_proxy(imr_v)%data,               &
        mr_proxy(imr_c)%data,               &
        mr_proxy(imr_r)%data,               &
        mr_proxy(imr_nc)%data,              &
        mr_proxy(imr_nr)%data,              &
        ndf_chi,                            &
        undf_chi,                           &
        map_chi,                            &
        chi_proxy(1)%data,                  &
        chi_proxy(2)%data,                  &
        chi_proxy(3)%data                   &
        )
    end do
    do imr = 1,nummr
       call mr_proxy(imr)%set_dirty()
    end do

  end subroutine invoke_initial_mr_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_rtheta_bd_kernel: Invoke the boundary part of the RHS of the theta equation
  subroutine invoke_rtheta_bd_kernel( r_theta_bd, theta, u, qr )

    use rtheta_bd_kernel_mod,   only : rtheta_bd_code
    use mesh_mod, only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,     only : stencil_dofmap_type, STENCIL_CROSS
    use reference_element_mod,  only : nfaces_h

    implicit none

    type (mesh_type), pointer            :: mesh => null()
    type( field_type ), intent( in )     :: theta, u
    type( field_type ), intent( inout )  :: r_theta_bd
    type( quadrature_type), intent( in ) :: qr

    type(stencil_dofmap_type), pointer :: cross_stencil_w2 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_wtheta
    integer                 :: undf_w2, undf_wtheta
    integer                 :: dim_w2, dim_wtheta
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer                 :: cross_stencil_w2_size

    integer, pointer        :: cross_stencil_wtheta_map(:,:,:) => null()
    integer                 :: cross_stencil_wtheta_size

    integer                 :: ff

    type( field_proxy_type )        :: r_theta_bd_proxy, u_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_wtheta_face(:,:,:,:,:),  &
                                      basis_w2_face(:,:,:,:,:)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: xp_f(:,:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    mesh => theta%get_mesh()

    r_theta_bd_proxy = r_theta_bd%get_proxy()
    theta_proxy      = theta%get_proxy()
    u_proxy          = u%get_proxy()

    cross_stencil_w2 => u_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w2_map => cross_stencil_w2%get_whole_dofmap()
    cross_stencil_w2_size = cross_stencil_w2%get_size()

    cross_stencil_wtheta => theta_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
    cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

    nlayers = theta_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    xp=>qr%get_xqp_h()
    zp=>qr%get_xqp_v()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal quad points in x and y
    nqp_h_1d = int(sqrt(real(nqp_h)))

    allocate(xp_f(nfaces_h, nqp_h_1d, 2))

    ndf_w2      = u_proxy%vspace%get_ndf( )
    dim_w2      = u_proxy%vspace%get_dim_space( )
    undf_w2     = u_proxy%vspace%get_undf()
    allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    ! Quadrature points on horizontal faces

    xp_f(1, :, :) = xp(1:nqp_h_1d, :)
    xp_f(1, :, 1) = 0.0_r_def

    xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f(2, :, 2) = 0.0_r_def

    xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f(3, :, 1) = 1.0_r_def

    xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f(4, :, 2) = 1.0_r_def

    ! Filling up the basis vector with value of the basis functions at the horizontal faces quadrature points

    do ff = 1, nfaces_h

      call u_proxy%vspace%compute_basis_function( &
        basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      call theta_proxy%vspace%compute_basis_function( &
        basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

    end do

    if(r_theta_bd_proxy%is_dirty(depth=2) ) call r_theta_bd_proxy%halo_exchange(depth=2)
    if(theta_proxy%is_dirty(depth=2) ) call theta_proxy%halo_exchange(depth=2)
    if(u_proxy%is_dirty(depth=2) ) call u_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()

    do cell = 1, mesh%get_last_halo_cell(1)

      call rtheta_bd_code(nlayers,                                         &
                          ndf_w2, undf_w2,                                 &
                          cross_stencil_w2_map(:,:,cell),                  &
                          cross_stencil_w2_size,                           &
                          ndf_wtheta, undf_wtheta,                         &
                          cross_stencil_wtheta_map(:,:,cell),              &
                          cross_stencil_wtheta_size,                       &
                          r_theta_bd_proxy%data,                           &
                          theta_proxy%data,                                &
                          u_proxy%data,                                    &
                          nqp_v, nqp_h_1d, wv,                             &
                          basis_w2_face, basis_wtheta_face,                &
                          adjacent_face(:,cell) )

    end do
    call r_theta_bd_proxy%set_dirty()

    deallocate(basis_w2_face, basis_wtheta_face)

  end subroutine invoke_rtheta_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_ru_bd_kernel: Invoke the boundary part of the RHS of the momentum equation
  subroutine invoke_ru_bd_kernel( r_u_bd, rho, theta, qr )

    use ru_bd_kernel_mod,       only : ru_bd_code
    use mesh_mod,               only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,     only : stencil_dofmap_type, STENCIL_CROSS
    use reference_element_mod,  only : nfaces_h

    implicit none

    type( mesh_type ), pointer           :: mesh => null()
    type( field_type ), intent( in )     :: rho, theta
    type( field_type ), intent( inout )  :: r_u_bd
    type( quadrature_type), intent( in ) :: qr

    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w2(:) => null()
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size

    integer, pointer        :: cross_stencil_wtheta_map(:,:,:) => null()
    integer                 :: cross_stencil_wtheta_size

    integer                 :: ff

    type( field_proxy_type )        :: r_u_bd_proxy, rho_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:), &
      basis_w3_face(:,:,:,:,:), &
      basis_wtheta_face(:,:,:,:,:)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: xp_f(:,:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:) => null()
    real(kind=r_def), pointer :: wv(:) => null()

    mesh => rho%get_mesh()

    r_u_bd_proxy = r_u_bd%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()

    cross_stencil_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    cross_stencil_wtheta => theta_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
    cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

    nlayers = rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal qp in x and y
    nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

    allocate(xp_f(nfaces_h, nqp_h_1d, 2))

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_w3  = rho_proxy%vspace%get_ndf( )
    dim_w3  = rho_proxy%vspace%get_dim_space( )
    undf_w3 = rho_proxy%vspace%get_undf()
    allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    ! Quadrature points on horizontal faces

    xp_f(1, :, :) = xp(1:nqp_h_1d, :)
    xp_f(1, :, 1) = 0.0_r_def

    xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f(2, :, 2) = 0.0_r_def

    xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f(3, :, 1) = 1.0_r_def

    xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f(4, :, 2) = 1.0_r_def

    ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

    do ff = 1, nfaces_h

      call r_u_bd_proxy%vspace%compute_basis_function( &
        basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      call rho_proxy%vspace%compute_basis_function( &
        basis_w3_face(ff,:,:,:,:), ndf_w3, nqp_h_1d, nqp_v, xp_f(ff, :,:), zp)

      call theta_proxy%vspace%compute_basis_function( &
        basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

    end do

    if(theta_proxy%is_dirty(depth=2) ) call theta_proxy%halo_exchange(depth=2)
    if(rho_proxy%is_dirty(depth=2) ) call rho_proxy%halo_exchange(depth=2)
    if(r_u_bd_proxy%is_dirty(depth=2) ) call r_u_bd_proxy%halo_exchange(depth=2)

    
    adjacent_face => mesh%get_adjacent_face()

    do cell = 1, mesh%get_last_halo_cell(1)

      map_w2 => r_u_bd_proxy%vspace%get_cell_dofmap( cell )

      call ru_bd_code(nlayers,                            &
                      ndf_w2, undf_w2,                    &
                      map_w2,                             &
                      ndf_w3, undf_w3,                    &
                      cross_stencil_w3_map(:,:,cell),     &
                      cross_stencil_w3_size,              &
                      ndf_wtheta, undf_wtheta,            &
                      cross_stencil_wtheta_map(:,:,cell), &
                      cross_stencil_wtheta_size,          &
                      r_u_bd_proxy%data,                  &
                      rho_proxy%data,                     &
                      theta_proxy%data,                   &
                      nqp_v, nqp_h_1d, wv,                &
                      basis_w2_face,                      &
                      basis_w3_face, basis_wtheta_face,   &
                      adjacent_face(:,cell))

    end do
    call r_u_bd_proxy%set_dirty()

    deallocate(basis_w3_face, basis_w2_face, basis_wtheta_face )

  end subroutine invoke_ru_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_exner_gradient_bd_kernel: Invoke the boundary integral for the exner_gradient kernel
  subroutine invoke_exner_gradient_bd_kernel( r_u_bd, exner, theta, qr )

    use exner_gradient_bd_kernel_mod, only : exner_gradient_bd_code
    use mesh_mod,                     only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,           only : stencil_dofmap_type, STENCIL_CROSS
    use reference_element_mod,        only : nfaces_h

    implicit none

    type( mesh_type ), pointer           :: mesh => null()
    type( field_type ), intent( in )     :: exner, theta
    type( field_type ), intent( inout )  :: r_u_bd
    type( quadrature_type), intent( in ) :: qr

    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w2(:) => null()
    integer, pointer        :: map_wtheta(:) => null()
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size


    integer                 :: ff

    type( field_proxy_type )        :: r_u_bd_proxy, exner_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:)
    real(kind=r_def), allocatable  :: basis_w3_face(:,:,:,:,:)
    real(kind=r_def), allocatable  :: basis_wtheta_face(:,:,:,:,:)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: xp_f(:,:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:) => null()
    real(kind=r_def), pointer :: wv(:) => null()

    mesh => exner%get_mesh()

    r_u_bd_proxy = r_u_bd%get_proxy()
    exner_proxy  = exner%get_proxy()
    theta_proxy  = theta%get_proxy()

    cross_stencil_w3 => exner_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    nlayers = exner_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal qp in x and y
    nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

    allocate(xp_f(nfaces_h, nqp_h_1d, 2))

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_w3  = exner_proxy%vspace%get_ndf( )
    dim_w3  = exner_proxy%vspace%get_dim_space( )
    undf_w3 = exner_proxy%vspace%get_undf()
    allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    ! Quadrature points on horizontal faces

    xp_f(1, :, :) = xp(1:nqp_h_1d, :)
    xp_f(1, :, 1) = 0.0_r_def

    xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f(2, :, 2) = 0.0_r_def

    xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f(3, :, 1) = 1.0_r_def

    xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f(4, :, 2) = 1.0_r_def

    ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

    do ff = 1, nfaces_h

      call r_u_bd_proxy%vspace%compute_basis_function( &
        basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      call exner_proxy%vspace%compute_basis_function( &
        basis_w3_face(ff,:,:,:,:), ndf_w3, nqp_h_1d, nqp_v, xp_f(ff, :,:), zp)

      call theta_proxy%vspace%compute_basis_function( &
        basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

    end do

    if( theta_proxy%is_dirty(depth=1) ) call theta_proxy%halo_exchange(depth=1)

    if(exner_proxy%is_dirty(depth=2) ) call exner_proxy%halo_exchange(depth=2)
    if(r_u_bd_proxy%is_dirty(depth=2) ) call r_u_bd_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()

    do cell = 1,mesh%get_last_halo_cell(1)

      map_w2 => r_u_bd_proxy%vspace%get_cell_dofmap( cell )
      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )

      call exner_gradient_bd_code( nlayers,                            &
                                   ndf_w2, undf_w2,                    &
                                   map_w2,                             &
                                   ndf_w3, undf_w3,                    &
                                   cross_stencil_w3_map(:,:,cell),     &
                                   cross_stencil_w3_size,              &
                                   ndf_wtheta, undf_wtheta,            &
                                   map_wtheta,                         &
                                   r_u_bd_proxy%data,                  &
                                   exner_proxy%data,                   &
                                   theta_proxy%data,                   &
                                   nqp_v, nqp_h_1d, wv,                &
                                   basis_w2_face,                      &
                                   basis_w3_face, basis_wtheta_face,   &
                                   adjacent_face(:,cell))

    end do
    call r_u_bd_proxy%set_dirty()

    deallocate(basis_w3_face, basis_w2_face, basis_wtheta_face, xp_f)

  end subroutine invoke_exner_gradient_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_pert_pressure_gradient_bd_kernel: Invoke the boundary part of pert_pressure_gradient kernel.
  subroutine invoke_pert_pressure_gradient_bd_kernel( r_u_bd, rho, rho_ref, theta, theta_ref, qr )

    use pert_pressure_gradient_bd_kernel_mod, only : pert_pressure_gradient_bd_code
    use mesh_mod,                             only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,                   only : stencil_dofmap_type, STENCIL_CROSS
    use reference_element_mod,                only : nfaces_h

    implicit none

    type( mesh_type ), pointer           :: mesh => null()
    type( field_type ), intent( in )     :: rho, theta, rho_ref, theta_ref
    type( field_type ), intent( inout )  :: r_u_bd
    type( quadrature_type), intent( in ) :: qr

    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w2(:) => null()
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size

    integer, pointer        :: cross_stencil_wtheta_map(:,:,:) => null()
    integer                 :: cross_stencil_wtheta_size

    integer                 :: ff

    type( field_proxy_type )        :: r_u_bd_proxy, rho_proxy, theta_proxy, rho_ref_proxy, theta_ref_proxy

    real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:),    &
                                      basis_w3_face(:,:,:,:,:),    &
                                      basis_wtheta_face(:,:,:,:,:)

    real(kind=r_def), pointer :: xp_f(:,:,:) => null()
    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:) => null()
    real(kind=r_def), pointer :: wv(:) => null()

    mesh => rho%get_mesh()

    r_u_bd_proxy = r_u_bd%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()
    rho_ref_proxy    = rho_ref%get_proxy()
    theta_ref_proxy  = theta_ref%get_proxy()

    cross_stencil_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    cross_stencil_wtheta => theta_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
    cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

    nlayers = rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal qp in x and y
    nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

    allocate(xp_f(nfaces_h, nqp_h_1d, 2))

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_w3  = rho_proxy%vspace%get_ndf( )
    dim_w3  = rho_proxy%vspace%get_dim_space( )
    undf_w3 = rho_proxy%vspace%get_undf()
    allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    ! Quadrature points on horizontal faces

    xp_f(1, :, :) = xp(1:nqp_h_1d, :)
    xp_f(1, :, 1) = 0.0_r_def

    xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f(2, :, 2) = 0.0_r_def

    xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f(3, :, 1) = 1.0_r_def

    xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f(4, :, 2) = 1.0_r_def

    ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

    do ff = 1, nfaces_h

      call r_u_bd_proxy%vspace%compute_basis_function( &
        basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      call rho_proxy%vspace%compute_basis_function( &
        basis_w3_face(ff,:,:,:,:), ndf_w3, nqp_h_1d, nqp_v, xp_f(ff, :,:), zp)

      call theta_proxy%vspace%compute_basis_function( &
        basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

    end do

    if(theta_proxy%is_dirty(depth=2) ) call theta_proxy%halo_exchange(depth=2)
    if(rho_proxy%is_dirty(depth=2) ) call rho_proxy%halo_exchange(depth=2)
    if(theta_ref_proxy%is_dirty(depth=2) ) call theta_ref_proxy%halo_exchange(depth=2)
    if(rho_ref_proxy%is_dirty(depth=2) ) call rho_ref_proxy%halo_exchange(depth=2)
    if(r_u_bd_proxy%is_dirty(depth=2) ) call r_u_bd_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()

    do cell = 1, mesh%get_last_halo_cell(1)

      map_w2 => r_u_bd_proxy%vspace%get_cell_dofmap( cell )
      call pert_pressure_gradient_bd_code( nlayers,                            &
                                           ndf_w2, undf_w2,                    &
                                           map_w2,                             &
                                           ndf_w3, undf_w3,                    &
                                           cross_stencil_w3_map(:,:,cell),     &
                                           cross_stencil_w3_size,              &
                                           ndf_wtheta, undf_wtheta,            &
                                           cross_stencil_wtheta_map(:,:,cell), &
                                           cross_stencil_wtheta_size,          &
                                           r_u_bd_proxy%data,                  &
                                           rho_proxy%data,                     &
                                           rho_ref_proxy%data,                 &
                                           theta_proxy%data,                   &
                                           theta_ref_proxy%data,               &
                                           nqp_v, nqp_h_1d, wv,                &
                                           basis_w2_face,                      &
                                           basis_w3_face, basis_wtheta_face,   &
                                           adjacent_face(:,cell))

    end do
    call r_u_bd_proxy%set_dirty()

    deallocate(basis_w3_face, basis_w2_face, basis_wtheta_face, xp_f)

  end subroutine invoke_pert_pressure_gradient_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_div_bd_kernel: Invoke the boundary part of the divergence of the lhs Helmholtz
  subroutine invoke_weighted_div_bd_kernel_type(div_star, theta, qr)

      use weighted_div_bd_kernel_mod, only : weighted_div_bd_code
      use mesh_mod,                   only : mesh_type
      use reference_element_mod,      only : nfaces_h
      use stencil_dofmap_mod,         only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type( mesh_type ), pointer           :: mesh => null()
      type(field_type), intent(in)         :: theta
      type(operator_type), intent(inout)   :: div_star
      type(quadrature_type), intent(in)    :: qr

      integer :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
      integer :: ndf_w2, ndf_w3, ndf_wtheta, undf_wtheta
      integer :: dim_w2, dim_w3, dim_wtheta

      integer,                   pointer :: cross_stencil_wtheta_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()
      integer                            :: cross_stencil_wtheta_size
      integer                 :: ff

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:), &
      basis_w3_face(:,:,:,:,:), &
      basis_wtheta_face(:,:,:,:,:)

      real(kind=r_def), pointer :: xp(:,:) => null()
      real(kind=r_def), pointer :: xp_f(:,:,:) => null()
      real(kind=r_def), pointer :: zp(:) => null()
      real(kind=r_def), pointer :: wh(:) => null()
      real(kind=r_def), pointer :: wv(:) => null()

      type(operator_proxy_type) ::div_star_proxy
      type(field_proxy_type)    ::theta_proxy

      integer, pointer        :: adjacent_face(:,:) => null()

      !
      ! Initialise field proxies
      !
      div_star_proxy = div_star%get_proxy()
      theta_proxy = theta%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = div_star_proxy%fs_from%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => div_star%get_mesh()

      !
      ! Initialise qr values
      !
      nqp_h=qr%get_nqp_h()
      nqp_v=qr%get_nqp_v()
      zp=>qr%get_xqp_v()
      xp=>qr%get_xqp_h()
      wh=>qr%get_wqp_h()
      wv=>qr%get_wqp_v()

      ! Assumes same number of horizontal qp in x and y
      nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

      allocate(xp_f(nfaces_h, nqp_h_1d, 2))

      ndf_w2      = div_star_proxy%fs_to%get_ndf( )
      dim_w2      = div_star_proxy%fs_to%get_dim_space( )
      allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

      ndf_w3  = div_star_proxy%fs_from%get_ndf( )
      dim_w3  = div_star_proxy%fs_from%get_dim_space( )
      allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

      ndf_wtheta      = theta_proxy%vspace%get_ndf( )
      dim_wtheta      = theta_proxy%vspace%get_dim_space( )
      undf_wtheta     = theta_proxy%vspace%get_undf()
      allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

      ! Quadrature points on horizontal faces

      xp_f(1, :, :) = xp(1:nqp_h_1d, :)
      xp_f(1, :, 1) = 0.0_r_def

      xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
      xp_f(2, :, 2) = 0.0_r_def

      xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
      xp_f(3, :, 1) = 1.0_r_def

      xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
      xp_f(4, :, 2) = 1.0_r_def

      ! Stencil maps
      cross_stencil_wtheta => theta_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
      cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

      ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

      do ff = 1, nfaces_h

        call div_star_proxy%fs_to%compute_basis_function( &
          basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

        call div_star_proxy%fs_from%compute_basis_function( &
          basis_w3_face(ff,:,:,:,:), ndf_w3, nqp_h_1d, nqp_v, xp_f(ff, :,:), zp)

        call theta_proxy%vspace%compute_basis_function( &
          basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      end do
      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=2)) call theta_proxy%halo_exchange(depth=2)      
      !

      adjacent_face => mesh%get_adjacent_face()
      do cell=1,mesh%get_last_halo_cell(1)        

        call weighted_div_bd_code(cell,                               &
                                  nlayers,                            &
                                  div_star_proxy%ncell_3d,            &
                                  div_star_proxy%local_stencil,       &
                                  theta_proxy%data,                   &
                                  ndf_w2, ndf_w3,                     &
                                  ndf_wtheta, undf_wtheta,            &
                                  cross_stencil_wtheta_map(:,:,cell), &
                                  cross_stencil_wtheta_size,          &
                                  nqp_v, nqp_h_1d, wv,                &
                                  basis_w2_face,                      &
                                  basis_w3_face, basis_wtheta_face,   &
                                  adjacent_face(:,cell))
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate (basis_w3_face, basis_w2_face, basis_wtheta_face)
      !
    end subroutine invoke_weighted_div_bd_kernel_type

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_proj_2theta_bd_kernel: Invoke the boundary part of the projection from wtheta to w2 operator of the lhs Helmholtz
  subroutine invoke_weighted_proj_2theta_bd_kernel_type(p2theta, theta, rho, qr)

      use weighted_proj_2theta_bd_kernel_mod, only : weighted_proj_2theta_bd_code
      use mesh_mod,                           only : mesh_type
      use reference_element_mod,              only : nfaces_h
      use stencil_dofmap_mod,                 only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type( mesh_type ), pointer           :: mesh => null()
      type(field_type), intent(in)         :: theta, rho
      type(operator_type), intent(inout)   :: p2theta
      type(quadrature_type), intent(in)    :: qr

      integer :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
      integer :: ndf_w2, ndf_w3, undf_w3, ndf_wtheta, undf_wtheta
      integer :: dim_w2, dim_w3, dim_wtheta

      integer,                   pointer :: cross_stencil_wtheta_map(:,:,:) => null(), &
                                            cross_stencil_w3_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null(), &
                                            cross_stencil_w3 => null()

      integer                            :: cross_stencil_w3_size, cross_stencil_wtheta_size


      integer                 :: ff

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:), basis_w3_face(:,:,:,:,:), &
      basis_wtheta_face(:,:,:,:,:)

      real(kind=r_def), pointer :: xp(:,:), xp_f(:,:,:) => null()
      real(kind=r_def), pointer :: zp(:) => null()
      real(kind=r_def), pointer :: wv(:) => null()

      type(operator_proxy_type) ::p2theta_proxy
      type(field_proxy_type)    ::theta_proxy, rho_proxy

      integer, pointer        :: adjacent_face(:,:) => null()

      !
      ! Initialise field proxies
      !
      p2theta_proxy = p2theta%get_proxy()
      theta_proxy = theta%get_proxy()
      rho_proxy   = rho%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = p2theta_proxy%fs_from%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => p2theta%get_mesh()

      !
      ! Initialise qr values
      !
      nqp_h=qr%get_nqp_h()
      nqp_v=qr%get_nqp_v()
      zp=>qr%get_xqp_v()
      xp=>qr%get_xqp_h()
      wv=>qr%get_wqp_v()

      ! Assumes same number of horizontal qp in x and y
      nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

      allocate(xp_f(nfaces_h, nqp_h_1d, 2))

      ndf_w2      = p2theta_proxy%fs_to%get_ndf( )
      dim_w2      = p2theta_proxy%fs_to%get_dim_space( )
      allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

      ndf_w3   = rho_proxy%vspace%get_ndf( )
      undf_w3  = rho_proxy%vspace%get_undf( )
      dim_w3   = rho_proxy%vspace%get_dim_space( )
      allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

      ndf_wtheta      = p2theta_proxy%fs_from%get_ndf( )
      dim_wtheta      = p2theta_proxy%fs_from%get_dim_space( )
      undf_wtheta     = p2theta_proxy%fs_from%get_undf()
      allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

      ! Quadrature points on horizontal faces

      xp_f(1, :, :) = xp(1:nqp_h_1d, :)
      xp_f(1, :, 1) = 0.0_r_def

      xp_f(2, :, :) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
      xp_f(2, :, 2) = 0.0_r_def

      xp_f(3, :, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
      xp_f(3, :, 1) = 1.0_r_def

      xp_f(4, :, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
      xp_f(4, :, 2) = 1.0_r_def

      ! Stencil maps
      cross_stencil_wtheta => p2theta_proxy%fs_from%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
      cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

      cross_stencil_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
      cross_stencil_w3_size = cross_stencil_w3%get_size()

      ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

      do ff = 1, nfaces_h

        call p2theta_proxy%fs_to%compute_basis_function( &
          basis_w2_face(ff,:,:,:,:), ndf_w2, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

        call p2theta_proxy%fs_from%compute_basis_function( &
          basis_wtheta_face(ff,:,:,:,:), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(ff, :,:), zp)

        call rho_proxy%vspace%compute_basis_function( &
          basis_w3_face(ff,:,:,:,:), ndf_w3, nqp_h_1d, nqp_v, xp_f(ff, :, :), zp)

      end do
      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=2)) call theta_proxy%halo_exchange(depth=2)
      if (rho_proxy%is_dirty(depth=2))   call rho_proxy%halo_exchange(depth=2)
      !

      adjacent_face => mesh%get_adjacent_face()
      do cell=1,mesh%get_last_halo_cell(1)

        call weighted_proj_2theta_bd_code(cell,                       &
                                  nlayers,                            &
                                  p2theta_proxy%ncell_3d,             &
                                  p2theta_proxy%local_stencil,        &
                                  theta_proxy%data,                   &
                                  rho_proxy%data,                     &
                                  ndf_w2,                             &
                                  ndf_wtheta, undf_wtheta,            &
                                  cross_stencil_wtheta_map(:,:,cell), &
                                  cross_stencil_wtheta_size,          &
                                  ndf_w3, undf_w3,                    &
                                  cross_stencil_w3_map(:,:,cell),     &
                                  cross_stencil_w3_size,              &
                                  nqp_h_1d, nqp_v, wv,                &
                                  basis_w2_face,                      &
                                  basis_w3_face, basis_wtheta_face,   &
                                  adjacent_face(:,cell))
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate (basis_w3_face, basis_wtheta_face, basis_w2_face)
      !
    end subroutine invoke_weighted_proj_2theta_bd_kernel_type


  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_proj_theta2_bd_kernel: Invoke the boundary part of the projection from w2 to wtheta operator of the lhs Helmholtz  
  subroutine invoke_weighted_proj_theta2_bd_kernel_type(ptheta2, theta, qr)

      use weighted_proj_theta2_bd_kernel_mod, only : weighted_proj_theta2_bd_code
      use mesh_mod,                           only : mesh_type
      use reference_element_mod,              only : nfaces_h
      use stencil_dofmap_mod,                 only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type( mesh_type ), pointer           :: mesh => null()
      type(field_type), intent(in)         :: theta
      type(operator_type), intent(inout)   :: ptheta2
      type(quadrature_type), intent(in)    :: qr

      integer :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
      integer :: ndf_w2, ndf_wtheta, undf_wtheta
      integer :: dim_w2, dim_wtheta

      integer,                   pointer :: cross_stencil_wtheta_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()

      integer                            :: cross_stencil_wtheta_size

      integer                 :: ff

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:,:), basis_wtheta_face(:,:,:,:,:)

      real(kind=r_def), pointer :: xp(:,:), xp_f(:,:,:) => null()
      real(kind=r_def), pointer :: zp(:) => null()
      real(kind=r_def), pointer :: wv(:) => null()

      type(operator_proxy_type) ::ptheta2_proxy
      type(field_proxy_type)    ::theta_proxy

      integer, pointer        :: adjacent_face(:,:) => null()

      !
      ! Initialise field proxies
      !
      ptheta2_proxy = ptheta2%get_proxy()
      theta_proxy = theta%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = ptheta2_proxy%fs_from%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => ptheta2%get_mesh()

      !
      ! Initialise qr values
      !
      nqp_h=qr%get_nqp_h()
      nqp_v=qr%get_nqp_v()
      zp=>qr%get_xqp_v()
      xp=>qr%get_xqp_h()
      wv=>qr%get_wqp_v()

      ! Assumes same number of horizontal qp in x and y
      nqp_h_1d = int(sqrt(real(nqp_h)))  ! use sqrt

      allocate(xp_f(nqp_h_1d, 2, nfaces_h))

      ndf_w2      = ptheta2_proxy%fs_from%get_ndf( )
      dim_w2      = ptheta2_proxy%fs_from%get_dim_space( )
      allocate(basis_w2_face(dim_w2,ndf_w2,nqp_h_1d,nqp_v,nfaces_h))

      ndf_wtheta      = ptheta2_proxy%fs_to%get_ndf( )
      dim_wtheta      = ptheta2_proxy%fs_to%get_dim_space( )
      undf_wtheta     = ptheta2_proxy%fs_to%get_undf()
      allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v,nfaces_h))

      ! Quadrature points on horizontal faces

      xp_f(:, :, 1) = xp(1:nqp_h_1d, :)
      xp_f(:, 1, 1) = 0.0_r_def

      xp_f(:, :, 2) = xp(1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
      xp_f(:, 2, 2) = 0.0_r_def

      xp_f(:, :, 3) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
      xp_f(:, 1, 3) = 1.0_r_def

      xp_f(:, :, 4) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
      xp_f(:, 2, 4) = 1.0_r_def

      ! Stencil maps
      cross_stencil_wtheta => ptheta2_proxy%fs_to%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
      cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

      ! Filling up the face basis vector with value of the basis functions at the horizontal faces quadrature points

      do ff = 1, nfaces_h

        call ptheta2_proxy%fs_from%compute_basis_function( &
          basis_w2_face(:,:,:,:,ff), ndf_w2, nqp_h_1d, nqp_v, xp_f(:, :, ff), zp)

        call ptheta2_proxy%fs_to%compute_basis_function( &
          basis_wtheta_face(:,:,:,:,ff), ndf_wtheta, nqp_h_1d, nqp_v, xp_f(:, :, ff), zp)

      end do
      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=2)) call theta_proxy%halo_exchange(depth=2)
      !
      adjacent_face => mesh%get_adjacent_face()
      do cell=1,mesh%get_last_halo_cell(1)

        call weighted_proj_theta2_bd_code(cell,                       &
                                  nlayers,                            &
                                  ptheta2_proxy%ncell_3d,             &
                                  ptheta2_proxy%local_stencil,        &
                                  theta_proxy%data,                   &
                                  ndf_w2,                             &
                                  ndf_wtheta, undf_wtheta,            &
                                  cross_stencil_wtheta_map(:,:,cell), &
                                  cross_stencil_wtheta_size,          &
                                  nqp_h_1d, nqp_v, wv,                &
                                  basis_w2_face,                      &
                                  basis_wtheta_face,                  &
                                  adjacent_face(:,cell))
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate (basis_wtheta_face, basis_w2_face)
      !
    end subroutine invoke_weighted_proj_theta2_bd_kernel_type


!-------------------------------------------------------------------------------    
  subroutine invoke_inner_prod(x,y,inner_prod)
    use omp_lib
    USE log_mod, ONLY : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x,y
    real(kind=r_def),    intent(out) :: inner_prod
    real(kind=r_def)                 :: inner_prod_local_tmp
    real(kind=r_def), allocatable, dimension(:) :: lsum
    type( scalar_type )              :: inner_prod_local
    type( field_proxy_type)          ::  x_p,y_p
    INTEGER                          :: i,undf
    integer                          :: thread_id, num_threads, pad_size

    x_p = x%get_proxy()
    y_p = y%get_proxy()

    undf = x_p%vspace%get_last_dof_owned()
    !sanity check
    if(undf /= y_p%vspace%get_last_dof_owned() ) then
      ! they are not on the same function space
      call log_event("Psy:inner_prod:x and y live on different w-spaces",LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    ! get the number of threads
    num_threads = omp_get_max_threads()
    pad_size = cache_block / storage_size(x_p%data(1))
    ! allocate the array, pad to avoid false sharing
    allocate( lsum( 0 : (num_threads - 1)*pad_size ) )
    lsum = 0.0_r_def

    !$omp parallel default(none) &
    !$omp& shared(lsum, x_p, y_p, undf, num_threads, pad_size) & 
    !$omp& private(i, thread_id)
    ! which thread am i?
    thread_id = omp_get_thread_num()*pad_size

    !$omp do schedule(static) 
    do i = 1,undf
      lsum(thread_id) = lsum(thread_id) + &
                             ( x_p%data(i) * y_p%data(i) )
    end do
    !$omp end do

    !$omp end parallel

    inner_prod_local_tmp = 0.0_r_def 
    do i = 0, num_threads -1
       inner_prod_local_tmp = inner_prod_local_tmp + lsum(i*pad_size)
    end do
    deallocate(lsum)
    inner_prod_local = scalar_type( inner_prod_local_tmp )

    ! Get the global sum of the inner products
    inner_prod = inner_prod_local%get_sum()
    
  end subroutine invoke_inner_prod

!-------------------------------------------------------------------------------   
!> invoke_X_innerproduct_X:  Calculates inner product of a vector by itself  
  subroutine invoke_X_innerproduct_X(inner_prod, x)
    use omp_lib
    USE log_mod, ONLY : log_event, LOG_LEVEL_ERROR
    implicit none
    real(kind=r_def),    intent(out) :: inner_prod
    type( field_type ),  intent(in ) :: x
    real(kind=r_def)                 :: inner_prod_local_tmp
    real(kind=r_def), allocatable, dimension(:) :: lsum
    type( scalar_type )              :: inner_prod_local
    type( field_proxy_type)          :: x_p
    INTEGER                          :: i,undf
    integer                          :: thread_id, num_threads, pad_size

    x_p = x%get_proxy()

    undf = x_p%vspace%get_last_dof_owned()

    ! get the number of threads
    num_threads = omp_get_max_threads()
    pad_size = cache_block / storage_size(x_p%data(1))
    ! allocate the array, pad to avoid false sharing
    allocate( lsum( 0 : (num_threads - 1)*pad_size ) )
    lsum = 0.0_r_def

    !$omp parallel default(none) &
    !$omp& shared(lsum, x_p, undf, num_threads, pad_size) & 
    !$omp& private(i, thread_id)
    ! which thread am i?
    thread_id = omp_get_thread_num()*pad_size

    !$omp do schedule(static) 
    do i = 1,undf
      lsum(thread_id) = lsum(thread_id) + &
                             ( x_p%data(i) * x_p%data(i) )
    end do
    !$omp end do

    !$omp end parallel

    inner_prod_local_tmp = 0.0_r_def 
    do i = 0, num_threads -1
       inner_prod_local_tmp = inner_prod_local_tmp + lsum(i*pad_size)
    end do
    deallocate(lsum)
    inner_prod_local = scalar_type( inner_prod_local_tmp )

    ! Get the global sum of the inner products
    inner_prod = inner_prod_local%get_sum()
    
  end subroutine invoke_X_innerproduct_X
  
!-------------------------------------------------------------------------------   
!> invoke_axpy:  (a * x + y) ; a-scalar, x,y-vector     
  subroutine invoke_axpy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray

    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none) &
    !$omp&  shared(field1_proxy,field2_proxy, field_res_proxy, &
    !$omp&  undf, scalar),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) + field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_axpy
  
!-------------------------------------------------------------------------------   
!> invoke_axmy:  (a * x - y) ; a-scalar, x,y-vector
  subroutine invoke_axmy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf, scalar),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) - field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_axmy
  
!-------------------------------------------------------------------------------   
!> invoke_copy_field_data: Copy the data from one field to another ( a = b )
  subroutine invoke_copy_field_data(field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:copy_field_data:field1 and field_res live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()
    
    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_copy_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_minus_field_data: Subtract values of field2 from values of field 
!> ( c = a - b )
  subroutine invoke_minus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) - field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_minus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_plus_field_data:  Add values of field2 to values of field1
!> ( c = a + b )
  subroutine invoke_plus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) + field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_plus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_set_field_scalar: Set all values in a field to a single value
  subroutine invoke_set_field_scalar(scalar, field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field_res_proxy = field_res%get_proxy()

    undf = field_res_proxy%vspace%get_undf()

    !$omp parallel do schedule(static), default(none), shared(field_res_proxy, undf, scalar),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = scalar
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      call field_res_proxy%set_clean(dplp)
    end do

  end subroutine invoke_set_field_scalar

!-------------------------------------------------------------------------------
!> invoke_divide_field: Divide the values of field1 by field2 and put result in
!>field_res
!> c = a/b
  subroutine invoke_divide_field(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:divide_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:divide_field:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)/field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()
    
    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_divide_field

!-------------------------------------------------------------------------------   
!> invoke_copy_scaled_field_data: Copy the scaled data from one field to another ( a = scalar*b )
  subroutine invoke_copy_scaled_field_data(scalar,field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    real( kind=r_def ), intent(in)     :: scalar
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:copy_scaled_field_data:field1 and field_res live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field_res_proxy, undf, scalar),  private(i)
    do i = 1,undf
       field_res_proxy%data(i) = scalar*field1_proxy%data(i)
    end do
    !$omp end parallel do
   
    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()
        
    do dplp = 1, depth
       if( field1_proxy%is_dirty(depth=dplp) ) then
          call field_res_proxy%set_dirty()
       else
          call field_res_proxy%set_clean(dplp)
       end if
    end do
   
  end subroutine invoke_copy_scaled_field_data


!-------------------------------------------------------------------------------   
!> invoke_sum_field: Sum all values of a field x
  subroutine invoke_sum_field( x, field_sum )
    use omp_lib
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x
    real(kind=r_def),    intent(out) :: field_sum
    real(kind=r_def)                 :: field_sum_local_tmp
    real(kind=r_def), allocatable, dimension(:) :: lsum
    type( scalar_type )              :: field_sum_local
    type( field_proxy_type)          :: x_p
    integer(kind=i_def)              :: df, undf
    integer(kind=i_def)              :: thread_id, num_threads, pad_size

    x_p = x%get_proxy()   

    undf = x_p%vspace%get_last_dof_owned()
    
    ! Calculate the local sum on this partition
    ! get the number of threads
    num_threads = omp_get_max_threads()
    pad_size = cache_block / storage_size(x_p%data(1))
    ! allocate the array, pad to avoid false sharing
    allocate( lsum( 0 : (num_threads - 1)*pad_size ) )
    lsum = 0.0_r_def

    !$omp parallel default(none) &
    !$omp& shared(x_p, lsum, undf, num_threads, pad_size) &
    !$omp& private(thread_id, df)
    thread_id = omp_get_thread_num()*pad_size

    !$omp do schedule(static) 
    do df = 1,undf
      lsum(thread_id) = lsum(thread_id) + x_p%data(df)
    end do
    !$omp end do
    
    !$omp end parallel

    field_sum_local_tmp = 0.0_r_def
    do df = 0, num_threads - 1
       field_sum_local_tmp = field_sum_local_tmp + lsum(df*pad_size)
    end do
    deallocate(lsum)
    field_sum_local = scalar_type( field_sum_local_tmp )

    ! Now call the global sum
    field_sum = field_sum_local%get_sum()

  end subroutine invoke_sum_field

!-------------------------------------------------------------------------------   
!> invoke_axpby:  z = (a * x + b * y) ; a,b-scalar, x,y-vector     
  subroutine invoke_axpby(a,x,b,y,z)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: x, y
    type( field_type ), intent(inout ) :: z
    real(kind=r_def),   intent(in )    :: a, b
    type( field_proxy_type)            :: x_proxy,y_proxy      &
                                        , z_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    x_proxy = x%get_proxy()
    y_proxy = y%get_proxy()
    z_proxy = z%get_proxy()

    !sanity check
    undf = x_proxy%vspace%get_undf()
    if(undf /= y_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpby:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= z_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpby:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(z_proxy,x_proxy, y_proxy, undf, a, b),  private(i)
    do i = 1,undf
      z_proxy%data(i) = (a * x_proxy%data(i)) + (b * y_proxy%data(i))
    end do
    !$omp end parallel do

    mesh => z%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( x_proxy%is_dirty(depth=dplp) .or. y_proxy%is_dirty(depth=dplp) ) then
        call z_proxy%set_dirty()
      else
        call z_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_axpby

!-------------------------------------------------------------------------------   
!> Non pointwise Kernels


!------------------------------------------------------------------------------- 
!> invoke_sample_flux_kernel: Retrieve values from flux kernel 
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
  subroutine invoke_sample_flux_kernel(flux, u, multiplicity, q )
    use sample_flux_kernel_mod, only: sample_flux_code
    use mesh_mod,               only: mesh_type

    implicit none

    type(field_type), intent(in)         :: u, multiplicity, q
    type(field_type), intent(inout)      :: flux

    type(field_proxy_type)          :: u_p, m_p, q_p, flux_p
    
    integer                 :: cell, nlayers
    integer                 :: ndf_f, ndf_q
    integer                 :: undf_f, undf_q
    integer                 :: dim_q
    integer                 :: df_f, df_q

    integer, pointer        :: map_f(:) => null()
    integer, pointer        :: map_q(:) => null()
    real(kind=r_def), pointer :: nodes_f(:,:) => null()

    real(kind=r_def), allocatable  :: basis_q(:,:,:)

    type(mesh_type), pointer :: mesh => null()

    u_p    = u%get_proxy()
    q_p    = q%get_proxy()
    m_p    = multiplicity%get_proxy()
    flux_p = flux%get_proxy()

    nlayers = flux_p%vspace%get_nlayers()

    ndf_f  = flux_p%vspace%get_ndf( )
    undf_f = flux_p%vspace%get_undf()
    nodes_f => flux_p%vspace%get_nodes()

    ndf_q  = q_p%vspace%get_ndf( )
    dim_q  = q_p%vspace%get_dim_space( )
    undf_q = q_p%vspace%get_undf()

    ! Evaluate the basis function
    allocate(basis_q(dim_q, ndf_q, ndf_f))
    do df_f = 1, ndf_f
      do df_q = 1, ndf_q
        basis_q(:,df_q,df_f) = q_p%vspace%call_function(BASIS,df_q,nodes_f(:,df_f))
      end do
    end do

    mesh => flux%get_mesh()
    if (flux_p%is_dirty(depth=1)) call flux_p%halo_exchange(depth=1)
    if (   m_p%is_dirty(depth=1)) call    m_p%halo_exchange(depth=1)
    if (   u_p%is_dirty(depth=1)) call    u_p%halo_exchange(depth=1)
    if (   q_p%is_dirty(depth=1)) call    q_p%halo_exchange(depth=1)


    do cell = 1, mesh%get_last_halo_cell(1)
       map_f => flux_p%vspace%get_cell_dofmap( cell )
       map_q => q_p%vspace%get_cell_dofmap( cell )
       call sample_flux_code(nlayers, &
                             flux_p%data, & 
                             m_p%data, &
                             u_p%data, &
                             q_p%data, &
                             ndf_f, & 
                             undf_f, &
                             map_f, &
                             ndf_q, &
                             undf_q, &
                             map_q, &
                             basis_q &
                            )
    end do

    call flux_p%set_dirty()

  end subroutine invoke_sample_flux_kernel

  !-------------------------------------------------------------------------------
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data. 
  !> PSyclone support is required - documented in #942
  subroutine invoke_nodal_coordinates_kernel(nodal_coords, chi )
    use nodal_coordinates_kernel_mod, only: nodal_coordinates_code
    use mesh_mod,                     only: mesh_type
    implicit none
    
    type(field_type), intent(inout)      :: nodal_coords(3)
    type(field_type), intent(in)         :: chi(3) 

    type(field_proxy_type) :: x_p(3), chi_p(3)
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x
    integer                 :: undf_chi, undf_x
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf( )
    undf_x = x_p(1)%vspace%get_undf()
    nodes_x => x_p(1)%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    mesh => nodal_coords(1)%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call nodal_coordinates_code(nlayers, &
                                   x_p(1)%data, &
                                   x_p(2)%data, &
                                   x_p(3)%data, &
                                   chi_p(1)%data, &
                                   chi_p(2)%data, &
                                   chi_p(3)%data, &
                                   ndf_x, undf_x, map_x, &
                                   ndf_chi, undf_chi, map_chi, &
                                   basis_chi &
                                  )
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_nodal_coordinates_kernel

  !-------------------------------------------------------------------------------   
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data.
  !> PSyclone support is required - documented in #942
  subroutine invoke_convert_hcurl_field(phys_field, comp_field, chi )
    use convert_hcurl_field_kernel_mod, only: convert_hcurl_field_code
    use mesh_mod,                       only: mesh_type
    implicit none
    
    type(field_type), intent(inout)      :: phys_field(3)
    type(field_type), intent(in)         :: chi(3), comp_field
 
    type(field_proxy_type) :: phys_p(3), chi_p(3), comp_p
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_comp, ndf_phys
    integer                 :: undf_chi, undf_comp
    integer                 :: diff_dim_chi, dim_comp
    integer                 :: df_comp, df_chi, df_phys

    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map(:) => null()
    real(kind=r_def), pointer :: nodes_phys(:,:) => null()

    real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:), basis_comp(:,:,:)
    integer(kind=i_def) :: i
    type(mesh_type), pointer :: mesh => null()

    do i = 1,3
      phys_p(i) = phys_field(i)%get_proxy()
      chi_p(i)  = chi(i)%get_proxy()
    end do
    comp_p = comp_field%get_proxy()

    ndf_phys   = phys_p(1)%vspace%get_ndf()
    nodes_phys => phys_p(1)%vspace%get_nodes()

    nlayers = comp_p%vspace%get_nlayers()

    ndf_comp  = comp_p%vspace%get_ndf()
    undf_comp = comp_p%vspace%get_undf()
    dim_comp  = comp_p%vspace%get_dim_space()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )

    ! Evaluate the diff basis function
    allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf_phys))
    do df_phys = 1, ndf_phys
      do df_chi = 1, ndf_chi
        diff_basis_chi(:,df_chi,df_phys) = chi_p(1)%vspace%call_function(DIFF_BASIS,df_chi,nodes_phys(:,df_phys))
      end do
    end do

    ! Evaluate the basis function
    allocate(basis_comp(dim_comp, ndf_comp, ndf_phys))
    do df_phys = 1, ndf_phys
      do df_comp = 1, ndf_comp
        basis_comp(:,df_comp,df_phys) = comp_p%vspace%call_function(BASIS,df_comp,nodes_phys(:,df_phys))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
      call chi_p(1)%halo_exchange(depth=1)
    end if
    if (chi_p(2)%is_dirty(depth=1)) then
      call chi_p(2)%halo_exchange(depth=1)
    end if
    if (chi_p(3)%is_dirty(depth=1)) then
      call chi_p(3)%halo_exchange(depth=1)
    end if

    if (phys_p(1)%is_dirty(depth=1)) then
      call phys_p(1)%halo_exchange(depth=1)
    end if
    if (phys_p(2)%is_dirty(depth=1)) then
      call phys_p(2)%halo_exchange(depth=1)
    end if
    if (phys_p(3)%is_dirty(depth=1)) then
      call phys_p(3)%halo_exchange(depth=1)
    end if

    if (comp_p%is_dirty(depth=1)) then
      call comp_p%halo_exchange(depth=1)
    end if
    mesh => phys_field(1)%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map     => comp_p%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call convert_hcurl_field_code(nlayers, &
                                     phys_p(1)%data, &
                                     phys_p(2)%data, &
                                     phys_p(3)%data, &
                                     comp_p%data, &
                                     chi_p(1)%data, &
                                     chi_p(2)%data, &
                                     chi_p(3)%data, &
                                     ndf_comp, undf_comp, map, &
                                     ndf_chi, undf_chi, map_chi, &
                                     basis_comp, &
                                     diff_basis_chi &
                                    )
    end do

    call phys_p(1)%set_dirty()
    call phys_p(2)%set_dirty()
    call phys_p(3)%set_dirty()

    deallocate(diff_basis_chi, basis_comp)
  end subroutine invoke_convert_hcurl_field

  !-------------------------------------------------------------------------------
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data.
  !> PSyclone support is required - documented in #942
  subroutine invoke_convert_hdiv_field(phys_field, comp_field, chi )
    use convert_hdiv_field_kernel_mod, only: convert_hdiv_field_code
    use mesh_mod,                      only: mesh_type
    implicit none
    
    type(field_type), intent(inout)      :: phys_field(3)
    type(field_type), intent(in)         :: chi(3), comp_field

    type(field_proxy_type) :: phys_p(3), chi_p(3), comp_p
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_comp
    integer                 :: undf_chi, undf_comp
    integer                 :: diff_dim_chi, dim_comp
    integer                 :: df_comp1, df_comp2, df_chi
    integer, pointer        :: map(:) => null()
    integer, pointer        :: map_chi(:) => null()
    real(kind=r_def), pointer :: nodes_phys(:,:) => null()

    real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:), basis_comp(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()
    do i = 1,3
      phys_p(i) = phys_field(i)%get_proxy()
      chi_p(i)  = chi(i)%get_proxy()
    end do
    comp_p = comp_field%get_proxy()

    nodes_phys => phys_p(1)%vspace%get_nodes()

    nlayers = comp_p%vspace%get_nlayers()

    ndf_comp  = comp_p%vspace%get_ndf()
    undf_comp = comp_p%vspace%get_undf()
    dim_comp  = comp_p%vspace%get_dim_space()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )

    ! Evaluate the diff basis function
    allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf_comp))
    do df_comp2 = 1, ndf_comp
      do df_chi = 1, ndf_chi
        diff_basis_chi(:,df_chi,df_comp2) = chi_p(1)%vspace%call_function(DIFF_BASIS,df_chi,nodes_phys(:,df_comp2))
      end do
    end do

    ! Evaluate the diff basis function
    allocate(basis_comp(dim_comp, ndf_comp, ndf_comp) )
    do df_comp2 = 1, ndf_comp
      do df_comp1 = 1, ndf_comp
        basis_comp(:,df_comp1,df_comp2) = comp_p%vspace%call_function(BASIS,df_comp1,nodes_phys(:,df_comp2))
      end do
    end do


    if (chi_p(1)%is_dirty(depth=1)) then
      call chi_p(1)%halo_exchange(depth=1)
    end if
    if (chi_p(2)%is_dirty(depth=1)) then
      call chi_p(2)%halo_exchange(depth=1)
    end if
    if (chi_p(3)%is_dirty(depth=1)) then
      call chi_p(3)%halo_exchange(depth=1)
    end if

    if (phys_p(1)%is_dirty(depth=1)) then
      call phys_p(1)%halo_exchange(depth=1)
    end if
    if (phys_p(2)%is_dirty(depth=1)) then
      call phys_p(2)%halo_exchange(depth=1)
    end if
    if (phys_p(3)%is_dirty(depth=1)) then
      call phys_p(3)%halo_exchange(depth=1)
    end if

    if (comp_p%is_dirty(depth=1)) then
      call comp_p%halo_exchange(depth=1)
    end if
    mesh => comp_field%get_mesh()
    do cell = 1, mesh%get_last_halo_cell(1)
       map     => comp_p%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call convert_hdiv_field_code(nlayers, &
                                    phys_p(1)%data, &
                                    phys_p(2)%data, &
                                    phys_p(3)%data, &
                                    comp_p%data, &
                                    chi_p(1)%data, &
                                    chi_p(2)%data, &
                                    chi_p(3)%data, &
                                    ndf_comp, undf_comp, map, &
                                    ndf_chi, undf_chi, map_chi, &
                                    basis_comp, &
                                    diff_basis_chi &
                                   )
    end do

    call phys_p(1)%set_dirty()
    call phys_p(2)%set_dirty()
    call phys_p(3)%set_dirty()

    deallocate(diff_basis_chi, basis_comp)
  end subroutine invoke_convert_hdiv_field
!-------------------------------------------------------------------------------   
  subroutine invoke_convert_cart2sphere_vector( field, coords)
    use coord_transform_mod, only: cart2sphere_vector
    implicit none
    type(field_type), intent(inout) :: field(3)
    type(field_type), intent(in)    :: coords(3)

    type(field_proxy_type) :: f_p(3), x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: vector_in(3), vector_out(3), xyz(3)

    do i = 1,3
      f_p(i) = field(i)%get_proxy()
      x_p(i) = coords(i)%get_proxy()
    end do
    undf = f_p(1)%vspace%get_undf()
    do df = 1, undf
      vector_in(:)  = (/ f_p(1)%data(df), f_p(2)%data(df), f_p(3)%data(df) /)              
      xyz(:)        = (/ x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df) /)
      vector_out(:) = cart2sphere_vector(xyz, vector_in)
      f_p(1)%data(df) = vector_out(1)
      f_p(2)%data(df) = vector_out(2)
      f_p(3)%data(df) = vector_out(3)
    end do

    call f_p(1)%set_dirty()
    call f_p(2)%set_dirty()
    call f_p(3)%set_dirty()

  end subroutine invoke_convert_cart2sphere_vector
!-------------------------------------------------------------------------------   
  subroutine invoke_pointwise_convert_xyz2llr( coords)
    use coord_transform_mod, only: xyz2llr
    implicit none
    type(field_type), intent(inout) :: coords(3)

    type(field_proxy_type) :: x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: llr(3)

    do i = 1,3
      x_p(i) = coords(i)%get_proxy()
    end do
    undf = x_p(1)%vspace%get_undf()
    do df = 1, undf
      call xyz2llr(x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), &
                   llr(1), llr(2), llr(3))
      x_p(1)%data(df) = llr(1)
      x_p(2)%data(df) = llr(2)
      x_p(3)%data(df) = llr(3)
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

  end subroutine invoke_pointwise_convert_xyz2llr

!------------------------------------------------------------------------------- 
  subroutine invoke_compute_dof_level_kernel(level)

  use compute_dof_level_kernel_mod, only: compute_dof_level_code
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(inout) :: level
  type(field_proxy_type) :: l_p
  integer :: cell, ndf, undf
  real(kind=r_def), pointer :: nodes(:,:) => null()
  integer, pointer :: map(:) => null()
  type(mesh_type), pointer :: mesh => null()
  l_p = level%get_proxy()
  undf = l_p%vspace%get_undf()
  ndf  = l_p%vspace%get_ndf()
  nodes => l_p%vspace%get_nodes( )

  mesh => level%get_mesh()
  do cell = 1,mesh%get_last_halo_cell(1)
    map => l_p%vspace%get_cell_dofmap(cell)
    call compute_dof_level_code(l_p%vspace%get_nlayers(),                 &
                                l_p%data,                                 &
                                ndf,                                      &
                                undf,                                     &
                                map,                                      &
                                nodes                                     &
                               )
  end do 
  call l_p%set_dirty()

  end subroutine invoke_compute_dof_level_kernel

!------------------------------------------------------------------------------- 
subroutine invoke_write_fields(nodal_coordinates, level, nodal_output, fspace_dimension, output_unit, fname)
  use constants_mod,      only: str_max_filename

  implicit none
  type(field_type), intent(in) :: nodal_coordinates(3), level, nodal_output(3)
  integer,          intent(in) :: fspace_dimension, output_unit
  character(str_max_filename), intent(in) :: fname
  type(field_proxy_type) :: x_p(3), l_p, n_p(3)

  integer :: df, undf, i


  do i = 1,3
    x_p(i) = nodal_coordinates(i)%get_proxy()
    n_p(i) = nodal_output(i)%get_proxy()
  end do
  l_p = level%get_proxy()

  undf = n_p(1)%vspace%get_last_dof_owned()

  open(OUTPUT_UNIT, file = trim(fname), status = "replace")   
  write(OUTPUT_UNIT,'(A)') 'x = [' 
  if ( fspace_dimension  == 1 ) then
    do df = 1,undf
      write(OUTPUT_UNIT,'(5e16.8)') x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), l_p%data(df), n_p(1)%data(df)
    end do
  else
    do df = 1,undf
      write(OUTPUT_UNIT,'(7e16.8)') x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), l_p%data(df), &
                                    n_p(1)%data(df), n_p(2)%data(df), n_p(3)%data(df)
    end do
  end if
  write(OUTPUT_UNIT,'(A)') '];'
  close(OUTPUT_UNIT)
  
end subroutine invoke_write_fields

!-------------------------------------------------------------------------------  
!> invoke_subgrid_coeffs: Invoke the calculation of subgrid rho coefficients
subroutine invoke_subgrid_coeffs(a0,a1,a2,rho,cell_orientation,direction,rho_approximation_stencil_extent,halo_depth_to_compute)

    use flux_direction_mod,        only: x_direction, y_direction
    use stencil_dofmap_mod,        only: stencil_dofmap_type, &
                                         STENCIL_1DX,         &
                                         STENCIL_1DY
    use subgrid_coeffs_kernel_mod, only: subgrid_coeffs_code
    use subgrid_config_mod,        only: rho_approximation
    use mesh_mod,                  only: mesh_type
    use log_mod,                   only: log_event, LOG_LEVEL_ERROR

    implicit none

    type( field_type ), intent( inout ) :: a0
    type( field_type ), intent( inout ) :: a1
    type( field_type ), intent( inout ) :: a2
    type( field_type ), intent( in )    :: rho
    type( field_type ), intent( in )    :: cell_orientation
    integer, intent(in)                 :: direction
    integer, intent(in)                 :: rho_approximation_stencil_extent
    integer, intent(in)                 :: halo_depth_to_compute

    type( field_proxy_type )            :: rho_proxy
    type( field_proxy_type )            :: a0_proxy
    type( field_proxy_type )            :: a1_proxy
    type( field_proxy_type )            :: a2_proxy
    type( field_proxy_type )            :: cell_orientation_proxy

    type(stencil_dofmap_type), pointer  :: map => null()
    integer, pointer                    :: stencil_map(:,:) => null()
    integer                             :: rho_stencil_size
    integer                 :: cell
    integer                 :: nlayers
    integer                 :: ndf_w3
    integer                 :: undf_w3
    type(mesh_type), pointer :: mesh => null()
    integer                  :: d
    logical                  :: swap
    integer                  :: ncells_to_iterate

    a0_proxy   = a0%get_proxy()
    a1_proxy   = a1%get_proxy()
    a2_proxy   = a2%get_proxy()
    rho_proxy  = rho%get_proxy()
    cell_orientation_proxy = cell_orientation%get_proxy()

    undf_w3 = rho_proxy%vspace%get_undf()
    ndf_w3  = rho_proxy%vspace%get_ndf()
    nlayers = rho_proxy%vspace%get_nlayers()

    ! Note stencil grid types are of the form:
    !                                   |5|
    !                                   |3|
    ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
    !                                   |2|
    !                                   |4|
    if (direction == x_direction) then
      map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_approximation_stencil_extent)
    elseif (direction == y_direction) then
      map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_approximation_stencil_extent)
    end if
    rho_stencil_size = map%get_size()

    swap = .false.
    do d = 1,rho_approximation_stencil_extent
      if (rho_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_proxy%halo_exchange(depth=rho_approximation_stencil_extent)

    mesh => a0%get_mesh()
    if (halo_depth_to_compute==0) then
      ncells_to_iterate = mesh%get_last_edge_cell()
    elseif (halo_depth_to_compute > 0) then
      ncells_to_iterate = mesh%get_last_halo_cell(halo_depth_to_compute)
    else
      call log_event( "Error: negative halo_depth_to_compute value in subgrid coeffs call", LOG_LEVEL_ERROR )
    endif

    !NOTE: The default looping limits for this type of field would be 
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! inorder to function correctly. See ticket #1058.
    ! The kernel loops over all core and some halo cells.
    do cell = 1, ncells_to_iterate

      stencil_map => map%get_dofmap(cell)

      call subgrid_coeffs_code( nlayers,                                  &
                                rho_approximation,                        &
                                undf_w3,                                  &
                                rho_proxy%data,                           &
                                cell_orientation_proxy%data,              &
                                ndf_w3,                                   &
                                rho_stencil_size,                         &
                                stencil_map,                              &
                                direction,                                &
                                a0_proxy%data,                            &
                                a1_proxy%data,                            &
                                a2_proxy%data                             &
                                )

    end do
    call a0_proxy%set_dirty()
    call a1_proxy%set_dirty()
    call a2_proxy%set_dirty()

  end subroutine invoke_subgrid_coeffs


!------------------------------------------------------------------------------- 
subroutine invoke_conservative_fluxes(    rho,          &
                                          dep_pts,      &
                                          u_piola,      &
                                          mass_flux,    &
                                          a0_coeffs,    &
                                          a1_coeffs,    &
                                          a2_coeffs,    &
                                          direction,    &
                                          stencil_extent )

  use conservative_flux_kernel_mod, only: conservative_flux_code
  use flux_direction_mod,           only: x_direction, y_direction
  use stencil_dofmap_mod,           only: stencil_dofmap_type, &
                                          STENCIL_1DX, STENCIL_1DY
  use timestepping_config_mod,      only: dt
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(in)      :: rho
  type(field_type), intent(in)      :: dep_pts
  type(field_type), intent(in)      :: u_piola
  type(field_type), intent(inout)   :: mass_flux
  type(field_type), intent(in)      :: a0_coeffs
  type(field_type), intent(in)      :: a1_coeffs
  type(field_type), intent(in)      :: a2_coeffs
  integer, intent(in)               :: direction
  integer, intent(in)               :: stencil_extent

  type( field_proxy_type )  :: mass_flux_proxy, dep_pts_proxy, rho_proxy,     &
                               u_piola_proxy
  type( field_proxy_type )  :: a0_coeffs_proxy, a1_coeffs_proxy, a2_coeffs_proxy

  type(stencil_dofmap_type), pointer  :: map => null()

  integer, pointer :: map_rho(:) => null()
  integer, pointer :: map_w2(:) => null()
  integer, pointer :: stencil_map(:,:) => null()
  integer          :: stencil_size

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: cell
  integer :: nlayers
  type(mesh_type), pointer :: mesh => null()
  integer                  :: d
  logical                  :: swap

  rho_proxy     = rho%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()
  u_piola_proxy = u_piola%get_proxy()

  ndf_w3  = rho_proxy%vspace%get_ndf()
  undf_w3 = rho_proxy%vspace%get_undf()

  ndf_w2  = dep_pts_proxy%vspace%get_ndf()
  undf_w2 = dep_pts_proxy%vspace%get_undf()

  a0_coeffs_proxy = a0_coeffs%get_proxy()
  a1_coeffs_proxy = a1_coeffs%get_proxy()
  a2_coeffs_proxy = a2_coeffs%get_proxy()
  mass_flux_proxy = mass_flux%get_proxy()

  nlayers = rho_proxy%vspace%get_nlayers()

  ! Note stencil grid types are of the form:
  !                                   |5|
  !                                   |3|
  ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
  !                                   |2|
  !                                   |4|
  if (direction == x_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,stencil_extent)
  elseif (direction == y_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,stencil_extent)
  end if
  stencil_size = map%get_size()

  swap = .false.
  do d = 1,stencil_extent
    if (rho_proxy%is_dirty(depth=d)) swap = .true.
  end do
  if ( swap ) call rho_proxy%halo_exchange(depth=stencil_extent)

  mesh => rho%get_mesh()
  !NOTE: The default looping limits for this type of field would be 
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! inorder to function correctly. See ticket #1058.
  ! The kernel loops over all core cells only.
  do cell = 1, mesh%get_last_edge_cell() 
      map_rho => rho_proxy%vspace%get_cell_dofmap( cell )
      map_w2 => dep_pts_proxy%vspace%get_cell_dofmap( cell )

      stencil_map => map%get_dofmap(cell)

      call conservative_flux_code( nlayers,                     &
                                   undf_w3,                     &
                                   ndf_w3,                      &
                                   map_rho,                     &
                                   rho_proxy%data,              &
                                   a0_coeffs_proxy%data,        &
                                   a1_coeffs_proxy%data,        &
                                   a2_coeffs_proxy%data,        &
                                   undf_w2,                     &
                                   ndf_w2,                      &
                                   map_w2,                      &
                                   mass_flux_proxy%data,        &
                                   dep_pts_proxy%data,          &
                                   u_piola_proxy%data,          &
                                   stencil_size,                &
                                   stencil_map,                 &
                                   direction,                   &
                                   dt )

  end do
  call a0_coeffs_proxy%set_dirty()
  call a1_coeffs_proxy%set_dirty()
  call a2_coeffs_proxy%set_dirty()


end subroutine invoke_conservative_fluxes

!------------------------------------------------------------------------------- 
subroutine invoke_fv_mass_fluxes( rho,            &
                                  dep_pts,        &
                                  u_piola,        &
                                  mass_flux,      &
                                  a0_coeffs,      &
                                  a1_coeffs,      &
                                  a2_coeffs,      &
                                  direction,      &
                                  stencil_extent )

  use fv_mass_flux_kernel_mod,      only: fv_mass_flux_code
  use flux_direction_mod,           only: x_direction, y_direction
  use stencil_dofmap_mod,           only: stencil_dofmap_type, &
                                          STENCIL_1DX, STENCIL_1DY
  use timestepping_config_mod,      only: dt
  use mesh_mod,                     only: mesh_type
  implicit none

  type(field_type), intent(in)      :: rho
  type(field_type), intent(in)      :: dep_pts
  type(field_type), intent(in)      :: u_piola
  type(field_type), intent(inout)   :: mass_flux
  type(field_type), intent(in)      :: a0_coeffs
  type(field_type), intent(in)      :: a1_coeffs
  type(field_type), intent(in)      :: a2_coeffs
  integer, intent(in)               :: direction
  integer, intent(in)               :: stencil_extent

  type( field_proxy_type )  :: mass_flux_proxy, dep_pts_proxy, rho_proxy,     &
                               u_piola_proxy
  type( field_proxy_type )  :: a0_coeffs_proxy, a1_coeffs_proxy, a2_coeffs_proxy

  type(stencil_dofmap_type), pointer  :: map => null()

  integer, pointer :: map_rho(:) => null()
  integer, pointer :: map_w2(:) => null()
  integer, pointer :: stencil_map(:,:) => null()
  integer          :: stencil_size

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: cell
  integer :: nlayers
  type(mesh_type), pointer :: mesh => null()
  integer                  :: d
  logical                  :: swap

  rho_proxy     = rho%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()
  u_piola_proxy = u_piola%get_proxy()

  ndf_w3  = rho_proxy%vspace%get_ndf()
  undf_w3 = rho_proxy%vspace%get_undf()

  ndf_w2  = dep_pts_proxy%vspace%get_ndf()
  undf_w2 = dep_pts_proxy%vspace%get_undf()

  a0_coeffs_proxy = a0_coeffs%get_proxy()
  a1_coeffs_proxy = a1_coeffs%get_proxy()
  a2_coeffs_proxy = a2_coeffs%get_proxy()
  mass_flux_proxy = mass_flux%get_proxy()

  nlayers = rho_proxy%vspace%get_nlayers()

  ! Note stencil grid types are of the form:
  !                                   |5|
  !                                   |3|
  ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
  !                                   |2|
  !                                   |4|
  if (direction == x_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,stencil_extent)
  elseif (direction == y_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,stencil_extent)
  end if
  stencil_size = map%get_size()

  swap = .false.
  do d = 1,stencil_extent
    if (rho_proxy%is_dirty(depth=d)) swap = .true.
  end do
  if ( swap ) call rho_proxy%halo_exchange(depth=stencil_extent)

  mesh => rho%get_mesh()
  ! NOTE: The default looping limits for this type of field would be 
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! in order to function correctly. See ticket #1058.
  ! The kernel loops over all core cells only.
  do cell = 1, mesh%get_last_edge_cell() 
      map_rho => rho_proxy%vspace%get_cell_dofmap( cell )
      map_w2 => dep_pts_proxy%vspace%get_cell_dofmap( cell )

      stencil_map => map%get_dofmap(cell)

      call fv_mass_flux_code(  nlayers,                     &
                               undf_w3,                     &
                               ndf_w3,                      &
                               map_rho,                     &
                               rho_proxy%data,              &
                               a0_coeffs_proxy%data,        &
                               a1_coeffs_proxy%data,        &
                               a2_coeffs_proxy%data,        &
                               undf_w2,                     &
                               ndf_w2,                      &
                               map_w2,                      &
                               mass_flux_proxy%data,        &
                               dep_pts_proxy%data,          &
                               u_piola_proxy%data,          &
                               stencil_size,                &
                               stencil_map,                 &
                               direction,                   &
                               dt )

  end do
  call a0_coeffs_proxy%set_dirty()
  call a1_coeffs_proxy%set_dirty()
  call a2_coeffs_proxy%set_dirty()

end subroutine invoke_fv_mass_fluxes

!-------------------------------------------------------------------------------
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
subroutine invoke_calc_departure_wind(u_departure_wind, u_piola, chi )
  use calc_departure_wind_kernel_mod, only: calc_departure_wind_code
  use mesh_mod,                       only: mesh_type
  implicit none

  type(field_type), intent(inout)      :: u_departure_wind
  type(field_type), intent(in)         :: chi(3), u_piola

  type(field_proxy_type) :: u_departure_wind_p, chi_p(3), u_piola_p

  integer                 :: cell, nlayers
  integer                 :: ndf_chi, ndf_u, ndf_udep
  integer                 :: undf_chi, undf_u
  integer                 :: dim_u, diff_dim_chi
  integer                 :: df_u, df_chi, df_udep

  integer, pointer        :: map_chi(:) => null()
  integer, pointer        :: map(:) => null()
  real(kind=r_def), pointer :: nodes_udep(:,:) => null()

  real(kind=r_def), allocatable  :: nodal_basis_u(:,:,:)
  real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:)
  integer :: ii
  type(mesh_type), pointer :: mesh => null()

  do ii = 1,3
    chi_p(ii)  = chi(ii)%get_proxy()
  end do
  u_piola_p = u_piola%get_proxy()
  u_departure_wind_p = u_departure_wind%get_proxy()

  nlayers = u_piola_p%vspace%get_nlayers()

  ndf_udep   = u_departure_wind_p%vspace%get_ndf()
  nodes_udep => u_departure_wind_p%vspace%get_nodes()

  ndf_u  = u_piola_p%vspace%get_ndf()
  undf_u = u_piola_p%vspace%get_undf()
  dim_u = u_piola_p%vspace%get_dim_space()

  ndf_chi  = chi_p(1)%vspace%get_ndf()
  undf_chi = chi_p(1)%vspace%get_undf()
  diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff()

  ! Evaluate the basis function
  allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf_u))
  do df_udep = 1, ndf_udep
    do df_chi = 1, ndf_chi
      diff_basis_chi(:,df_chi,df_udep) = chi_p(1)%vspace%call_function(DIFF_BASIS,df_chi,nodes_udep(:,df_udep))
    end do
  end do

  ! Evaluate the basis function
  allocate(nodal_basis_u(dim_u, ndf_u, ndf_u))
  do df_udep = 1, ndf_udep
    do df_u = 1, ndf_u
      nodal_basis_u(:,df_u,df_udep) = u_piola_p%vspace%call_function(BASIS,df_u,nodes_udep(:,df_udep))
    end do
  end do


  if (u_piola_p%is_dirty(depth=1)) call u_piola_p%halo_exchange(depth=1)
  if (chi_p(1)%is_dirty(depth=1))  call chi_p(1)%halo_exchange(depth=1)
  if (chi_p(2)%is_dirty(depth=1))  call chi_p(2)%halo_exchange(depth=1)
  if (chi_p(3)%is_dirty(depth=1))  call chi_p(3)%halo_exchange(depth=1)

  mesh => u_piola%get_mesh()
  !NOTE: The default looping limits for this type of field would be 
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! inorder to function correctly. See ticket #1058.
  ! The kernel loops over all core and all halo cells.
  do cell = 1,mesh%get_ncells_2d()
     map     => u_piola_p%vspace%get_cell_dofmap( cell )
     map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
     call calc_departure_wind_code( nlayers,                                  &
                                    u_departure_wind_p%data,                  &
                                    u_piola_p%data,                           &
                                    chi_p(1)%data,                            &
                                    chi_p(2)%data,                            &
                                    chi_p(3)%data,                            &
                                    ndf_u, undf_u, map, nodal_basis_u,        &
                                    ndf_chi, undf_chi, map_chi,               &
                                    diff_basis_chi                            &
                                     )
  end do
  deallocate(nodal_basis_u)
  deallocate(diff_basis_chi)
  call u_departure_wind_p%set_dirty()
end subroutine invoke_calc_departure_wind
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> invoke_calc_deppts: Invoke the calculation of departure points in 1D
subroutine invoke_calc_deppts(  u_n,                  &
                                u_np1,                &
                                dep_pts,              &
                                cell_orientation,     &
                                direction,            &
                                dep_pt_method,        &
                                dep_pt_stencil_extent )

  use calc_departure_point_kernel_mod,  only : calc_departure_point_code
  use stencil_dofmap_mod,               only : stencil_dofmap_type, &
                                               STENCIL_1DX, &
                                               STENCIL_1DY
  use flux_direction_mod,               only : x_direction, y_direction
  use mesh_mod,                         only : mesh_type
  implicit none

  type( field_type ), intent( in )    :: u_n
  type( field_type ), intent( in )    :: u_np1
  type( field_type ), intent( inout ) :: dep_pts
  type( field_type ), intent( in )    :: cell_orientation
  integer, intent(in)                 :: direction
  integer, intent(in)                 :: dep_pt_method
  integer, intent(in)                 :: dep_pt_stencil_extent

  type( field_proxy_type )        :: u_n_proxy
  type( field_proxy_type )        :: u_np1_proxy
  type( field_proxy_type )        :: dep_pts_proxy
  type( field_proxy_type )        :: cell_orientation_proxy
  type(stencil_dofmap_type), pointer  :: map=>null()
  type(stencil_dofmap_type), pointer  :: map_w3=>null()

  integer, pointer        :: stencil_map_w2(:,:) => null()
  integer, pointer        :: stencil_map_w3(:,:) => null()
  integer                 :: transport_stencil_size

  integer                 :: cell
  integer                 :: nlayers
  integer                 :: ndf_w2
  integer                 :: undf_w2
  integer                 :: ndf_w3
  integer                 :: undf_w3
  type(mesh_type), pointer :: mesh => null()
  integer                  :: d
  logical                  :: swap
  
  u_n_proxy    = u_n%get_proxy()
  u_np1_proxy  = u_np1%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()
  cell_orientation_proxy = cell_orientation%get_proxy()

  ndf_w2  = u_n_proxy%vspace%get_ndf()
  undf_w2 = u_n_proxy%vspace%get_undf()

  ndf_w3  =   cell_orientation_proxy%vspace%get_ndf()
  undf_w3 =   cell_orientation_proxy%vspace%get_undf()

  if (direction == x_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,dep_pt_stencil_extent)
    map_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,dep_pt_stencil_extent)
  elseif (direction == y_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,dep_pt_stencil_extent)
    map_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,dep_pt_stencil_extent)
  endif
  transport_stencil_size = map%get_size()

  nlayers = u_n_proxy%vspace%get_nlayers()

  swap = .false.
  do d = 1,dep_pt_stencil_extent
    if (u_n_proxy%is_dirty(depth=d)) swap = .true.
  end do
  if ( swap ) call u_n_proxy%halo_exchange(depth=dep_pt_stencil_extent)
  swap = .false.
  do d = 1,dep_pt_stencil_extent
    if (u_np1_proxy%is_dirty(depth=d)) swap = .true.
  end do
  if ( swap ) call u_np1_proxy%halo_exchange(depth=dep_pt_stencil_extent)

  mesh => u_n%get_mesh()
  !NOTE: The default looping limits for this type of field would be 
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! inorder to function correctly. See ticket #1058.
  ! The kernel loops over all core cells only.
  do cell=1,mesh%get_last_edge_cell()

    stencil_map_w2 => map%get_dofmap(cell)
    stencil_map_w3 => map_w3%get_dofmap(cell)

    call calc_departure_point_code( nlayers,                      &
                                    dep_pts_proxy%data,           &
                                    transport_stencil_size,       &
                                    undf_w2,                      &
                                    ndf_w2,                       &
                                    stencil_map_w2,               &
                                    undf_w3,                      &
                                    ndf_w3,                       &
                                    stencil_map_w3,               &
                                    cell_orientation_proxy%data,  &
                                    u_n_proxy%data,               &
                                    u_np1_proxy%data,             &
                                    direction,                    &
                                    dep_pt_method )

  end do
  call dep_pts_proxy%set_dirty()

end subroutine invoke_calc_deppts

!-------------------------------------------------------------------------------   
!> invoke_inc_axpy: x = a * x + y ; a-scalar, x,y-vector     
  subroutine invoke_inc_axpy(scalar, field1, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    real(kind=r_def),   intent(in)     :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inc_axpy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, undf, scalar),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = scalar * field1_proxy%data(i) + field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_inc_axpy
  
!-------------------------------------------------------------------------------   
!> invoke_increment_field: y = y + x ; x,y-vector     
  subroutine invoke_increment_field(x, y)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: y
    type( field_type ), intent(in)     :: x
    type( field_proxy_type)            :: x_proxy, y_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    x_proxy = x%get_proxy()
    y_proxy = y%get_proxy()

    !sanity check
    undf = x_proxy%vspace%get_undf()
    if(undf /= y_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:increment_field:field1 and field2 live on different spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(y_proxy, x_proxy, undf),  private(i)
    do i = 1,undf
      y_proxy%data(i) = y_proxy%data(i) + x_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => y%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( y_proxy%is_dirty(depth=dplp) .or. &
          x_proxy%is_dirty(depth=dplp) ) then
        call y_proxy%set_dirty()
      else
        call y_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_increment_field
  
!-------------------------------------------------------------------------------   

!> invoke_divide_field_data: Divide the values of field1 by field2
!> c = a/b
  subroutine invoke_divide_field_data(field1, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:divide_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, undf),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = field1_proxy%data(i)/field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_divide_field_data

!-------------------------------------------------------------------------------   
!> invoke_multiply_field_data:  z =  x * y     
  subroutine invoke_multiply_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy,     &
                                          field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:multiply_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:multiply_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) * field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_multiply_field_data
!-------------------------------------------------------------------------------   
!> invoke_inc_axpby: x = a * x + b * y ; a,b-scalar, x,y-vector     
  subroutine invoke_inc_axpby(scalar1, field1, scalar2, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    real(kind=r_def),   intent(in)     :: scalar1, scalar2
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inc_axpby:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, scalar1, scalar2, undf),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = scalar1 * field1_proxy%data(i) &
                           + scalar2 * field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_inc_axpby
  
!-------------------------------------------------------------------------------   
!> invoke_inc_xpby: x = x + b * y ; b-scalar, x,y-vector     
  subroutine invoke_inc_xpby(field1, scalar2, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    real(kind=r_def),   intent(in)     :: scalar2
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inc_axpby:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, scalar2, undf),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = field1_proxy%data(i) &
                           + scalar2*field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_inc_xpby
  
!-------------------------------------------------------------------------------   
!> invoke_scale_field_data: scale data from one field ( a = scalar*a )
  subroutine invoke_scale_field_data(scalar,field)
    implicit none
    real( kind=r_def ), intent(in)     :: scalar
    type( field_type ), intent(in )    :: field
    type( field_proxy_type)            :: field_proxy
    integer(kind=i_def)                :: i,undf

    field_proxy = field%get_proxy()

    undf = field_proxy%vspace%get_undf()

    !$omp parallel do schedule(static), default(none), shared(field_proxy, undf, scalar),  private(i)
    do i = 1,undf
       field_proxy%data(i) = scalar*field_proxy%data(i)
    end do
    !$omp end parallel do
   
  end subroutine invoke_scale_field_data


!-------------------------------------------------------------------------------   

!-------------------------------------------------------------------------------
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
subroutine invoke_sample_poly_flux( flux, wind, density, stencil_extent )

  use sample_poly_flux_kernel_mod, only: sample_poly_flux_code
  use stencil_dofmap_mod,          only: stencil_dofmap_type, STENCIL_CROSS
  use mesh_mod,                    only: mesh_type
  implicit none

  type(field_type), intent(inout)      :: flux
  type(field_type), intent(in)         :: wind
  type(field_type), intent(in)         :: density
  integer, intent(in)                  :: stencil_extent

  type( field_proxy_type )  :: flux_proxy, wind_proxy, density_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()

  integer, pointer :: map_w2(:,:)        => null()
  integer, pointer :: stencil_map(:,:,:) => null()

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: df_w2, df_wind
  integer :: cell
  integer :: nlayers
  integer :: stencil_size 
  integer :: d
  logical :: swap
  type(mesh_type), pointer :: mesh => null()
  real(kind=r_def), pointer :: nodes_w2(:,:) => null()
  real(kind=r_def), allocatable :: basis_w2(:,:,:)
  integer :: dim_w2

  flux_proxy    = flux%get_proxy()
  wind_proxy    = wind%get_proxy()
  density_proxy = density%get_proxy()

  ndf_w3  = density_proxy%vspace%get_ndf()
  undf_w3 = density_proxy%vspace%get_undf()

  ndf_w2  = flux_proxy%vspace%get_ndf()
  undf_w2 = flux_proxy%vspace%get_undf()
  dim_w2  = flux_proxy%vspace%get_dim_space()
  nodes_w2 => flux_proxy%vspace%get_nodes()

  ! Evaluate the basis function
  allocate (basis_w2(dim_w2, ndf_w2, ndf_w2))
  do df_w2 = 1, ndf_w2
    do df_wind = 1, ndf_w2
      basis_w2(:,df_wind,df_w2) = wind_proxy%vspace%call_function(BASIS,df_wind,nodes_w2(:,df_w2))
    end do
  end do

  nlayers = flux_proxy%vspace%get_nlayers()

  stencil => density_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                     stencil_extent)
  stencil_size = stencil%get_size()
  stencil_map => stencil%get_whole_dofmap()
  map_w2 => flux_proxy%vspace%get_whole_dofmap()

  if(wind_proxy%is_dirty(depth=1) ) then
    call wind_proxy%halo_exchange(depth=1)
  end if
  swap = .false.
  do d = 1,stencil_extent+1
    if(density_proxy%is_dirty(depth=d) ) swap = .true.
  end do
  if ( swap )  call density_proxy%halo_exchange(depth=stencil_extent+1)
  if(flux_proxy%is_dirty(depth=1) ) then
    call flux_proxy%halo_exchange(depth=1)
  end if

  mesh => flux%get_mesh()
  do cell = 1,mesh%get_last_halo_cell(1)

      call sample_poly_flux_code( nlayers,                     &
                                  flux_proxy%data,             &
                                  wind_proxy%data,             &
                                  density_proxy%data,          &
                                  ndf_w2,                      &
                                  undf_w2,                     &
                                  map_w2(:,cell),              &
                                  basis_w2,                    &
                                  ndf_w3,                      &
                                  undf_w3,                     &
                                  stencil_size,                &
                                  stencil_map(:,:,cell)        &
                                  )

  end do

  call flux_proxy%set_dirty()

end subroutine invoke_sample_poly_flux
!-------------------------------------------------------------------------------
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
subroutine invoke_sample_poly_adv( adv, tracer, wind, stencil_extent )

  use sample_poly_adv_kernel_mod, only: sample_poly_adv_code
  use stencil_dofmap_mod,          only: stencil_dofmap_type, STENCIL_CROSS
  use mesh_mod,                    only: mesh_type

  implicit none

  type(field_type), intent(inout)      :: adv
  type(field_type), intent(in)         :: wind
  type(field_type), intent(in)         :: tracer
  integer, intent(in)                  :: stencil_extent

  type( field_proxy_type )  :: adv_proxy, wind_proxy, tracer_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()

  integer, pointer :: map_w2(:,:)        => null()
  integer, pointer :: stencil_map(:,:,:) => null()

  integer :: undf_wt, ndf_wt
  integer :: undf_w2, ndf_w2
  integer :: ndf_adv
  integer :: df_adv, df_w2
  integer :: cell
  integer :: nlayers
  integer :: stencil_size
  integer :: d
  logical :: swap
  type(mesh_type), pointer :: mesh => null()
  real(kind=r_def), pointer :: nodes_adv(:,:) => null()

  real(kind=r_def), allocatable :: basis_w2(:,:,:)
  integer :: dim_w2

  adv_proxy    = adv%get_proxy()
  wind_proxy   = wind%get_proxy()
  tracer_proxy = tracer%get_proxy()

  ndf_adv = adv_proxy%vspace%get_ndf()
  nodes_adv => adv_proxy%vspace%get_nodes()

  ndf_wt  = tracer_proxy%vspace%get_ndf()
  undf_wt = tracer_proxy%vspace%get_undf()

  ndf_w2  = wind_proxy%vspace%get_ndf()
  undf_w2 = wind_proxy%vspace%get_undf()
  dim_w2  = wind_proxy%vspace%get_dim_space()

  ! Evaluate the basis function
  allocate (basis_w2(dim_w2, ndf_w2, ndf_adv))
  do df_adv = 1, ndf_adv
    do df_w2 = 1, ndf_w2
      basis_w2(:,df_w2,df_adv) = wind_proxy%vspace%call_function(BASIS,df_w2,nodes_adv(:,df_adv))
    end do
  end do

  nlayers = adv_proxy%vspace%get_nlayers()

  stencil => tracer_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                    stencil_extent)
  stencil_size = stencil%get_size()
  stencil_map => stencil%get_whole_dofmap()
  map_w2 => wind_proxy%vspace%get_whole_dofmap()

  if(wind_proxy%is_dirty(depth=1) ) then
    call wind_proxy%halo_exchange(depth=1)
  end if
  swap = .false.
  do d = 1,stencil_extent+1
    if(tracer_proxy%is_dirty(depth=d) ) swap = .true.
  end do
  if ( swap ) call tracer_proxy%halo_exchange(depth=stencil_extent+1)

  mesh => adv%get_mesh()
  do cell = 1,mesh%get_last_edge_cell()

      call sample_poly_adv_code( nlayers,                     &
                                 adv_proxy%data,              &
                                 tracer_proxy%data,           &
                                 wind_proxy%data,             &
                                 ndf_wt,                      &
                                 undf_wt,                     &
                                 ndf_w2,                      &
                                 undf_w2,                     &
                                 map_w2(:,cell),              &
                                 basis_w2,                    &
                                 stencil_size,                &
                                 stencil_map(:,:,cell)        &
                                 )

  end do

  call adv_proxy%set_dirty()

end subroutine invoke_sample_poly_adv
!------------------------------------------------------------------------------- 
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
    subroutine invoke_compute_tri_precon_kernel(tri_precon, theta, rho, chi, m3_inv )
      use compute_tri_precon_kernel_mod, only: compute_tri_precon_code
      use mesh_mod, only: mesh_type
      type(field_type), intent(inout) :: tri_precon(3)
      type(field_type), intent(in) :: theta, rho, chi(3)
      type(operator_type), intent(in) :: m3_inv

      integer, pointer :: map_w3(:) => null(), map_wtheta(:) => null(), map_any_space_1_chi(:) => null()
      real(kind=r_def), pointer :: nodes_w3(:,:) => null()
      integer :: cell
      integer :: ndf_w3, undf_w3, ndf_wtheta, undf_wtheta, ndf_any_space_1_chi, undf_any_space_1_chi
      integer :: df_w3, df_chi
      type(mesh_type), pointer :: mesh => null()
      integer :: nlayers
      type(field_proxy_type) :: tri_precon_proxy(3), theta_proxy, rho_proxy, chi_proxy(3)
      type(operator_proxy_type) :: m3_inv_proxy
      real(kind=r_def), allocatable :: diff_basis_chi(:,:,:)
      integer :: diff_dim_chi
      !
      ! Initialise field proxies
      !
      tri_precon_proxy(1) = tri_precon(1)%get_proxy()
      tri_precon_proxy(2) = tri_precon(2)%get_proxy()
      tri_precon_proxy(3) = tri_precon(3)%get_proxy()
      theta_proxy = theta%get_proxy()
      rho_proxy = rho%get_proxy()
      chi_proxy(1) = chi(1)%get_proxy()
      chi_proxy(2) = chi(2)%get_proxy()
      chi_proxy(3) = chi(3)%get_proxy()
      m3_inv_proxy = m3_inv%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = tri_precon_proxy(1)%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => tri_precon(1)%get_mesh()
      !
      ! Initialise sizes and allocate any basis arrays for w3
      !
      ndf_w3 = tri_precon_proxy(1)%vspace%get_ndf()
      undf_w3 = tri_precon_proxy(1)%vspace%get_undf()
      nodes_w3 => tri_precon_proxy(1)%vspace%get_nodes()
      !
      ! Initialise sizes and allocate any basis arrays for w0
      !
      ndf_wtheta = theta_proxy%vspace%get_ndf()
      undf_wtheta = theta_proxy%vspace%get_undf()
      !
      ! Initialise sizes and allocate any basis arrays for any_space_1_chi
      !
      ndf_any_space_1_chi = chi_proxy(1)%vspace%get_ndf()
      undf_any_space_1_chi = chi_proxy(1)%vspace%get_undf()

      ! Compute nodal basis functions
      diff_dim_chi  = chi_proxy(1)%vspace%get_dim_space_diff( )

      ! Evaluate the basis function
      allocate( diff_basis_chi(diff_dim_chi, ndf_any_space_1_chi, ndf_w3) )
      do df_w3 = 1, ndf_w3
        do df_chi = 1, ndf_any_space_1_chi
          diff_basis_chi(:,df_chi,df_w3) = chi_proxy(1)%vspace%call_function(DIFF_BASIS,df_chi,nodes_w3(:,df_w3))
        end do
      end do

      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=1)) then
        call theta_proxy%halo_exchange(depth=1)
      end if 
      !
      if (rho_proxy%is_dirty(depth=1)) then
        call rho_proxy%halo_exchange(depth=1)
      end if 
      !
      if (chi_proxy(1)%is_dirty(depth=1)) then
        CALL chi_proxy(1)%halo_exchange(depth=1)
      end if 
      !
      if (chi_proxy(2)%is_dirty(depth=1)) then
        CALL chi_proxy(2)%halo_exchange(depth=1)
      end if 
      !
      if (chi_proxy(3)%is_dirty(depth=1)) then
        CALL chi_proxy(3)%halo_exchange(depth=1)
      end if 
      !
      do cell=1,mesh%get_last_edge_cell()
        !
        map_w3 => tri_precon_proxy(1)%vspace%get_cell_dofmap(cell)
        map_wtheta => theta_proxy%vspace%get_cell_dofmap(cell)
        map_any_space_1_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)
        !
        CALL compute_tri_precon_code(cell, nlayers, &
                                     tri_precon_proxy(1)%data, &
                                     tri_precon_proxy(2)%data, &
                                     tri_precon_proxy(3)%data, &
                                     theta_proxy%data, &
                                     rho_proxy%data, &
                                     chi_proxy(1)%data, &
                                     chi_proxy(2)%data, &
                                     chi_proxy(3)%data, &
                                     m3_inv_proxy%ncell_3d, &
                                     m3_inv_proxy%local_stencil, &
                                     ndf_w3, &
                                     undf_w3, &
                                     map_w3, &
                                     ndf_wtheta, &
                                     undf_wtheta, &
                                     map_wtheta, &
                                     ndf_any_space_1_chi, &
                                     undf_any_space_1_chi, &
                                     map_any_space_1_chi, &
                                     diff_basis_chi)
      end do 
      !
      ! Set halos dirty for fields modified in the above loop
      !
      call tri_precon_proxy(1)%set_dirty()
      call tri_precon_proxy(2)%set_dirty()
      call tri_precon_proxy(3)%set_dirty()
       
       
    end subroutine invoke_compute_tri_precon_kernel

!------------------------------------------------------------------------------- 
!> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
!> (not quadrature) will be removed and the functionality will be implemented via
!> PSY using kernal meta data. 
!> PSyclone support is required - documented in #942
!> invoke_initial_buoyancy_kernel: invoke the buoyancy initialization
  subroutine invoke_initial_buoyancy_kernel( b, chi )

    use initial_buoyancy_kernel_mod, only : initial_buoyancy_code
    use mesh_mod,                    only : mesh_type

    implicit none

    type( field_type ), intent( inout ) :: b
    type( field_type ), intent( in )    :: chi(3)

    integer          :: cell
    integer          :: ndf_wt, undf_wt, &
                        ndf_chi, undf_chi, dim_chi, &
                        df_wt, df_chi

    integer, pointer        :: map_wt(:) => null()
    integer, pointer        :: map_chi(:)    => null()
    real(kind=r_def), pointer :: nodes_wt(:,:) => null()


    type( field_proxy_type )        :: b_proxy
    type( field_proxy_type )        :: chi_proxy(3)

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    type(mesh_type),  pointer :: mesh => null()

    b_proxy  = b%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    mesh => b%get_mesh()

    ndf_wt  = b_proxy%vspace%get_ndf( )
    undf_wt = b_proxy%vspace%get_undf( )
    nodes_wt => b_proxy%vspace%get_nodes()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi  = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate( basis_chi(dim_chi, ndf_chi, ndf_wt) )
    do df_wt = 1, ndf_wt
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_wt) = chi_proxy(1)%vspace%call_function(BASIS,df_chi,nodes_wt(:,df_wt))
      end do
    end do


    if (chi_proxy(1)%is_dirty(depth=1)) then
      call chi_proxy(1)%halo_exchange(depth=1)
    end if
    if (chi_proxy(2)%is_dirty(depth=1)) then
      call chi_proxy(2)%halo_exchange(depth=1)
    end if
    if (chi_proxy(3)%is_dirty(depth=1)) then
      call chi_proxy(3)%halo_exchange(depth=1)
    end if
 
    do cell = 1,mesh%get_last_edge_cell()

      map_wt => b_proxy%vspace%get_cell_dofmap( cell )
      map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      call initial_buoyancy_code(       &
        b_proxy%vspace%get_nlayers(),   &
        b_proxy%data,                   &
        chi_proxy(1)%data,              &
        chi_proxy(2)%data,              &
        chi_proxy(3)%data,              &
        ndf_wt,                         &
        undf_wt,                        &
        map_wt,                         &
        ndf_chi,                        &
        undf_chi,                       &
        map_chi,                        &
        basis_chi                       &
        )
    end do

    call b_proxy%set_dirty()
  end subroutine invoke_initial_buoyancy_kernel
  !------------------------------------------------------------------------------
! Needs stencil support and reformatting of kernel to use stencils
  subroutine invoke_viscosity( theta_inc, theta_n, u_inc, u_n, chi )
   
    use viscosity_kernel_mod, only : viscosity_code
    use mesh_mod,             only : mesh_type 
    use stencil_dofmap_mod,   only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( in ) :: theta_inc, theta_n, u_inc, u_n, chi(3)
    type( mesh_type ), pointer       :: mesh => null()
    integer                 :: cell, nlayers
    integer                 :: ndf_w0, ndf_w2, ndf_chi
    integer                 :: undf_w0, undf_w2, undf_chi

    type( field_proxy_type )  :: theta_inc_proxy, theta_n_proxy, &
                                 u_inc_proxy, u_n_proxy, &
                                 chi_proxy(3)

    integer, pointer        :: map_chi(:)    => null()

    type(stencil_dofmap_type), pointer :: cross_stencil_w2 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_w0 => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer                 :: cross_stencil_w2_size

    integer, pointer        :: cross_stencil_w0_map(:,:,:) => null()
    integer                 :: cross_stencil_w0_size

    theta_inc_proxy = theta_inc%get_proxy()
    theta_n_proxy   = theta_n%get_proxy()
    u_inc_proxy = u_inc%get_proxy()
    u_n_proxy   = u_n%get_proxy()

    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    nlayers = theta_inc_proxy%vspace%get_nlayers()
    ndf_w0  = theta_inc_proxy%vspace%get_ndf( )
    undf_w0 = theta_inc_proxy%vspace%get_undf()

    ndf_w2  = u_inc_proxy%vspace%get_ndf( )
    undf_w2 = u_inc_proxy%vspace%get_undf()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf()

    mesh => theta_inc%get_mesh()

    if (theta_n_proxy%is_dirty(depth=1)) CALL theta_n_proxy%halo_exchange(depth=2)
    if (    u_n_proxy%is_dirty(depth=1)) CALL     u_n_proxy%halo_exchange(depth=2)

    cross_stencil_w2 => u_inc_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w2_map => cross_stencil_w2%get_whole_dofmap()
    cross_stencil_w2_size = cross_stencil_w2%get_size()

    cross_stencil_w0 => theta_inc_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w0_map => cross_stencil_w0%get_whole_dofmap()
    cross_stencil_w0_size = cross_stencil_w0%get_size()
    
    do cell = 1, mesh%get_last_halo_cell(1)
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

       call viscosity_code(nlayers,                   &
                           theta_inc_proxy%data,      &
                           theta_n_proxy%data,        &
                           u_inc_proxy%data,          &
                           u_n_proxy%data,            &
                           chi_proxy(1)%data,         &
                           chi_proxy(2)%data,         &
                           chi_proxy(3)%data,         &
                           ndf_w0, undf_w0,           &
                           cross_stencil_w0_map(:,:,cell), &
                           cross_stencil_w0_size,          &
                           ndf_w2, undf_w2,           &
                           cross_stencil_w2_map(:,:,cell), &
                           cross_stencil_w2_size,          &
                           ndf_chi, undf_chi, map_chi &
                          )
    end do
     
   call theta_inc_proxy%set_dirty()
   call     u_inc_proxy%set_dirty()

  end subroutine invoke_viscosity

!> invoke_raise_field: Raise a field to an integer exponent
!> x => x^a
  subroutine invoke_raise_field(x, a)
    implicit none
    type( field_type ), intent(inout)  :: x
    integer(i_def),     intent(in)     :: a
    type( field_proxy_type)            :: xp
    integer(kind=i_def)                :: i,undf

    xp = x%get_proxy()

    undf = xp%vspace%get_last_dof_owned()
    !$omp parallel do schedule(static), default(none), shared(xp, undf, a),  private(i)
    do i = 1,undf
      xp%data(i) = xp%data(i)**a
    end do
    !$omp end parallel do

    call xp%set_dirty()
  end subroutine invoke_raise_field

!> invoke_raise_field: Raise a field to a real exponent
!> x => x^a
  subroutine invoke_real_raise_field(x, a)
    implicit none
    type( field_type ), intent(inout)  :: x
    real(r_def),        intent(in)     :: a
    type( field_proxy_type)            :: xp
    integer(kind=i_def)                :: i,undf

    xp = x%get_proxy()

    undf = xp%vspace%get_last_dof_owned()
    !$omp parallel do schedule(static), default(none), shared(xp, undf, a),  private(i)
    do i = 1,undf
      xp%data(i) = xp%data(i)**a
    end do
    !$omp end parallel do

    call xp%set_dirty()
  end subroutine invoke_real_raise_field

!-------------------------------------------------------------------------------
! Implmented in #965, kernel requires stencil support. Note that the w2_field is
! required to obtain the W2 stencil_cross which is used in determining
! orientation of cells in the halo
  subroutine invoke_mpi_calc_cell_orientation(w2_field,cell_orientation)

    use mesh_mod, only: mesh_type
    use stencil_dofmap_mod,               only : stencil_dofmap_type, &
                                                 STENCIL_CROSS
    use calc_cell_orientation_kernel_mod, only : calc_cell_orientation_code

    implicit none

    type(field_type), intent(in)      :: w2_field
    type(field_type), intent(inout)   :: cell_orientation

    integer                 :: cell, nlayers
    integer                 :: ndf_w3
    integer                 :: undf_w3
    integer                 :: ndf_w2
    integer, pointer        :: map_w3(:) => null()

    type(field_proxy_type) :: cell_orientation_proxy
    type(field_proxy_type) :: w2_field_proxy

    type(mesh_type), pointer           :: mesh => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_w2 => null()
    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size

    mesh => w2_field%get_mesh()

    cell_orientation_proxy = cell_orientation%get_proxy()
    w2_field_proxy = w2_field%get_proxy()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = w2_field_proxy%vspace%get_ndf( )

    ! Obtain the stencil for core cells only
    cross_stencil_w2 => w2_field_proxy%vspace%get_stencil_dofmap(             &
                                        STENCIL_CROSS, mesh%get_halo_depth())
    cross_stencil_w2_map => cross_stencil_w2%get_whole_dofmap()

    cross_stencil_w3 => cell_orientation_proxy%vspace%get_stencil_dofmap(     &
                                        STENCIL_CROSS, mesh%get_halo_depth())
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only
      map_w3 => cell_orientation_proxy%vspace%get_cell_dofmap(cell)

      call calc_cell_orientation_code(  nlayers,                              &
                                        cell_orientation_proxy%data,          &
                                        undf_w3,                              &
                                        ndf_w3,                               &
                                        map_w3,                               &
                                        ndf_w2,                               &
                                        cross_stencil_w3_size,                &
                                        cross_stencil_w2_map(:,:,cell),       &
                                        cross_stencil_w3_map(:,:,cell) )
    end do

  end subroutine invoke_mpi_calc_cell_orientation

!=============================================================================!
! #999 Requires psyclone support for enforce_operator_bc_code,
! being implemented as part of issue #22
! #1001 will implement algortihm layer calls to the kernel
  subroutine invoke_enforce_operator_bc_kernel_type(op)
    use enforce_operator_bc_kernel_mod, only: enforce_operator_bc_code
    use mesh_mod, only: mesh_type
    type(operator_type), intent(inout) :: op
    integer, pointer                   :: boundary_dofs(:,:) => null()
    integer                            :: cell, ncell_3d
    integer                            :: ndf1, ndf2
    type(mesh_type), pointer           :: mesh => null()
    integer                            :: nlayers
    type(operator_proxy_type)          :: op_proxy
    !
    ! Initialise operator proxies
    !
    op_proxy = op%get_proxy()
    !
    ! Initialise number of layers
    !
    nlayers = op_proxy%fs_to%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => op%get_mesh()
    !
    ! Get size of operator array
    !
    ncell_3d = op_proxy%ncell_3d
    ndf1 = op_proxy%fs_to%get_ndf()
    ndf2 = op_proxy%fs_from%get_ndf()
    !
    ! Pull out boundary array
    !
    boundary_dofs => op_proxy%fs_to%get_boundary_dofs()
    !
    ! Call kernels and communication routines
    !
    do cell=1,mesh%get_last_halo_cell(1)
      !
      call enforce_operator_bc_code(cell, &
                                    nlayers, &
                                    op_proxy%local_stencil, &
                                    ncell_3d, &
                                    ndf1, &
                                    ndf2, &
                                    boundary_dofs)
    end do 
  end subroutine invoke_enforce_operator_bc_kernel_type

!-------------------------------------------------------------------------------   
!> invoke_inc_aX_times_Y: x = a * x * y ; a-scalar, x,y-vector     
  subroutine invoke_inc_aX_times_Y(scalar, field1, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    real(kind=r_def),   intent(in)     :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inc_aX_times_Y:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, undf, scalar),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = scalar * field1_proxy%data(i) * field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_inc_aX_times_Y
!-------------------------------------------------------------------------------   
!> invoke_inc_X_times_Y: x = x * y ; x,y-vector     
  subroutine invoke_inc_X_times_Y(field1, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)  :: field1
    type( field_type ), intent(in)     :: field2
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inc_X_times_Y:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy, field2_proxy, undf),  private(i)
    do i = 1,undf
      field1_proxy%data(i) = field1_proxy%data(i) * field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => field1%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field1_proxy%set_dirty()
      else
        call field1_proxy%set_clean(dplp)
      end if
    end do
  end subroutine invoke_inc_X_times_Y

  !-------------------------------------------------------------------------------
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data. 
  !> PSyclone support is required - documented in #942
  subroutine invoke_get_height_kernel(height, chi )
    use get_height_kernel_mod, only: get_height_code
    use mesh_mod,              only: mesh_type
    implicit none
    
    type(field_type), intent(inout)      :: height
    type(field_type), intent(in)         :: chi(3) 

    type(field_proxy_type) :: x_p, chi_p(3)
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x
    integer                 :: undf_chi, undf_x
    integer                 :: dim_chi
    integer                 :: df_x, df_chi

    integer, pointer        :: map_chi(:) => null()
    integer, pointer        :: map_x(:) => null()
    real(kind=r_def), pointer :: nodes_x(:,:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    integer :: i
    type(mesh_type), pointer :: mesh => null()

    x_p   = height%get_proxy()
    do i = 1,3
      chi_p(i) = chi(i)%get_proxy()
    end do

    nlayers = x_p%vspace%get_nlayers()

    ndf_x  = x_p%vspace%get_ndf( )
    undf_x = x_p%vspace%get_undf()
    nodes_x => x_p%vspace%get_nodes()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    dim_chi = chi_p(1)%vspace%get_dim_space( )

    ! Evaluate the basis function
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))
    do df_x = 1, ndf_x
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
      end do
    end do

    if (chi_p(1)%is_dirty(depth=1)) then
       call chi_p(1)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(2)%is_dirty(depth=1)) then
       call chi_p(2)%halo_exchange(depth=1)
    end if
      !
    if (chi_p(3)%is_dirty(depth=1)) then
       call chi_p(3)%halo_exchange(depth=1)
    end if
    mesh => height%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
       map_x   => x_p%vspace%get_cell_dofmap( cell )
       map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
       call get_height_code(nlayers,                    &
                            x_p%data,                   &
                            chi_p(1)%data,              &
                            chi_p(2)%data,              &
                            chi_p(3)%data,              &
                            ndf_x, undf_x, map_x,       &
                            ndf_chi, undf_chi, map_chi, &
                            basis_chi                   &
                            )
    end do

    call x_p%set_dirty()

    deallocate(basis_chi)
  end subroutine invoke_get_height_kernel
  
  !-------------------------------------------------------------------------------
  !> In #937, the evaluator was removed in the PSy-lite layer. In #938, evaluator
  !> (not quadrature) will be removed and the functionality will be implemented via
  !> PSY using kernal meta data. 
  !> PSyclone support is required - documented in #942
  subroutine invoke_mpi_detj_at_w2(detj_at_w2, chi)
    use calc_detj_at_w2_kernel_mod, only : calc_detj_at_w2_code
    use mesh_mod, only                   : mesh_type

    implicit none

    type(field_type), intent(inout) :: detj_at_w2
    type(field_type), intent(in)    :: chi(3)

    type(field_proxy_type)          :: detj_at_w2_p, chi_p(3)
    type(mesh_type)                 :: mesh

    integer                         :: cell, nlayers
    integer                         :: ndf_chi, ndf_w2
    integer                         :: undf_chi, undf_w2
    integer                         :: dim_u, diff_dim_chi
    integer, pointer                :: map_chi(:) => null()
    integer, pointer                :: map_w2(:)  => null()


    real(kind=r_def), allocatable   :: nodal_basis_u(:,:,:)
    real(kind=r_def), allocatable   :: diff_basis_chi(:,:,:)
    real(kind=r_def), pointer       :: nodes(:,:) => null()
    integer                         :: ii, df_chi, df_w2

    do ii = 1,3
      chi_p(ii)  = chi(ii)%get_proxy()
    end do
    detj_at_w2_p = detj_at_w2%get_proxy()

    nlayers = detj_at_w2_p%vspace%get_nlayers()

    ndf_w2  = detj_at_w2_p%vspace%get_ndf( )
    undf_w2 = detj_at_w2_p%vspace%get_undf()
    dim_u = detj_at_w2_p%vspace%get_dim_space()
    allocate(nodal_basis_u(dim_u, ndf_w2, ndf_w2))

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )
    allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf_w2))

    nodes => detj_at_w2_p%vspace%get_nodes( )

    do df_w2 = 1, ndf_w2
      do df_chi = 1, ndf_chi
        diff_basis_chi(:,df_chi,df_w2) = chi_p(1)%vspace%call_function(DIFF_BASIS,df_chi,nodes(:,df_w2))
      end do
    end do

    do df_w2 = 1, ndf_w2
      do df_chi = 1, ndf_w2
        nodal_basis_u(:,df_chi,df_w2) =  detj_at_w2_p%vspace%call_function(BASIS,df_chi,nodes(:,df_w2))
      end do
    end do

    mesh = detj_at_w2%get_mesh()

    do cell = 1, mesh%get_last_halo_cell(1)
      map_w2  => detj_at_w2_p%vspace%get_cell_dofmap( cell )
      map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )

      call calc_detj_at_w2_code( nlayers,       &
         detj_at_w2_p%data,                     &
         chi_p(1)%data,                         &
         chi_p(2)%data,                         &
         chi_p(3)%data,                         &
         ndf_w2, undf_w2, map_w2,               &
         ndf_chi, undf_chi, map_chi,            &
         diff_basis_chi                         )

    end do

    deallocate(nodal_basis_u)
    deallocate(diff_basis_chi)

  end subroutine invoke_mpi_detj_at_w2

end module psykal_lite_mod
