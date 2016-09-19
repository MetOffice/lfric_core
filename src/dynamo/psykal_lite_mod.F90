!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_mod

  use field_mod,          only : field_type, field_proxy_type 
  use scalar_mod,         only : scalar_type
  use operator_mod,       only : operator_type, operator_proxy_type
  use quadrature_mod,     only : quadrature_type
  use constants_mod,      only : r_def, i_def, cache_block
  use mesh_mod,           only : mesh_type
  use evaluator_xyoz_mod, only : quadrature_xyoz_type, &
                                 quadrature_xyoz_proxy_type
  use evaluator_xoyoz_mod
  use evaluator_xyz_mod
  use evaluate_function_mod, only : BASIS, DIFF_BASIS

  implicit none
  public

contains

    !------------------------------------------------------------------------------------------
    !> Computation of nodal basis function for coordinates chi not currently supported PSyClone,
    !> will be introduced in modification of quadrature strategy, see ticket #723.
    SUBROUTINE invoke_compute_geopotential_kernel_type(geopotential, chi)
      USE compute_geopotential_kernel_mod, ONLY: compute_geopotential_code
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(inout) :: geopotential, chi(3)
      INTEGER, pointer :: map_w0(:), map_chi(:) => null()
      INTEGER cell
      REAL(KIND=r_def), allocatable :: basis_chi(:,:,:)
      INTEGER ndf_w0, undf_w0, ndf_chi, undf_chi, dim_chi
      TYPE(mesh_type) mesh
      INTEGER nlayers
      TYPE(field_proxy_type) geopotential_proxy, chi_proxy(3)

      real(kind=r_def), pointer :: nodes(:,:) => null()

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
      mesh = geopotential%get_mesh()
      !
      ! Initialise sizes and allocate any basis arrays for w0
      !
      ndf_w0 = geopotential_proxy%vspace%get_ndf()
      undf_w0 = geopotential_proxy%vspace%get_undf()
      !
      ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
      undf_chi  = chi_proxy(1)%vspace%get_undf( )
      dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

      allocate( basis_chi(dim_chi, ndf_chi, ndf_w0) )

      nodes => geopotential_proxy%vspace%get_nodes( )

      call chi_proxy(1)%vspace%compute_nodal_basis_function( basis_chi, ndf_chi, ndf_w0, nodes)

      DO cell=1, geopotential_proxy%vspace%get_ncell()

        map_w0 => geopotential_proxy%vspace%get_cell_dofmap(cell)
        map_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)

        CALL compute_geopotential_code(nlayers, geopotential_proxy%data, chi_proxy(1)%data, chi_proxy(2)%data, chi_proxy(3)%data, ndf_w0, undf_w0, map_w0, ndf_chi, undf_chi, map_chi, basis_chi)
      END DO
    END SUBROUTINE invoke_compute_geopotential_kernel_type

  !------------------------------------------------------------------------------
  !> invoke_initial_theta_kernel: invoke the potential temperature initialization for a generic space
  subroutine invoke_initial_theta_kernel( theta, chi )

    use initial_theta_kernel_mod, only : initial_theta_code

    implicit none

    type( field_type ), intent( inout ) :: theta
    type( field_type ), intent( in ) :: chi(3)

    integer          :: cell
    integer          :: ndf_wtheta, undf_wtheta, &
        ndf_chi, undf_chi, dim_chi

    integer, pointer        :: map_wtheta(:) => null()
    integer, pointer        :: map_chi(:)    => null()

    type( field_proxy_type )        :: theta_proxy
    type( field_proxy_type )        :: chi_proxy(3)

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    real(kind=r_def), pointer :: nodes(:,:) => null()

    theta_proxy  = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ndf_wtheta  = theta_proxy%vspace%get_ndf( )
    undf_wtheta  = theta_proxy%vspace%get_undf( )

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi  = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

    allocate( basis_chi(dim_chi, ndf_chi, ndf_wtheta) )

    nodes => theta_proxy%vspace%get_nodes( )

    call chi_proxy(1)%vspace%compute_nodal_basis_function( basis_chi, ndf_chi, ndf_wtheta, nodes)

    do cell = 1, theta_proxy%vspace%get_ncell()

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
  end subroutine invoke_initial_theta_kernel

  !-------------------------------------------------------------------------------
  !> Invoke_rtheta_bd_kernel: Invoke the boundary part of the RHS of the theta equation
  subroutine invoke_rtheta_bd_kernel( r_theta_bd, theta, f, rho, qr )

    use rtheta_bd_kernel_mod,   only : rtheta_bd_code
    use mesh_mod, only : mesh_type ! Work around for intel_v15 failues on the Cray

    implicit none

    type (mesh_type)                 :: mesh
    type( field_type ), intent( in ) :: theta, f, rho
    type( field_type ), intent( inout ) :: r_theta_bd
    type( quadrature_type), intent( in ) :: qr

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_wtheta, ndf_w3
    integer                 :: undf_w2, undf_wtheta, undf_w3
    integer                 :: dim_w2, dim_wtheta, dim_w3
    integer, pointer        :: map_w2(:), map_w3(:), map_wtheta(:) => null()
    integer, pointer        :: map_wtheta_W(:), map_w3_W(:), map_W2_W(:) => null()
    integer, pointer        :: map_wtheta_S(:), map_w3_S(:), map_W2_S(:) => null()
    integer, pointer        :: map_wtheta_E(:), map_w3_E(:), map_W2_E(:) => null()
    integer, pointer        :: map_wtheta_N(:), map_w3_N(:), map_W2_N(:) => null()
!    integer, pointer        :: orientation_w2(:) => null()

    type( field_proxy_type )        :: r_theta_bd_proxy, f_proxy, theta_proxy, rho_proxy

    real(kind=r_def), allocatable  :: basis_wtheta_f_W(:,:,:,:),   &
      basis_wtheta_f_S(:,:,:,:),     &
      basis_wtheta_f_E(:,:,:,:),     &
      basis_wtheta_f_N(:,:,:,:),     &
      basis_wtheta_face(:,:,:,:,:),  &
      basis_w2_f_W(:,:,:,:),         &
      basis_w2_f_S(:,:,:,:),         &
      basis_w2_f_E(:,:,:,:),         &
      basis_w2_f_N(:,:,:,:),         &
      basis_w2_face(:,:,:,:,:),      &
      basis_w3_f_W(:,:,:,:),         &
      basis_w3_f_S(:,:,:,:),         &
      basis_w3_f_E(:,:,:,:),         &
      basis_w3_f_N(:,:,:,:),         &
      basis_w3_face(:,:,:,:,:),      &
      xp_f_W(:,:), xp_f_S(:,:), xp_f_E(:,:), xp_f_N(:,:)



    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    mesh = theta%get_mesh()

    r_theta_bd_proxy = r_theta_bd%get_proxy()
    theta_proxy  = theta%get_proxy()
    f_proxy  = f%get_proxy()
    rho_proxy = rho%get_proxy()

    nlayers = theta_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    xp=>qr%get_xqp_h()
    zp=>qr%get_xqp_v()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal quad points in x and y
    nqp_h_1d = int(nqp_h**(1./2.))

    ndf_w2      = f_proxy%vspace%get_ndf( )
    dim_w2      = f_proxy%vspace%get_dim_space( )
    undf_w2     = f_proxy%vspace%get_undf()
    allocate(basis_w2_f_W(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_S(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_E(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_N(dim_w2,ndf_w2 ,nqp_h_1d,nqp_v))
    allocate(basis_w2_face(4,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_w3      = rho_proxy%vspace%get_ndf( )
    dim_w3      = rho_proxy%vspace%get_dim_space( )
    undf_w3     = rho_proxy%vspace%get_undf()
    allocate(basis_w3_f_W(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_S(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_E(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_N(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_face(4,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_f_W(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_S(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_E(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_N(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_face(4,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    allocate(xp_f_W(nqp_h_1d, 2))
    allocate(xp_f_S(nqp_h_1d, 2))
    allocate(xp_f_E(nqp_h_1d, 2))
    allocate(xp_f_N(nqp_h_1d, 2))

    ! Quadrature points on horizontal faces

    xp_f_W(:, :) = xp(1:nqp_h_1d, :)
    xp_f_W(:, 1) = 0.0_r_def

    xp_f_S(:, :) = xp( 1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f_S(:, 2) = 0.0_r_def

    xp_f_E(:, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f_E(:, 1) = 1.0_r_def

    xp_f_N(:, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f_N(:, 2) = 1.0_r_def

    ! Filling up the basis vector with value of the basis functions at the horizontal faces quadrature points

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_W, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_W, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_S, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_S, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_E, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_E, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_N, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_N, zp)

    basis_wtheta_face(1,:,:,:,:) = basis_wtheta_f_W
    basis_wtheta_face(2,:,:,:,:) = basis_wtheta_f_S
    basis_wtheta_face(3,:,:,:,:) = basis_wtheta_f_E
    basis_wtheta_face(4,:,:,:,:) = basis_wtheta_f_N

    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_W, ndf_w3, nqp_h_1d, nqp_v, xp_f_W, zp)

    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_S, ndf_w3, nqp_h_1d, nqp_v, xp_f_S, zp)

    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_E, ndf_w3, nqp_h_1d, nqp_v, xp_f_E, zp)

    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_N, ndf_w3, nqp_h_1d, nqp_v, xp_f_N, zp)

    basis_w3_face(1,:,:,:,:) = basis_w3_f_W
    basis_w3_face(2,:,:,:,:) = basis_w3_f_S
    basis_w3_face(3,:,:,:,:) = basis_w3_f_E
    basis_w3_face(4,:,:,:,:) = basis_w3_f_N

    call f_proxy%vspace%compute_basis_function( &
      basis_w2_f_W, ndf_w2, nqp_h_1d, nqp_v, xp_f_W, zp)

    call f_proxy%vspace%compute_basis_function( &
      basis_w2_f_S, ndf_w2, nqp_h_1d, nqp_v, xp_f_S, zp)

    call f_proxy%vspace%compute_basis_function( &
      basis_w2_f_E, ndf_w2, nqp_h_1d, nqp_v, xp_f_E, zp)

    call f_proxy%vspace%compute_basis_function( &
      basis_w2_f_N, ndf_w2, nqp_h_1d, nqp_v, xp_f_N, zp)

    basis_w2_face(1,:,:,:,:) = basis_w2_f_W
    basis_w2_face(2,:,:,:,:) = basis_w2_f_S
    basis_w2_face(3,:,:,:,:) = basis_w2_f_E
    basis_w2_face(4,:,:,:,:) = basis_w2_f_N

    do cell = 1, theta_proxy%vspace%get_ncell()

      map_w2 => f_proxy%vspace%get_cell_dofmap( cell )
      map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )

      ! Getting dofmaps on neighbouring cells for w3 and wtheta

      map_w2_W => f_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(1, cell) )
      map_w2_S => f_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(2, cell) )
      map_w2_E => f_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(3, cell) )
      map_w2_N => f_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(4, cell) )

      map_w3_W => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(1, cell) )
      map_w3_S => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(2, cell) )
      map_w3_E => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(3, cell) )
      map_w3_N => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(4, cell) )

      map_wtheta_W => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(1, cell) )
      map_wtheta_S => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(2, cell) )
      map_wtheta_E => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(3, cell) )
      map_wtheta_N => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(4, cell) )

      call rtheta_bd_code( nlayers,                      &
        ndf_w2, undf_w2,                                 &
        map_w2,                                          &
        map_w2_W, map_w2_S,                              &
        map_w2_E, map_w2_N,                              &
        r_theta_bd_proxy%data,                           &
        ndf_w3, undf_w3,                                 &
        map_w3,                                          &
        map_w3_W, map_w3_S,                              &
        map_w3_E, map_w3_N,                              &
        rho_proxy%data,                                  &
        theta_proxy%data,                                &
        f_proxy%data,                                    &
        ndf_wtheta, undf_wtheta, map_wtheta,             &
        map_wtheta_W, map_wtheta_S,                      &
        map_wtheta_E, map_wtheta_N,                      &
        nqp_v, nqp_h_1d, wv,                             &
        basis_w2_face, basis_w3_face, basis_wtheta_face  &
        )

    end do

    deallocate(basis_w2_f_W, basis_w2_f_S, basis_w2_f_E, basis_w2_f_N, basis_w2_face, &
      basis_w3_f_W, basis_w3_f_S, basis_w3_f_E, basis_w3_f_N, basis_w3_face,          &
      basis_wtheta_f_W, basis_wtheta_f_S, basis_wtheta_f_E, basis_wtheta_f_N, basis_wtheta_face, &
      xp_f_W, xp_f_S, xp_f_E, xp_f_N )

  end subroutine invoke_rtheta_bd_kernel

  !-------------------------------------------------------------------------------
  !> Invoke_ru_bd_kernel: Invoke the boundary part of the RHS of the momentum equation
  subroutine invoke_ru_bd_kernel( r_u_bd, rho, theta, qr )

    use ru_bd_kernel_mod,   only : ru_bd_code
    use mesh_mod,           only : mesh_type ! Work around for intel_v15 failues on the Cray

    implicit none

    type( field_type ), intent( inout ) :: r_u_bd
    type( field_type ), intent( in ) :: rho, theta
    type( quadrature_type), intent( in ) :: qr
    type( mesh_type )                    :: mesh

    integer                 :: cell, nlayers, nqp_h, nqp_v, nqp_h_1d
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w3(:), map_w2(:), map_wtheta(:) => null()
    integer, pointer        :: map_w3_W(:), map_wtheta_W(:) => null()
    integer, pointer        :: map_w3_S(:), map_wtheta_S(:) => null()
    integer, pointer        :: map_w3_E(:), map_wtheta_E(:) => null()
    integer, pointer        :: map_w3_N(:), map_wtheta_N(:) => null()

    type( field_proxy_type )        :: r_u_bd_proxy, rho_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_w2_f_W(:,:,:,:), &
      basis_w2_f_S(:,:,:,:), &
      basis_w2_f_E(:,:,:,:), &
      basis_w2_f_N(:,:,:,:), &
      basis_w3_f_W(:,:,:,:), &
      basis_wtheta_f_W(:,:,:,:), &
      basis_w3_f_S(:,:,:,:), &
      basis_wtheta_f_S(:,:,:,:), &
      basis_w3_f_E(:,:,:,:), &
      basis_wtheta_f_E(:,:,:,:), &
      basis_w3_f_N(:,:,:,:), &
      basis_wtheta_f_N(:,:,:,:), &
      basis_w2_face(:,:,:,:,:), &
      basis_w3_face(:,:,:,:,:), &
      basis_wtheta_face(:,:,:,:,:)

    real(kind=r_def), pointer :: xp(:,:), xp_f_W(:,:), xp_f_S(:,:), xp_f_E(:,:), xp_f_N(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    mesh = rho%get_mesh()

    r_u_bd_proxy = r_u_bd%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()

    nlayers = rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ! Assumes same number of horizontal qp in x and y
    nqp_h_1d = int(nqp_h**(1./2.))  ! use sqrt

    allocate(xp_f_W(nqp_h_1d, 2))
    allocate(xp_f_S(nqp_h_1d, 2))
    allocate(xp_f_E(nqp_h_1d, 2))
    allocate(xp_f_N(nqp_h_1d, 2))


    ndf_w3  = rho_proxy%vspace%get_ndf( )
    dim_w3  = rho_proxy%vspace%get_dim_space( )
    undf_w3 = rho_proxy%vspace%get_undf()
    allocate(basis_w3_f_W(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_S(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_E(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_f_N(dim_w3,ndf_w3,nqp_h_1d,nqp_v))
    allocate(basis_w3_face(4,dim_w3,ndf_w3,nqp_h_1d,nqp_v))

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_f_W(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_S(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_E(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_f_N(dim_w2,ndf_w2,nqp_h_1d,nqp_v))
    allocate(basis_w2_face(4,dim_w2,ndf_w2,nqp_h_1d,nqp_v))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_f_W(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_S(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_E(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_f_N(dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))
    allocate(basis_wtheta_face(4,dim_wtheta,ndf_wtheta,nqp_h_1d,nqp_v))

    allocate(xp_f_W(nqp_h_1d, 2))
    allocate(xp_f_S(nqp_h_1d, 2))
    allocate(xp_f_E(nqp_h_1d, 2))
    allocate(xp_f_N(nqp_h_1d, 2))

    ! Quadrature points on horizontal faces

    xp_f_W(:, :) = xp(1:nqp_h_1d, :)
    xp_f_W(:, 1) = 0.0_r_def

    xp_f_S(:, :) = xp( 1:nqp_h - nqp_h_1d + 1:nqp_h_1d, :)
    xp_f_S(:, 2) = 0.0_r_def

    xp_f_E(:, :) = xp(nqp_h - nqp_h_1d + 1:nqp_h, :)
    xp_f_E(:, 1) = 1.0_r_def

    xp_f_N(:, :) = xp(nqp_h_1d:nqp_h:nqp_h_1d, :)
    xp_f_N(:, 2) = 1.0_r_def

    ! Filling up the basis vector with value of the basis functions at the horizontal faces quadrature points

    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_W, ndf_w3, nqp_h_1d, nqp_v, xp_f_W, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_W, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_W, zp)

    call r_u_bd_proxy%vspace%compute_basis_function( &
      basis_w2_f_W, ndf_w2, nqp_h_1d, nqp_v, xp_f_W, zp)


    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_S, ndf_w3, nqp_h_1d, nqp_v, xp_f_S, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_S, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_S, zp)

    call r_u_bd_proxy%vspace%compute_basis_function( &
      basis_w2_f_S, ndf_w2, nqp_h_1d, nqp_v, xp_f_S, zp)


    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_E, ndf_w3, nqp_h_1d, nqp_v, xp_f_E, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_E, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_E, zp)

    call r_u_bd_proxy%vspace%compute_basis_function( &
      basis_w2_f_E, ndf_w2, nqp_h_1d, nqp_v, xp_f_E, zp)


    call rho_proxy%vspace%compute_basis_function( &
      basis_w3_f_N, ndf_w3, nqp_h_1d, nqp_v, xp_f_N, zp)

    call theta_proxy%vspace%compute_basis_function( &
      basis_wtheta_f_N, ndf_wtheta, nqp_h_1d, nqp_v, xp_f_N, zp)

    call r_u_bd_proxy%vspace%compute_basis_function( &
      basis_w2_f_N, ndf_w2, nqp_h_1d, nqp_v, xp_f_N, zp)

    ! Grouping face basis vectors in a indexed higher-dimensional vector

    basis_w3_face(1,:,:,:,:) = basis_w3_f_W
    basis_w3_face(2,:,:,:,:) = basis_w3_f_S
    basis_w3_face(3,:,:,:,:) = basis_w3_f_E
    basis_w3_face(4,:,:,:,:) = basis_w3_f_N

    basis_wtheta_face(1,:,:,:,:) = basis_wtheta_f_W
    basis_wtheta_face(2,:,:,:,:) = basis_wtheta_f_S
    basis_wtheta_face(3,:,:,:,:) = basis_wtheta_f_E
    basis_wtheta_face(4,:,:,:,:) = basis_wtheta_f_N

    basis_w2_face(1,:,:,:,:) = basis_w2_f_W
    basis_w2_face(2,:,:,:,:) = basis_w2_f_S
    basis_w2_face(3,:,:,:,:) = basis_w2_f_E
    basis_w2_face(4,:,:,:,:) = basis_w2_f_N

    do cell = 1, r_u_bd_proxy%vspace%get_ncell()

      map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
      map_w2 => r_u_bd_proxy%vspace%get_cell_dofmap( cell )
      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )

      ! Getting dofmaps on neighbouring cells for w3 and wtheta

      map_w3_W => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(1, cell) )
      map_wtheta_W => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(1, cell) )

      map_w3_S => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(2,cell) )
      map_wtheta_S => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(2, cell) )

      map_w3_E => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(3, cell) )
      map_wtheta_E => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(3, cell) )

      map_w3_N => rho_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(4, cell) )
      map_wtheta_N => theta_proxy%vspace%get_cell_dofmap( mesh%get_cell_next(4, cell) )

      call ru_bd_code( nlayers,           &
        ndf_w2, undf_w2,                  &
        map_w2,                           &
        r_u_bd_proxy%data,                &
        ndf_w3, undf_w3,                  &
        map_w3,                           &
        map_w3_W, map_w3_S,               &
        map_w3_E, map_w3_N,               &
        rho_proxy%data,                   &
        theta_proxy%data,                 &
        ndf_wtheta, undf_wtheta,          &
        map_wtheta,                       &
        map_wtheta_W, map_wtheta_S,       &
        map_wtheta_E, map_wtheta_N,       &
        nqp_v, nqp_h_1d, wv,              &
        basis_w2_face,                    &
        basis_w3_face, basis_wtheta_face  &
        )

    end do

    deallocate(basis_w3_f_W, basis_w3_f_S, basis_w3_f_E, basis_w3_f_N, basis_w3_face, &
      basis_w2_f_W, basis_w2_f_S, basis_w2_f_E, basis_w2_f_N, basis_w2_face, &
      basis_wtheta_f_W, basis_wtheta_f_S, basis_wtheta_f_E, basis_wtheta_f_N, basis_wtheta_face)

  end subroutine invoke_ru_bd_kernel

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
    type(mesh_type),  pointer          :: mesh => null()

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
!> invoke_multiply_field: Compute y = a*x for scalar a and fields y and x
  subroutine invoke_multiply_field(a, x, y)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in)    :: x
    type( field_type ), intent(inout) :: y
    real(kind=r_def),   intent(in)    :: a
    type( field_proxy_type)           :: x_proxy, y_proxy
    integer(kind=i_def)               :: i,undf
    integer(kind=i_def)               :: depth, dplp
    type(mesh_type), pointer          :: mesh => null()

    x_proxy = x%get_proxy()
    y_proxy = y%get_proxy()

    undf = x_proxy%vspace%get_undf()
    if(undf /= y_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:multiply_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    end if
    !$omp parallel do schedule(static), default(none), shared(y_proxy,x_proxy, undf, a),  private(i)
    do i = 1,undf
      y_proxy%data(i) = a*x_proxy%data(i)
    end do
    !$omp end parallel do

    mesh => y%get_mesh()
    depth = mesh%get_halo_depth()
    
    do dplp = 1, depth
      if( x_proxy%is_dirty(depth=dplp)  ) then
        call y_proxy%set_dirty()
      else
        call y_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_multiply_field

!-------------------------------------------------------------------------------   
!> Non pointwise Kernels


!------------------------------------------------------------------------------- 
!> invoke_sample_flux_kernel: Retrieve values from flux kernel  
  subroutine invoke_sample_flux_kernel(flux, u, multiplicity, q)
    use sample_flux_kernel_mod, only: sample_flux_code
    implicit none

    type(field_type), intent(in)    :: u, multiplicity, q
    type(field_type), intent(inout) :: flux

    type(field_proxy_type)          :: u_p, m_p, q_p, flux_p
    
    integer                 :: cell, nlayers
    integer                 :: ndf_f, ndf_q
    integer                 :: undf_f, undf_q
    integer                 :: dim_q
    integer, pointer        :: map_f(:), map_q(:) => null()

    real(kind=r_def), allocatable  :: basis_q(:,:,:)

    real(kind=r_def), pointer :: nodes(:,:) => null()

    u_p    = u%get_proxy()
    q_p    = q%get_proxy()
    m_p    = multiplicity%get_proxy()
    flux_p = flux%get_proxy()

    nlayers = flux_p%vspace%get_nlayers()

    ndf_f  = flux_p%vspace%get_ndf( )
    undf_f = flux_p%vspace%get_undf()

    ndf_q  = q_p%vspace%get_ndf( )
    dim_q  = q_p%vspace%get_dim_space( )
    undf_q = q_p%vspace%get_undf()
    allocate(basis_q(dim_q, ndf_q, ndf_f))

    nodes => flux_p%vspace%get_nodes( )
    call q_p%vspace%compute_nodal_basis_function(basis_q, ndf_q, ndf_f, nodes)

    if(flux_p%is_dirty(depth=1) ) then
      call flux_p%halo_exchange(depth=1)
    end if
    if(u_p%is_dirty(depth=1) ) then
      call u_p%halo_exchange(depth=1)
    end if
    if(m_p%is_dirty(depth=1) ) then
      call m_p%halo_exchange(depth=1)
    end if
    if(q_p%is_dirty(depth=1) ) then
      call q_p%halo_exchange(depth=1)
    end if

    do cell = 1, flux_p%vspace%get_ncell()
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

  subroutine invoke_nodal_coordinates_kernel(nodal_coords, chi)
    use nodal_coordinates_kernel_mod, only: nodal_coordinates_code
    implicit none
    
    type(field_type), intent(inout) :: nodal_coords(3)
    type(field_type), intent(in)    :: chi(3)
 
    type(field_proxy_type) :: x_p(3), chi_p(3)
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf_x
    integer                 :: undf_chi, undf_x
    integer                 :: dim_chi
    integer, pointer        :: map_chi(:), map_x(:) => null()

    real(kind=r_def), allocatable  :: basis_chi(:,:,:)
    real(kind=r_def), pointer :: nodes(:,:) => null()
    integer :: i

    do i = 1,3
      x_p(i)   = nodal_coords(i)%get_proxy()
      chi_p(i) = chi(i)%get_proxy()
    end do

    nlayers = x_p(1)%vspace%get_nlayers()

    ndf_x  = x_p(1)%vspace%get_ndf( )
    undf_x = x_p(1)%vspace%get_undf()

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()

    dim_chi = chi_p(1)%vspace%get_dim_space( )
    allocate(basis_chi(dim_chi, ndf_chi, ndf_x))

    nodes => x_p(1)%vspace%get_nodes( )
    call chi_p(1)%vspace%compute_nodal_basis_function(basis_chi, ndf_chi, ndf_x, nodes)    

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

    do cell = 1, x_p(1)%vspace%get_ncell()
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

  subroutine invoke_convert_hcurl_field(phys_field, comp_field, chi)
    use convert_hcurl_field_kernel_mod, only: convert_hcurl_field_code
    implicit none
    
    type(field_type), intent(inout) :: phys_field(3)
    type(field_type), intent(in)    :: chi(3), comp_field
 
    type(field_proxy_type) :: phys_p(3), chi_p(3), comp_p
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf
    integer                 :: undf_chi, undf
    integer                 :: diff_dim_chi, dim
    integer, pointer        :: map_chi(:), map(:) => null()

    real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:), basis(:,:,:)
    real(kind=r_def), pointer :: nodes(:,:) => null()
    integer(kind=i_def) :: i

    do i = 1,3
      phys_p(i) = phys_field(i)%get_proxy()
      chi_p(i)  = chi(i)%get_proxy()
    end do
    comp_p = comp_field%get_proxy()

    nlayers = comp_p%vspace%get_nlayers()

    ndf  = comp_p%vspace%get_ndf( )
    undf = comp_p%vspace%get_undf()
    dim  = comp_p%vspace%get_dim_space( )
    allocate(basis(dim, ndf, ndf) )

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )
    allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf))

    nodes => comp_p%vspace%get_nodes( )
    call chi_p(1)%vspace%compute_nodal_diff_basis_function(diff_basis_chi, ndf_chi, ndf, nodes)    
    call comp_p%vspace%compute_nodal_basis_function(basis, ndf, ndf, nodes) 

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

    do cell = 1, comp_p%vspace%get_ncell()
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
                                     ndf, undf, map, &
                                     ndf_chi, undf_chi, map_chi, &
                                     basis, &
                                     diff_basis_chi &
                                    )
    end do

    call phys_p(1)%set_dirty()
    call phys_p(2)%set_dirty()
    call phys_p(3)%set_dirty()

    deallocate(diff_basis_chi, basis)
  end subroutine invoke_convert_hcurl_field

!-------------------------------------------------------------------------------   

  subroutine invoke_convert_hdiv_field(phys_field, comp_field, chi)
    use convert_hdiv_field_kernel_mod, only: convert_hdiv_field_code
    implicit none
    
    type(field_type), intent(inout) :: phys_field(3)
    type(field_type), intent(in)    :: chi(3), comp_field
 
    type(field_proxy_type) :: phys_p(3), chi_p(3), comp_p
   
    integer                 :: cell, nlayers
    integer                 :: ndf_chi, ndf
    integer                 :: undf_chi, undf
    integer                 :: diff_dim_chi, dim
    integer, pointer        :: map_chi(:), map(:) => null()

    real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:), basis(:,:,:)
    real(kind=r_def), pointer :: nodes(:,:) => null()
    integer :: i

    do i = 1,3
      phys_p(i) = phys_field(i)%get_proxy()
      chi_p(i)  = chi(i)%get_proxy()
    end do
    comp_p = comp_field%get_proxy()

    nlayers = comp_p%vspace%get_nlayers()

    ndf  = comp_p%vspace%get_ndf( )
    undf = comp_p%vspace%get_undf()
    dim  = comp_p%vspace%get_dim_space( )
    allocate(basis(dim, ndf, ndf) )

    ndf_chi  = chi_p(1)%vspace%get_ndf( )
    undf_chi = chi_p(1)%vspace%get_undf()
    diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )
    allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf))

    nodes => comp_p%vspace%get_nodes( )
    call chi_p(1)%vspace%compute_nodal_diff_basis_function(diff_basis_chi, ndf_chi, ndf, nodes)    
    call comp_p%vspace%compute_nodal_basis_function(basis, ndf, ndf, nodes) 

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

    do cell = 1, comp_p%vspace%get_ncell()
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
                                    ndf, undf, map, &
                                    ndf_chi, undf_chi, map_chi, &
                                    basis, &
                                    diff_basis_chi &
                                   )
    end do

    call phys_p(1)%set_dirty()
    call phys_p(2)%set_dirty()
    call phys_p(3)%set_dirty()

    deallocate(diff_basis_chi, basis)
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

  implicit none

  type(field_type), intent(inout) :: level
  type(field_proxy_type) :: l_p
  integer :: cell, ncell, ndf, undf
  real(kind=r_def), pointer :: nodes(:,:) => null()
  integer, pointer :: map(:) => null()

  l_p = level%get_proxy()
  undf = l_p%vspace%get_undf()
  ndf  = l_p%vspace%get_ndf()
  nodes => l_p%vspace%get_nodes( )
 
  ncell = l_p%vspace%get_ncell()

  do cell = 1,ncell
    map => l_p%vspace%get_cell_dofmap(cell)
    call compute_dof_level_code(l_p%vspace%get_nlayers(),                 &
                                l_p%data,                                 &
                                ndf,                                      &
                                undf,                                     &
                                map,                                      &
                                nodes                                     &
                               )
  end do 

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
subroutine invoke_subgrid_coeffs(a0,a1,a2,rho,direction,rho_stencil_extent)

    use flux_direction_mod,        only: x_direction, y_direction
    use stencil_dofmap_mod,        only: stencil_dofmap_type, &
                                         STENCIL_1DX,         &
                                         STENCIL_1DY
    use subgrid_coeffs_kernel_mod, only: subgrid_coeffs_code
    use subgrid_config_mod,        only: rho_approximation

    implicit none

    type( field_type ), intent( inout ) :: a0
    type( field_type ), intent( inout ) :: a1
    type( field_type ), intent( inout ) :: a2
    type( field_type ), intent( in )    :: rho
    integer, intent(in)                 :: direction
    integer, intent(in)                 :: rho_stencil_extent

    type( field_proxy_type )            :: rho_proxy
    type( field_proxy_type )            :: a0_proxy
    type( field_proxy_type )            :: a1_proxy
    type( field_proxy_type )            :: a2_proxy

    type(stencil_dofmap_type), pointer  :: map => null()
    integer, pointer                    :: stencil_map(:,:) => null()
    integer                             :: rho_stencil_size
    integer                 :: cell
    integer                 :: nlayers
    integer                 :: ndf_w3
    integer                 :: undf_w3

    a0_proxy   = a0%get_proxy()
    a1_proxy   = a1%get_proxy()
    a2_proxy   = a2%get_proxy()
    rho_proxy  = rho%get_proxy()

    undf_w3 = rho_proxy%vspace%get_undf()
    ndf_w3  = rho_proxy%vspace%get_ndf()
    nlayers = rho_proxy%vspace%get_nlayers()

    ! Note stencil grid types are of the form:
    !                                   |5|
    !                                   |3|
    ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
    !                                   |2|
    !                                   |4|
    if (direction .EQ. x_direction) then
      map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_stencil_extent)
    elseif (direction .EQ. y_direction) then
      map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_stencil_extent)
    end if
    rho_stencil_size = map%get_size()

    do cell = 1, rho_proxy%vspace%get_ncell()

      stencil_map => map%get_dofmap(cell)

      call subgrid_coeffs_code( nlayers,                                  &
                                rho_approximation,                        &
                                undf_w3,                                  &
                                rho_proxy%data,                           &
                                ndf_w3,                                   &
                                rho_stencil_size,                         &
                                stencil_map,                              &
                                a0_proxy%data,                            &
                                a1_proxy%data,                            &
                                a2_proxy%data                             &
                                )

    end do

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
  if (direction .EQ. x_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,stencil_extent)
  elseif (direction .EQ. y_direction) then
    map => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,stencil_extent)
  end if
  stencil_size = map%get_size()

  do cell = 1, rho_proxy%vspace%get_ncell()
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

end subroutine invoke_conservative_fluxes

!-------------------------------------------------------------------------------
subroutine invoke_calc_departure_wind(u_departure_wind, u_piola, chi)
  use calc_departure_wind_kernel_mod, only: calc_departure_wind_code
  implicit none

  type(field_type), intent(inout) :: u_departure_wind
  type(field_type), intent(in)    :: chi(3), u_piola

  type(field_proxy_type) :: u_departure_wind_p, chi_p(3), u_piola_p

  integer                 :: cell, nlayers
  integer                 :: ndf_chi, ndf
  integer                 :: undf_chi, undf
  integer                 :: dim_u, diff_dim_chi
  integer, pointer        :: map_chi(:), map(:) => null()

  real(kind=r_def), allocatable  :: nodal_basis_u(:,:,:)
  real(kind=r_def), allocatable  :: diff_basis_chi(:,:,:)
  real(kind=r_def), pointer :: nodes(:,:) => null()
  integer :: ii

  do ii = 1,3
    chi_p(ii)  = chi(ii)%get_proxy()
  end do
  u_piola_p = u_piola%get_proxy()
  u_departure_wind_p = u_departure_wind%get_proxy()

  nlayers = u_piola_p%vspace%get_nlayers()

  ndf  = u_piola_p%vspace%get_ndf( )
  undf = u_piola_p%vspace%get_undf()
  dim_u = u_piola_p%vspace%get_dim_space()
  allocate(nodal_basis_u(dim_u, ndf, ndf))

  ndf_chi  = chi_p(1)%vspace%get_ndf( )
  undf_chi = chi_p(1)%vspace%get_undf()
  diff_dim_chi = chi_p(1)%vspace%get_dim_space_diff( )
  allocate(diff_basis_chi(diff_dim_chi, ndf_chi, ndf))

  nodes => u_piola_p%vspace%get_nodes( )
  call chi_p(1)%vspace%compute_nodal_diff_basis_function(diff_basis_chi, ndf_chi, ndf, nodes)
  call u_piola_p%vspace%compute_nodal_basis_function(nodal_basis_u, ndf, ndf, nodes)

  do cell = 1, u_piola_p%vspace%get_ncell()
     map     => u_piola_p%vspace%get_cell_dofmap( cell )
     map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )
     call calc_departure_wind_code( nlayers,                                  &
                                    u_departure_wind_p%data,                  &
                                    u_piola_p%data,                           &
                                    chi_p(1)%data,                            &
                                    chi_p(2)%data,                            &
                                    chi_p(3)%data,                            &
                                    ndf, undf, map, nodal_basis_u,            &
                                    ndf_chi, undf_chi, map_chi,               &
                                    diff_basis_chi                            &
                                     )
  end do
  deallocate(nodal_basis_u)
  deallocate(diff_basis_chi)
end subroutine invoke_calc_departure_wind
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> invoke_calc_deppts: Invoke the calculation of departure points in 1D
subroutine invoke_calc_deppts(u_n,u_np1,dep_pts,direction,dep_pt_method)

  use calc_departure_point_kernel_mod,  only : calc_departure_point_code
  use stencil_dofmap_mod,               only : stencil_dofmap_type, &
                                               STENCIL_1DX, &
                                               STENCIL_1DY
  use flux_direction_mod,               only : x_direction, y_direction
  use subgrid_config_mod,               only : transport_stencil_extent

  implicit none

  type( field_type ), intent( in )    :: u_n
  type( field_type ), intent( in )    :: u_np1
  type( field_type ), intent( inout ) :: dep_pts
  integer, intent(in)                 :: direction
  integer, intent(in)                 :: dep_pt_method

  type( field_proxy_type )        :: u_n_proxy
  type( field_proxy_type )        :: u_np1_proxy
  type( field_proxy_type )        :: dep_pts_proxy
  type(stencil_dofmap_type), pointer  :: map=>null()

  integer, pointer        :: stencil_map_w2(:,:) => null()
  integer                 :: transport_stencil_size

  integer                 :: cell
  integer                 :: nlayers
  integer                 :: ndf_w2
  integer                 :: undf_w2

  u_n_proxy    = u_n%get_proxy()
  u_np1_proxy  = u_np1%get_proxy()
  dep_pts_proxy = dep_pts%get_proxy()

  ndf_w2  = u_n_proxy%vspace%get_ndf()
  undf_w2 = u_n_proxy%vspace%get_undf()

  if (direction .EQ. x_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,transport_stencil_extent)
  elseif (direction .EQ. y_direction) then
    map => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,transport_stencil_extent)
  endif
  transport_stencil_size = map%get_size()

  nlayers = u_n_proxy%vspace%get_nlayers()

  do cell=1,u_n_proxy%vspace%get_ncell()

    stencil_map_w2 => map%get_dofmap(cell)

    call calc_departure_point_code( nlayers,                      &
                                    dep_pts_proxy%data,           &
                                    transport_stencil_size,       &
                                    undf_w2,                      &
                                    ndf_w2,                       &
                                    stencil_map_w2,               &
                                    u_n_proxy%data,               &
                                    u_np1_proxy%data,             &
                                    direction,                    &
                                    dep_pt_method )

  end do

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
    type(mesh_type)                    :: mesh

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

    mesh = field1%get_mesh()
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
    type(mesh_type)                    :: mesh

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

    mesh = y%get_mesh()
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

!> invoke_divide_field: Divide the values of field1 by field2
!> c = a/b
  subroutine invoke_divide_field_data(field1, field2)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    use mesh_mod,only : mesh_type
    implicit none
    type( field_type ), intent(inout)   ::field1
    type( field_type ), intent(in)     :: field2
    type( field_proxy_type)            :: field1_proxy,field2_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type)                    :: mesh

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

    mesh = field1%get_mesh()
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
    type(mesh_type)                    :: mesh

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

    mesh = field_res%get_mesh()
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
    type(mesh_type)                    :: mesh

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

    mesh = field1%get_mesh()
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
    type(mesh_type)                    :: mesh

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

    mesh = field1%get_mesh()
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

    !$omp parallel do schedule(static), default(none), shared(field_proxy, undf, scalar),  private(i)
    do i = 1,undf
       field_proxy%data(i) = scalar*field_proxy%data(i)
    end do
    !$omp end parallel do
   
  end subroutine invoke_scale_field_data


!-------------------------------------------------------------------------------   

!-------------------------------------------------------------------------------
subroutine invoke_set_boundary_kernel(field, bc)
  use set_boundary_kernel_mod, only: set_boundary_code
  use mesh_mod,                only: mesh_type ! Work around for intel_v15 failues on the Cray

  implicit none

  type(field_type), intent(inout) :: field
  real(kind=r_def), intent(in)    :: bc

  integer, pointer :: boundary_dofs(:,:) => null()
  integer, pointer :: map(:) => null()
  integer :: cell,  nlayers
  integer :: ndf, undf
  type(mesh_type)        :: mesh
  type(field_proxy_type) :: f_proxy

  f_proxy = field%get_proxy()
  nlayers = f_proxy%vspace%get_nlayers()
  mesh = field%get_mesh()
  undf = f_proxy%vspace%get_undf()
  ndf  = f_proxy%vspace%get_ndf()
      
  if (f_proxy%is_dirty(depth=1)) then
    call f_proxy%halo_exchange(depth=1)
  end if 

  boundary_dofs => f_proxy%vspace%get_boundary_dofs()

  do cell=1,mesh%get_last_halo_cell(1)
         
    map => f_proxy%vspace%get_cell_dofmap(cell)
         
    CALL set_boundary_code(nlayers, &
                           f_proxy%data, &
                           ndf, &
                           undf, &
                           map, &
                           boundary_dofs, &
                           bc)
         
  end do 
       
  call f_proxy%set_dirty()
       
end subroutine invoke_set_boundary_kernel

!> invoke_times_field: times the values of field1 by field2 and put result in
!>field_res
!> c = a/b
  subroutine invoke_times_field(field1,field2,field_res)
    use log_mod,  only : log_event, LOG_LEVEL_ERROR
    use mesh_mod, only: mesh_type ! Work around for intel_v15 failues on the Cray
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type)                    :: mesh

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:times_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:times_field:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)*field2_proxy%data(i)
    end do
    !$omp end parallel do

    mesh = field_res%get_mesh()
    depth = mesh%get_halo_depth()
    
    do dplp = 1, depth
      if( field1_proxy%is_dirty(depth=dplp) .or. &
          field2_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_times_field

!-------------------------------------------------------------------------------   

!-------------------------------------------------------------------------------   
subroutine invoke_sample_poly_flux( flux, wind, density, stencil_extent) 

  use sample_poly_flux_kernel_mod, only: sample_poly_flux_code
  use stencil_dofmap_mod,          only: stencil_dofmap_type, STENCIL_CROSS
  use log_mod, only : log_event, LOG_LEVEL_ERROR

  implicit none

  type(field_type), intent(inout)   :: flux
  type(field_type), intent(in)      :: wind
  type(field_type), intent(in)      :: density
  integer, intent(in)               :: stencil_extent

  type( field_proxy_type )  :: flux_proxy, wind_proxy, density_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()

  integer, pointer :: map_w2(:,:)        => null()
  integer, pointer :: stencil_map(:,:,:) => null()

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: cell
  integer :: nlayers
  integer :: stencil_size 
  
  real(kind=r_def), allocatable :: basis_w2(:,:,:)
  integer :: dim_w2
  real(kind=r_def), pointer :: nodes(:,:) => null()

  flux_proxy    = flux%get_proxy()
  wind_proxy    = wind%get_proxy()
  density_proxy = density%get_proxy()

  ndf_w3  = density_proxy%vspace%get_ndf()
  undf_w3 = density_proxy%vspace%get_undf()

  ndf_w2  = flux_proxy%vspace%get_ndf()
  undf_w2 = flux_proxy%vspace%get_undf()
  dim_w2  = flux_proxy%vspace%get_dim_space()
  allocate (basis_w2(dim_w2, ndf_w2, ndf_w2))
  nodes => flux_proxy%vspace%get_nodes( )
  call flux_proxy%vspace%compute_nodal_basis_function(basis_w2, ndf_w2, ndf_w2, nodes)

  nlayers = flux_proxy%vspace%get_nlayers()

  stencil => density_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                     stencil_extent)
  stencil_size = stencil%get_size()
  stencil_map => stencil%get_whole_dofmap()
  map_w2 => flux_proxy%vspace%get_whole_dofmap()

  if(wind_proxy%is_dirty(depth=1) ) then
    call wind_proxy%halo_exchange(depth=1)
  end if
  ! This needs extending for larger haloes
  if(density_proxy%is_dirty(depth=1) ) then
    call density_proxy%halo_exchange(depth=1)
  end if

  do cell = 1, flux_proxy%vspace%get_ncell()

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
subroutine invoke_sample_poly_adv( adv, wind, tracer, stencil_extent) 

  use sample_poly_adv_kernel_mod, only: sample_poly_adv_code
  use stencil_dofmap_mod,          only: stencil_dofmap_type, STENCIL_CROSS
  use log_mod, only : log_event, LOG_LEVEL_ERROR

  implicit none

  type(field_type), intent(inout)   :: adv
  type(field_type), intent(in)      :: wind
  type(field_type), intent(in)      :: tracer
  integer, intent(in)               :: stencil_extent

  type( field_proxy_type )  :: adv_proxy, wind_proxy, tracer_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()

  integer, pointer :: map_w2(:,:)        => null()
  integer, pointer :: stencil_map(:,:,:) => null()

  integer :: undf_wt, ndf_wt
  integer :: undf_w2, ndf_w2
  integer :: cell
  integer :: nlayers
  integer :: stencil_size
  
  real(kind=r_def), allocatable :: basis_w2(:,:,:)
  integer :: dim_w2
  real(kind=r_def), pointer :: nodes(:,:) => null()

  adv_proxy    = adv%get_proxy()
  wind_proxy   = wind%get_proxy()
  tracer_proxy = tracer%get_proxy()

  ndf_wt  = tracer_proxy%vspace%get_ndf()
  undf_wt = tracer_proxy%vspace%get_undf()

  ndf_w2  = wind_proxy%vspace%get_ndf()
  undf_w2 = wind_proxy%vspace%get_undf()
  dim_w2  = wind_proxy%vspace%get_dim_space()
  allocate (basis_w2(dim_w2, ndf_w2, ndf_wt))
  nodes => tracer_proxy%vspace%get_nodes( )
  call wind_proxy%vspace%compute_nodal_basis_function(basis_w2, ndf_w2, ndf_wt, nodes)

  nlayers = adv_proxy%vspace%get_nlayers()

  stencil => tracer_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                    stencil_extent)
  stencil_size = stencil%get_size()
  stencil_map => stencil%get_whole_dofmap()
  map_w2 => wind_proxy%vspace%get_whole_dofmap()

  if(wind_proxy%is_dirty(depth=1) ) then
    call wind_proxy%halo_exchange(depth=1)
  end if
  ! This needs extending for larger haloes
  if(tracer_proxy%is_dirty(depth=1) ) then
    call tracer_proxy%halo_exchange(depth=1)
  end if

  do cell = 1, adv_proxy%vspace%get_ncell()

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
! Needs differential nodal basis for chi space which is not currently
! supported by psyclone. Dynamo tickets #723, #734 address this as new types of
! quadrature that will then need to be supported by psyclone
    subroutine invoke_compute_tri_precon_kernel(tri_precon, theta, rho, chi)
      use compute_tri_precon_kernel_mod, only: compute_tri_precon_code
      use mesh_mod, only: mesh_type
      type(field_type), intent(inout) :: tri_precon(3)
      type(field_type), intent(in) :: theta, rho, chi(3)
      integer, pointer :: map_w3(:) => null(), map_w0(:) => null(), map_any_space_1_chi(:) => null()
      integer :: cell
      integer :: ndf_w3, undf_w3, ndf_w0, undf_w0, ndf_any_space_1_chi, undf_any_space_1_chi
      type(mesh_type), pointer :: mesh => null()
      integer :: nlayers
      type(field_proxy_type) :: tri_precon_proxy(3), theta_proxy, rho_proxy, chi_proxy(3)
      real(kind=r_def), pointer :: nodes(:,:) => null()
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
      !
      ! Initialise sizes and allocate any basis arrays for w0
      !
      ndf_w0 = theta_proxy%vspace%get_ndf()
      undf_w0 = theta_proxy%vspace%get_undf()
      !
      ! Initialise sizes and allocate any basis arrays for any_space_1_chi
      !
      ndf_any_space_1_chi = chi_proxy(1)%vspace%get_ndf()
      undf_any_space_1_chi = chi_proxy(1)%vspace%get_undf()

      ! Compute nodal basis functions
      diff_dim_chi  = chi_proxy(1)%vspace%get_dim_space_diff( )
      allocate( diff_basis_chi(diff_dim_chi, ndf_any_space_1_chi, ndf_w3) )
      nodes => rho_proxy%vspace%get_nodes( )
      call chi_proxy(1)%vspace%compute_nodal_diff_basis_function( diff_basis_chi, ndf_any_space_1_chi, ndf_w3, nodes)
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
        map_w0 => theta_proxy%vspace%get_cell_dofmap(cell)
        map_any_space_1_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)
        !
        CALL compute_tri_precon_code(nlayers, &
                                     tri_precon_proxy(1)%data, &
                                     tri_precon_proxy(2)%data, &
                                     tri_precon_proxy(3)%data, &
                                     theta_proxy%data, &
                                     rho_proxy%data, &
                                     chi_proxy(1)%data, &
                                     chi_proxy(2)%data, &
                                     chi_proxy(3)%data, &
                                     ndf_w3, &
                                     undf_w3, &
                                     map_w3, &
                                     ndf_w0, &
                                     undf_w0, &
                                     map_w0, &
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

end module psykal_lite_mod
