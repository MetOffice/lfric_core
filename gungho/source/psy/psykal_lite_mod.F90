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

  use field_mod,                    only : field_type, field_proxy_type 
  use scalar_mod,                   only : scalar_type
  use operator_mod,                 only : operator_type, operator_proxy_type
  use constants_mod,                only : r_def, i_def, cache_block
  use mesh_mod,                     only : mesh_type
  use function_space_mod,           only : BASIS, DIFF_BASIS

  ! The following modules are not currently implemented as part of the main
  ! code but they do have unit tests so need to be declared here so they are
  ! built ready for unit testing.

  use quadrature_xoyoz_mod
  use quadrature_xyz_mod
  use quadrature_xyoz_mod, only : quadrature_xyoz_type, &
                                  quadrature_xyoz_proxy_type
  use quadrature_face_mod, only : quadrature_face_type, &
                                  quadrature_face_proxy_type

!  use quadrature_xyoz_mod,          only : quadrature_xyoz_type
!  use quadrature_rule_gaussian_mod, only : quadrature_rule_gaussian_type


  implicit none
  public

contains



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
  subroutine invoke_initial_mr_kernel( theta, rho, mr, chi )

    use initial_mr_kernel_mod, only : initial_mr_code
    use mr_indices_mod,        only : imr_v, imr_c, imr_r, imr_nc, imr_nr, nummr
    use mesh_mod,              only : mesh_type
    implicit none

    type( field_type ), intent( inout ) :: mr(nummr)
    type( field_type ), intent( in ) :: theta, rho
    type( field_type ), intent( in ) :: chi(3)

    integer          :: cell
    integer          :: ndf_wtheta, undf_wtheta, &
                        ndf_w3, undf_w3, &
                        ndf_chi, undf_chi
    integer, pointer :: map_wtheta(:) => null()
    integer, pointer :: map_w3(:)    => null()
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
    rho_proxy = rho%get_proxy()

    ndf_w3  = rho_proxy%vspace%get_ndf( )
    undf_w3 = rho_proxy%vspace%get_undf( )

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
      map_w3     => rho_proxy%vspace%get_cell_dofmap( cell )
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
        ndf_w3,                             &
        undf_w3,                            &
        map_w3,                             &
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

    use reference_element_mod,  only : reference_element_type
    use rtheta_bd_kernel_mod,   only : rtheta_bd_code
    use mesh_mod, only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,     only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( in )     :: theta, u
    type( field_type ), intent( inout )  :: r_theta_bd
    type(quadrature_face_type), intent( in ) :: qr

    type (mesh_type),              pointer :: mesh => null()
    class(reference_element_type), pointer :: reference_element =>null()
    type(stencil_dofmap_type),     pointer :: cross_stencil_w2 => null()
    type(stencil_dofmap_type),     pointer :: cross_stencil_wtheta => null()

    type(quadrature_face_proxy_type) :: qr_proxy

    integer                 :: cell, nlayers, nqp
    integer                 :: ndf_w2, ndf_wtheta
    integer                 :: undf_w2, undf_wtheta
    integer                 :: dim_w2, dim_wtheta
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer                 :: cross_stencil_w2_size

    integer, pointer        :: cross_stencil_wtheta_map(:,:,:) => null()
    integer                 :: cross_stencil_wtheta_size

    type( field_proxy_type )        :: r_theta_bd_proxy, u_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_wtheta_face(:,:,:,:),  &
                                      basis_w2_face(:,:,:,:)

    real(kind=r_def), pointer :: wqp(:,:) => null()

    real(r_def), allocatable :: normal_to_face(:,:)
    real(r_def), allocatable :: out_face_normal(:,:)

    integer(i_def) :: nfaces_h

    mesh              => theta%get_mesh()
    reference_element => mesh%get_reference_element()

    nfaces_h = reference_element%get_number_horizontal_faces()

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

    qr_proxy = qr%get_quadrature_proxy()
    nqp=qr_proxy%np_xyz
    wqp=>qr_proxy%weights_xyz

    ndf_w2      = u_proxy%vspace%get_ndf( )
    dim_w2      = u_proxy%vspace%get_dim_space( )
    undf_w2     = u_proxy%vspace%get_undf()
    allocate(basis_w2_face(dim_w2,ndf_w2,nqp,nfaces_h))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp,nfaces_h))

    call qr%compute_function(BASIS, u_proxy%vspace,     dim_w2,     ndf_w2,     basis_w2_face)
    call qr%compute_function(BASIS, theta_proxy%vspace, dim_wtheta, ndf_wtheta, basis_wtheta_face)

    if(r_theta_bd_proxy%is_dirty(depth=2) ) call r_theta_bd_proxy%halo_exchange(depth=2)
    if(theta_proxy%is_dirty(depth=2) ) call theta_proxy%halo_exchange(depth=2)
    if(u_proxy%is_dirty(depth=2) ) call u_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()
    call reference_element%get_normals_to_faces( normal_to_face )
    call reference_element%get_normals_to_out_faces( out_face_normal )

    do cell = 1, mesh%get_last_halo_cell(1)

      call rtheta_bd_code( nlayers,                            &
                           ndf_w2, undf_w2,                    &
                           cross_stencil_w2_map(:,:,cell),     &
                           cross_stencil_w2_size,              &
                           ndf_wtheta, undf_wtheta,            &
                           cross_stencil_wtheta_map(:,:,cell), &
                           cross_stencil_wtheta_size,          &
                           r_theta_bd_proxy%data,              &
                           theta_proxy%data,                   &
                           u_proxy%data,                       &
                           nqp, wqp,                           &
                           basis_w2_face, basis_wtheta_face,   &
                           adjacent_face(:,cell),              &
                           normal_to_face,                     &
                           out_face_normal )

    end do
    call r_theta_bd_proxy%set_dirty()

    deallocate( basis_wtheta_face, basis_w2_face )

  end subroutine invoke_rtheta_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_ru_bd_kernel: Invoke the boundary part of the RHS of the momentum equation
  subroutine invoke_ru_bd_kernel( r_u_bd, exner, theta, qr )

    use reference_element_mod,  only : reference_element_type
    use ru_bd_kernel_mod,       only : ru_bd_code
    use mesh_mod,               only : mesh_type ! Work around for intel_v15 failues on the Cray
    use stencil_dofmap_mod,     only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( in )     :: exner, theta
    type( field_type ), intent( inout )  :: r_u_bd
    type( quadrature_face_type ), intent( in ) :: qr

    type(mesh_type),               pointer           :: mesh    => null()
    class(reference_element_type), pointer :: reference_element => null()

    type(quadrature_face_proxy_type) :: qr_proxy
    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()

    integer                 :: cell, nlayers, nqp
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w2(:,:) => null()
    integer, pointer        :: map_wtheta(:,:) => null()
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size

    type( field_proxy_type )        :: r_u_bd_proxy, exner_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:), &
                                      basis_w3_face(:,:,:,:), &
                                      basis_wtheta_face(:,:,:,:)

    real(kind=r_def), pointer :: wqp(:,:) => null()  

    integer(i_def)           :: nfaces_h
    real(r_def), allocatable :: out_face_normal(:,:)

    mesh => exner%get_mesh()
    reference_element => mesh%get_reference_element()
    nfaces_h = reference_element%get_number_horizontal_faces()

    r_u_bd_proxy = r_u_bd%get_proxy()
    exner_proxy  = exner%get_proxy()
    theta_proxy  = theta%get_proxy()

    cross_stencil_w3 => exner_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    nlayers = exner_proxy%vspace%get_nlayers()
    qr_proxy = qr%get_quadrature_proxy()
    nqp = qr_proxy%np_xyz
    wqp => qr_proxy%weights_xyz

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_face(dim_w2,ndf_w2,nqp,nfaces_h))

    ndf_w3  = exner_proxy%vspace%get_ndf( )
    dim_w3  = exner_proxy%vspace%get_dim_space( )
    undf_w3 = exner_proxy%vspace%get_undf()
    allocate(basis_w3_face(dim_w3,ndf_w3,nqp,nfaces_h))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp,nfaces_h))

    call qr%compute_function(BASIS, r_u_bd_proxy%vspace, dim_w2,     ndf_w2,     basis_w2_face)
    call qr%compute_function(BASIS, exner_proxy%vspace,  dim_w3,     ndf_w3,     basis_w3_face)
    call qr%compute_function(BASIS, theta_proxy%vspace,  dim_wtheta, ndf_wtheta, basis_wtheta_face)

    if(theta_proxy%is_dirty(depth=2) ) call theta_proxy%halo_exchange(depth=2)
    if(exner_proxy%is_dirty(depth=2) ) call exner_proxy%halo_exchange(depth=2)
    if(r_u_bd_proxy%is_dirty(depth=2) ) call r_u_bd_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()
    call reference_element%get_normals_to_out_faces( out_face_normal )

    map_w2     => r_u_bd_proxy%vspace%get_whole_dofmap( )
    map_wtheta => theta_proxy%vspace%get_whole_dofmap( )
    do cell = 1, mesh%get_last_halo_cell(1)


      call ru_bd_code( nlayers,                           &
                       ndf_w2, undf_w2,                   &
                       map_w2(:,cell),                    &
                       ndf_w3, undf_w3,                   &
                       cross_stencil_w3_map(:,:,cell),    &
                       cross_stencil_w3_size,             &
                       ndf_wtheta, undf_wtheta,           &
                       map_wtheta(:,cell),                &
                       r_u_bd_proxy%data,                 &
                       exner_proxy%data,                  &
                       theta_proxy%data,                  &
                       nqp, wqp,                          &
                       basis_w2_face,                     &
                       basis_w3_face, basis_wtheta_face,  &
                       adjacent_face(:,cell),             &
                       out_face_normal )

    end do
    call r_u_bd_proxy%set_dirty()

    deallocate( basis_wtheta_face, basis_w3_face, basis_w2_face )

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
    use reference_element_mod,        only : reference_element_type

    implicit none

    type( field_type ), intent( in )     :: exner, theta
    type( field_type ), intent( inout )  :: r_u_bd
    type( quadrature_face_type ), intent( in ) :: qr

    type(mesh_type),               pointer :: mesh              => null()
    class(reference_element_type), pointer :: reference_element => null()

    type(quadrature_face_proxy_type) :: qr_proxy
    type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()

    integer                 :: cell, nlayers, nqp
    integer                 :: ndf_w2, ndf_w3, ndf_wtheta
    integer                 :: undf_w2, undf_w3, undf_wtheta
    integer                 :: dim_w2, dim_w3, dim_wtheta
    integer, pointer        :: map_w2(:) => null()
    integer, pointer        :: map_wtheta(:) => null()
    integer, pointer        :: adjacent_face(:,:) => null()

    integer, pointer        :: cross_stencil_w3_map(:,:,:) => null()
    integer                 :: cross_stencil_w3_size

    type( field_proxy_type )        :: r_u_bd_proxy, exner_proxy, theta_proxy

    real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:)
    real(kind=r_def), allocatable  :: basis_w3_face(:,:,:,:)
    real(kind=r_def), allocatable  :: basis_wtheta_face(:,:,:,:)

    real(kind=r_def), pointer :: wqp(:,:) => null()

    integer(i_def)           :: nfaces_h
    real(r_def), allocatable :: out_face_normal(:,:)

    mesh => exner%get_mesh()
    reference_element => mesh%get_reference_element()
    nfaces_h = reference_element%get_number_horizontal_faces()

    r_u_bd_proxy = r_u_bd%get_proxy()
    exner_proxy  = exner%get_proxy()
    theta_proxy  = theta%get_proxy()

    cross_stencil_w3 => exner_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
    cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
    cross_stencil_w3_size = cross_stencil_w3%get_size()

    nlayers = exner_proxy%vspace%get_nlayers()
    qr_proxy = qr%get_quadrature_proxy()
    nqp =  qr_proxy%np_xyz
    wqp => qr_proxy%weights_xyz

    ndf_w2      = r_u_bd_proxy%vspace%get_ndf( )
    dim_w2      = r_u_bd_proxy%vspace%get_dim_space( )
    undf_w2     = r_u_bd_proxy%vspace%get_undf()
    allocate(basis_w2_face(nfaces_h,dim_w2,ndf_w2,nqp))

    ndf_w3  = exner_proxy%vspace%get_ndf( )
    dim_w3  = exner_proxy%vspace%get_dim_space( )
    undf_w3 = exner_proxy%vspace%get_undf()
    allocate(basis_w3_face(nfaces_h,dim_w3,ndf_w3,nqp))

    ndf_wtheta      = theta_proxy%vspace%get_ndf( )
    dim_wtheta      = theta_proxy%vspace%get_dim_space( )
    undf_wtheta     = theta_proxy%vspace%get_undf()
    allocate(basis_wtheta_face(nfaces_h,dim_wtheta,ndf_wtheta,nqp))

    call qr%compute_function(BASIS, r_u_bd_proxy%vspace, dim_w2,     ndf_w2,     basis_w2_face)
    call qr%compute_function(BASIS, exner_proxy%vspace,  dim_w3,     ndf_w3,     basis_w3_face)
    call qr%compute_function(BASIS, theta_proxy%vspace,  dim_wtheta, ndf_wtheta, basis_wtheta_face)

    if( theta_proxy%is_dirty(depth=1) ) call theta_proxy%halo_exchange(depth=1)


    if(exner_proxy%is_dirty(depth=2) ) call exner_proxy%halo_exchange(depth=2)
    if(r_u_bd_proxy%is_dirty(depth=2) ) call r_u_bd_proxy%halo_exchange(depth=2)

    adjacent_face => mesh%get_adjacent_face()
    call reference_element%get_normals_to_out_faces( out_face_normal )

    do cell = 1,mesh%get_last_halo_cell(1)

      map_w2 => r_u_bd_proxy%vspace%get_cell_dofmap( cell )
      map_wtheta => theta_proxy%vspace%get_cell_dofmap( cell )

      call exner_gradient_bd_code( nlayers,                          &
                                   ndf_w2, undf_w2,                  &
                                   map_w2,                           &
                                   ndf_w3, undf_w3,                  &
                                   cross_stencil_w3_map(:,:,cell),   &
                                   cross_stencil_w3_size,            &
                                   ndf_wtheta, undf_wtheta,          &
                                   map_wtheta,                       &
                                   r_u_bd_proxy%data,                &
                                   exner_proxy%data,                 &
                                   theta_proxy%data,                 &
                                   nqp, wqp,                         &
                                   basis_w2_face,                    &
                                   basis_w3_face, basis_wtheta_face, &
                                   adjacent_face(:,cell),            &
                                   out_face_normal )

    end do
    call r_u_bd_proxy%set_dirty()

    deallocate(basis_w3_face, basis_w2_face, basis_wtheta_face)


  end subroutine invoke_exner_gradient_bd_kernel

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_div_bd_kernel: Invoke the boundary part of the divergence of the lhs Helmholtz
  subroutine invoke_weighted_div_bd_kernel_type(div_star, theta, s, qr)

      use weighted_div_bd_kernel_mod, only : weighted_div_bd_code
      use mesh_mod,                   only : mesh_type
      use reference_element_mod,      only : reference_element_type
      use stencil_dofmap_mod,         only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type(field_type),           intent(in)    :: theta
      type(operator_type),        intent(inout) :: div_star
      type(quadrature_face_type), intent(in)    :: qr
      real(r_def),                intent(in)    :: s

      type(mesh_type),               pointer :: mesh              => null()
      class(reference_element_type), pointer :: reference_element => null()

      type(quadrature_face_proxy_type) :: qr_proxy
      integer :: cell, nlayers, nqp
      integer :: ndf_w2, ndf_w3, ndf_wtheta, undf_wtheta
      integer :: dim_w2, dim_w3, dim_wtheta

      integer,                   pointer :: cross_stencil_wtheta_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()
      integer                            :: cross_stencil_wtheta_size

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:), &
                                        basis_w3_face(:,:,:,:), &
                                        basis_wtheta_face(:,:,:,:)

      real(kind=r_def), pointer :: wqp(:,:) => null()

      type(operator_proxy_type) ::div_star_proxy
      type(field_proxy_type)    ::theta_proxy

      integer(i_def) :: nfaces_h

      integer(i_def), pointer     :: adjacent_face(:,:) => null()
      real(r_def),    allocatable :: out_face_normal(:,:)

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
      mesh              => div_star%get_mesh()
      reference_element => mesh%get_reference_element()

      nfaces_h = reference_element%get_number_horizontal_faces()

      !
      ! Initialise qr values
      !
      qr_proxy = qr%get_quadrature_proxy()
      nqp =  qr_proxy%np_xyz
      wqp => qr_proxy%weights_xyz

      ndf_w2      = div_star_proxy%fs_to%get_ndf( )
      dim_w2      = div_star_proxy%fs_to%get_dim_space( )
      allocate(basis_w2_face(dim_w2,ndf_w2,nqp,nfaces_h))

      ndf_w3  = div_star_proxy%fs_from%get_ndf( )
      dim_w3  = div_star_proxy%fs_from%get_dim_space( )
      allocate(basis_w3_face(dim_w3,ndf_w3,nqp,nfaces_h))

      ndf_wtheta      = theta_proxy%vspace%get_ndf( )
      dim_wtheta      = theta_proxy%vspace%get_dim_space( )
      undf_wtheta     = theta_proxy%vspace%get_undf()
      allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp,nfaces_h))

      call qr%compute_function(BASIS, div_star_proxy%fs_to,   dim_w2,     ndf_w2,     basis_w2_face)
      call qr%compute_function(BASIS, div_star_proxy%fs_from, dim_w3,     ndf_w3,     basis_w3_face)
      call qr%compute_function(BASIS, theta_proxy%vspace,     dim_wtheta, ndf_wtheta, basis_wtheta_face)

      ! Stencil maps
      cross_stencil_wtheta => theta_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
      cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=2)) call theta_proxy%halo_exchange(depth=2)
      !

      adjacent_face => mesh%get_adjacent_face()
      call reference_element%get_normals_to_out_faces( out_face_normal )

      do cell=1,mesh%get_last_halo_cell(1)

        call weighted_div_bd_code(                                           &
                            cell,                                            &
                            nlayers,                                         &
                            div_star_proxy%ncell_3d,                         &
                            div_star_proxy%local_stencil,                    &
                            theta_proxy%data,                                &
                            s,                                               &
                            ndf_w2, ndf_w3,                                  &
                            ndf_wtheta, undf_wtheta,                         &
                            cross_stencil_wtheta_map(:,:,cell),              &
                            cross_stencil_wtheta_size,                       &
                            nqp, wqp,                                        &
                            basis_w2_face,                                   &
                            basis_w3_face, basis_wtheta_face,                &
                            adjacent_face(:,cell),                           &
                            reference_element%get_number_horizontal_faces(), &
                            out_face_normal )
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate( basis_wtheta_face, basis_w3_face, basis_w2_face )
      !

      nullify( mesh,                     &
               cross_stencil_wtheta_map, &
               cross_stencil_wtheta,     &
               wqp, adjacent_face )

    end subroutine invoke_weighted_div_bd_kernel_type

  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_proj_2theta_bd_kernel: Invoke the boundary part of the projection from wtheta to w2 operator of the lhs Helmholtz
  subroutine invoke_weighted_proj_2theta_bd_kernel_type(p2theta, exner, s, qr)

      use weighted_proj_2theta_bd_kernel_mod, only : weighted_proj_2theta_bd_code
      use mesh_mod,                           only : mesh_type
      use reference_element_mod,              only : reference_element_type
      use stencil_dofmap_mod,                 only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type(field_type), intent(in)         :: exner
      real(r_def), intent(in)              :: s
      type(operator_type), intent(inout)   :: p2theta
      type(quadrature_face_type), intent(in)    :: qr

      type(mesh_type),               pointer :: mesh              => null()
      class(reference_element_type), pointer :: reference_element => null()

      type(quadrature_face_proxy_type) :: qr_proxy
      integer :: cell, nlayers, nqp
      integer :: ndf_w2, ndf_w3, undf_w3, ndf_wtheta
      integer :: dim_w2, dim_w3, dim_wtheta

      integer,                   pointer :: cross_stencil_w3_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_w3 => null()
      integer                            :: cross_stencil_w3_size

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:), &
                                        basis_w3_face(:,:,:,:), &
                                        basis_wtheta_face(:,:,:,:)

      real(kind=r_def), pointer :: wqp(:,:) => null()

      type(operator_proxy_type) :: p2theta_proxy
      type(field_proxy_type)    :: exner_proxy

      integer(i_def) :: nfaces_h

      integer(i_def), pointer     :: adjacent_face(:,:) => null()
      real(r_def),    allocatable :: out_face_normal(:,:)

      !
      ! Initialise field proxies
      !
      p2theta_proxy = p2theta%get_proxy()
      exner_proxy = exner%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = p2theta_proxy%fs_from%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh              => p2theta%get_mesh()
      reference_element => mesh%get_reference_element()

      nfaces_h = reference_element%get_number_horizontal_faces()

      !
      ! Initialise qr values
      !
      qr_proxy = qr%get_quadrature_proxy()
      nqp =  qr_proxy%np_xyz
      wqp => qr_proxy%weights_xyz

      ndf_w2      = p2theta_proxy%fs_to%get_ndf( )
      dim_w2      = p2theta_proxy%fs_to%get_dim_space( )
      allocate(basis_w2_face(dim_w2,ndf_w2,nqp,nfaces_h))

      ndf_w3   = exner_proxy%vspace%get_ndf( )
      undf_w3  = exner_proxy%vspace%get_undf( )
      dim_w3   = exner_proxy%vspace%get_dim_space( )
      allocate(basis_w3_face(dim_w3,ndf_w3,nqp,nfaces_h))

      ndf_wtheta      = p2theta_proxy%fs_from%get_ndf( )
      dim_wtheta      = p2theta_proxy%fs_from%get_dim_space( )
      allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp,nfaces_h))

      call qr%compute_function(BASIS, p2theta_proxy%fs_to,   dim_w2,     ndf_w2,     basis_w2_face)
      call qr%compute_function(BASIS, exner_proxy%vspace,    dim_w3,     ndf_w3,     basis_w3_face)
      call qr%compute_function(BASIS, p2theta_proxy%fs_from, dim_wtheta, ndf_wtheta, basis_wtheta_face)

      ! Stencil maps
      cross_stencil_w3 => exner_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_w3_map => cross_stencil_w3%get_whole_dofmap()
      cross_stencil_w3_size = cross_stencil_w3%get_size()

      !
      ! Call kernels and communication routines
      !
      if (exner_proxy%is_dirty(depth=2)) call exner_proxy%halo_exchange(depth=2)
      !

      adjacent_face => mesh%get_adjacent_face()
      call reference_element%get_normals_to_out_faces( out_face_normal )

      do cell=1,mesh%get_last_halo_cell(1)

        call weighted_proj_2theta_bd_code( cell,                             &
                                           nlayers,                          &
                                           p2theta_proxy%ncell_3d,           &
                                           p2theta_proxy%local_stencil,      &
                                           exner_proxy%data,                 &
                                           s,                                &
                                           ndf_w2,                           &
                                           ndf_wtheta,                       &
                                           ndf_w3, undf_w3,                  &
                                           cross_stencil_w3_map(:,:,cell),   &
                                           cross_stencil_w3_size,            &
                                           nqp, wqp,                         &
                                           basis_w2_face,                    &
                                           basis_w3_face, basis_wtheta_face, &
                                           adjacent_face(:,cell),            &
                                           out_face_normal )
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate( basis_wtheta_face, basis_w3_face, basis_w2_face )
      !
    end subroutine invoke_weighted_proj_2theta_bd_kernel_type


  !-------------------------------------------------------------------------------
  !> Computation of 2d quadrature on faces not currently supported PSyClone,
  !> will be introduced in #793  by modifying quadrature tools developed in ticket #761.
  !> Kernel requires mesh information (adjacent_face) that will be implemented
  !> in lfric:#986 + psylcone:#18
  !> Invoke_weighted_proj_theta2_bd_kernel: Invoke the boundary part of the projection from w2 to wtheta operator of the lhs Helmholtz  
  subroutine invoke_weighted_proj_theta2_bd_kernel_type(ptheta2, theta, s, qr)

      use weighted_proj_theta2_bd_kernel_mod, only : weighted_proj_theta2_bd_code
      use mesh_mod,                           only : mesh_type
      use reference_element_mod,              only : reference_element_type
      use stencil_dofmap_mod,                 only : stencil_dofmap_type, STENCIL_CROSS

      implicit none

      type(field_type), intent(in)         :: theta
      type(operator_type), intent(inout)   :: ptheta2
      real(r_def), intent(in)              :: s
      type(quadrature_face_type), intent(in)    :: qr

      type(mesh_type),               pointer :: mesh              => null()
      class(reference_element_type), pointer :: reference_element => null()

      type(quadrature_face_proxy_type) :: qr_proxy
      integer :: cell, nlayers, nqp
      integer :: ndf_w2, ndf_wtheta, undf_wtheta
      integer :: dim_w2, dim_wtheta

      integer,                   pointer :: cross_stencil_wtheta_map(:,:,:) => null()
      type(stencil_dofmap_type), pointer :: cross_stencil_wtheta => null()

      integer                            :: cross_stencil_wtheta_size

      real(kind=r_def), allocatable  :: basis_w2_face(:,:,:,:), basis_wtheta_face(:,:,:,:)

      real(kind=r_def), pointer :: wqp(:,:) => null()

      type(operator_proxy_type) ::ptheta2_proxy
      type(field_proxy_type)    ::theta_proxy

      integer(i_def) :: nfaces_h

      integer(i_def), pointer     :: adjacent_face(:,:) => null()
      real(r_def),    allocatable :: normal_to_face(:,:)
      real(r_def),    allocatable :: out_face_normal(:,:)

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
      mesh              => ptheta2%get_mesh()
      reference_element => mesh%get_reference_element()

      nfaces_h = reference_element%get_number_horizontal_faces()

      !
      ! Initialise qr values
      !
      qr_proxy = qr%get_quadrature_proxy()
      nqp =  qr_proxy%np_xyz
      wqp => qr_proxy%weights_xyz

      ndf_w2      = ptheta2_proxy%fs_from%get_ndf( )
      dim_w2      = ptheta2_proxy%fs_from%get_dim_space( )
      allocate(basis_w2_face(dim_w2,ndf_w2,nqp,nfaces_h))

      ndf_wtheta      = ptheta2_proxy%fs_to%get_ndf( )
      dim_wtheta      = ptheta2_proxy%fs_to%get_dim_space( )
      undf_wtheta     = ptheta2_proxy%fs_to%get_undf()
      allocate(basis_wtheta_face(dim_wtheta,ndf_wtheta,nqp,nfaces_h))

      call qr%compute_function(BASIS, ptheta2_proxy%fs_from, dim_w2,     ndf_w2,     basis_w2_face)
      call qr%compute_function(BASIS, ptheta2_proxy%fs_to,   dim_wtheta, ndf_wtheta, basis_wtheta_face)

      ! Stencil maps
      cross_stencil_wtheta => ptheta2_proxy%fs_to%get_stencil_dofmap(STENCIL_CROSS, 1)
      cross_stencil_wtheta_map => cross_stencil_wtheta%get_whole_dofmap()
      cross_stencil_wtheta_size = cross_stencil_wtheta%get_size()

      !
      ! Call kernels and communication routines
      !
      if (theta_proxy%is_dirty(depth=2)) call theta_proxy%halo_exchange(depth=2)
      !
      adjacent_face => mesh%get_adjacent_face()
      call reference_element%get_normals_to_faces( normal_to_face )
      call reference_element%get_normals_to_out_faces( out_face_normal )

      do cell=1,mesh%get_last_halo_cell(1)

        call weighted_proj_theta2_bd_code( cell,                             &
                                         nlayers,                            &
                                         ptheta2_proxy%ncell_3d,             &
                                         ptheta2_proxy%local_stencil,        &
                                         theta_proxy%data,                   &
                                         s,                                  &
                                         ndf_w2,                             &
                                         ndf_wtheta, undf_wtheta,            &
                                         cross_stencil_wtheta_map(:,:,cell), &
                                         cross_stencil_wtheta_size,          &
                                         nqp, wqp,                           &
                                         basis_w2_face,                      &
                                         basis_wtheta_face,                  &
                                         adjacent_face(:,cell),              &
                                         normal_to_face, out_face_normal )
      end do
      !
      ! Deallocate basis arrays
      !
      deallocate( basis_wtheta_face, basis_w2_face )
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
  !> io_mod uses this routine. However, because io_mod is not a algorithm its currently
  !> not clear if it should call into psy. #1253 will address this point and remove
  !> once decided.
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

    type(stencil_dofmap_type), pointer  :: map_x_w3 => null()
    type(stencil_dofmap_type), pointer  :: map_y_w3 => null()
    integer, pointer                    :: map_w3(:) => null()
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
    map_x_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_approximation_stencil_extent)
    map_y_w3 => rho_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_approximation_stencil_extent)

    rho_stencil_size = map_x_w3%get_size()

    mesh => a0%get_mesh()

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_proxy%halo_exchange(depth=mesh%get_halo_depth())

    if (halo_depth_to_compute==0) then
      ncells_to_iterate = mesh%get_last_edge_cell()
    elseif (halo_depth_to_compute > 0) then
      ncells_to_iterate = mesh%get_last_halo_cell(halo_depth_to_compute)
    else
      call log_event( "Error: negative halo_depth_to_compute value in subgrid coeffs call", LOG_LEVEL_ERROR )
    endif

    !NOTE: The default looping limits for this type of field would be 
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and some halo cells.
    do cell = 1, ncells_to_iterate

      map_w3 => rho_proxy%vspace%get_cell_dofmap(cell)

      if (direction == x_direction) then
        if (nint(cell_orientation_proxy%data(map_w3(1))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1))) == 4) then
          stencil_map => map_y_w3%get_dofmap(cell)
        else
          stencil_map => map_x_w3%get_dofmap(cell)
        end if
      elseif (direction == y_direction) then
        if (nint(cell_orientation_proxy%data(map_w3(1))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1))) == 4) then
          stencil_map => map_x_w3%get_dofmap(cell)
        else
          stencil_map => map_y_w3%get_dofmap(cell)
        end if
      end if

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
!> invoke_subgrid_coeffs: Invoke the calculation of subgrid rho coefficients
!>                        The routine also includes a special type of halo
!>                        exchange where the values in the halos need to be
!>                        corrected. This is due to 1D direction updates of the
!>                        density field and the panels of the cubed-sphere having
!>                        different orientation.
  subroutine invoke_subgrid_coeffs_conservative( a0_x,                             &
                                                 a1_x,                             &
                                                 a2_x,                             &
                                                 a0_y,                             &
                                                 a1_y,                             &
                                                 a2_y,                             &
                                                 rho_x,                            &
                                                 rho_y,                            &
                                                 rho_x_halos_corrected,            &
                                                 rho_y_halos_corrected,            &
                                                 cell_orientation,                 &
                                                 rho_approximation_stencil_extent, &
                                                 dep_pt_stencil_extent  )

    use flux_direction_mod,               only: x_direction, y_direction
    use stencil_dofmap_mod,               only: stencil_dofmap_type, &
                                                STENCIL_1DX,         &
                                                STENCIL_1DY
    use subgrid_coeffs_kernel_mod,        only: subgrid_coeffs_code
    use subgrid_config_mod,               only: rho_approximation
    use mesh_mod,                         only: mesh_type
    use log_mod,                          only: log_event, LOG_LEVEL_ERROR
    use cosmic_halo_correct_x_kernel_mod, only: cosmic_halo_correct_x_code
    use cosmic_halo_correct_y_kernel_mod, only: cosmic_halo_correct_y_code

    implicit none

    type( field_type ), intent( inout ) :: a0_x
    type( field_type ), intent( inout ) :: a1_x
    type( field_type ), intent( inout ) :: a2_x
    type( field_type ), intent( inout ) :: a0_y
    type( field_type ), intent( inout ) :: a1_y
    type( field_type ), intent( inout ) :: a2_y
    type( field_type ), intent( in )    :: rho_x
    type( field_type ), intent( in )    :: rho_y
    type( field_type ), intent( inout ) :: rho_x_halos_corrected
    type( field_type ), intent( inout ) :: rho_y_halos_corrected
    type( field_type ), intent( in )    :: cell_orientation
    integer, intent(in)                 :: rho_approximation_stencil_extent
    integer, intent(in)                 :: dep_pt_stencil_extent

    type( field_proxy_type )            :: rho_x_proxy
    type( field_proxy_type )            :: rho_y_proxy
    type( field_proxy_type )            :: rho_x_halos_corrected_proxy
    type( field_proxy_type )            :: rho_y_halos_corrected_proxy
    type( field_proxy_type )            :: a0_x_proxy
    type( field_proxy_type )            :: a1_x_proxy
    type( field_proxy_type )            :: a2_x_proxy
    type( field_proxy_type )            :: a0_y_proxy
    type( field_proxy_type )            :: a1_y_proxy
    type( field_proxy_type )            :: a2_y_proxy
    type( field_proxy_type )            :: cell_orientation_proxy

    type(stencil_dofmap_type), pointer  :: map_x_w3 => null()
    type(stencil_dofmap_type), pointer  :: map_y_w3 => null()
    type(mesh_type), pointer            :: mesh => null()

    integer, pointer                    :: map_w3(:,:) => null()
    integer, pointer                    :: stencil_map(:,:) => null()
    integer                             :: rho_stencil_size
    integer                             :: cell
    integer                             :: nlayers
    integer                             :: ndf_w3
    integer                             :: undf_w3
    integer                             :: d
    logical                             :: swap
    integer                             :: ncells_to_iterate


    a0_x_proxy   = a0_x%get_proxy()
    a1_x_proxy   = a1_x%get_proxy()
    a2_x_proxy   = a2_x%get_proxy()
    a0_y_proxy   = a0_y%get_proxy()
    a1_y_proxy   = a1_y%get_proxy()
    a2_y_proxy   = a2_y%get_proxy()

    rho_x_proxy  = rho_x%get_proxy()
    rho_y_proxy  = rho_y%get_proxy()
    rho_x_halos_corrected_proxy = rho_x_halos_corrected%get_proxy()
    rho_y_halos_corrected_proxy = rho_y_halos_corrected%get_proxy()

    cell_orientation_proxy = cell_orientation%get_proxy()

    undf_w3 = rho_x_proxy%vspace%get_undf()
    ndf_w3  = rho_x_proxy%vspace%get_ndf()
    nlayers = rho_x_proxy%vspace%get_nlayers()

    mesh => rho_x%get_mesh()

    map_w3 => rho_x_proxy%vspace%get_whole_dofmap()

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_x_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_x_proxy%halo_exchange(depth=mesh%get_halo_depth())

    swap = .false.
    do d = 1,mesh%get_halo_depth()
      if (rho_y_proxy%is_dirty(depth=d)) swap = .true.
    end do
    if ( swap ) call rho_y_proxy%halo_exchange(depth=mesh%get_halo_depth())

    ! Loop over all core and halo cells.
    do cell=1,mesh%get_ncells_2d()

      call cosmic_halo_correct_x_code(  nlayers,                            &
                                        rho_x_halos_corrected_proxy%data,   &
                                        rho_x_proxy%data,                   &
                                        rho_y_proxy%data,                   &
                                        cell_orientation_proxy%data,        &
                                        ndf_w3,                             &
                                        undf_w3,                            &
                                        map_w3(:,cell))
    end do

    ! Loop over all core and halo cells.
    do cell=1,mesh%get_ncells_2d()

      call cosmic_halo_correct_y_code(  nlayers,                            &
                                        rho_y_halos_corrected_proxy%data,   &
                                        rho_x_proxy%data,                   &
                                        rho_y_proxy%data,                   &
                                        cell_orientation_proxy%data,        &
                                        ndf_w3,                             &
                                        undf_w3,                            &
                                        map_w3(:,cell))
    end do

    map_x_w3 => rho_x_proxy%vspace%get_stencil_dofmap(STENCIL_1DX,rho_approximation_stencil_extent)
    map_y_w3 => rho_y_proxy%vspace%get_stencil_dofmap(STENCIL_1DY,rho_approximation_stencil_extent)

    rho_stencil_size = map_x_w3%get_size()

    if (dep_pt_stencil_extent==0) then
      ncells_to_iterate = mesh%get_last_edge_cell()
    elseif (dep_pt_stencil_extent > 0) then
      ncells_to_iterate = mesh%get_last_halo_cell(dep_pt_stencil_extent)
    else
      call log_event( "Error: negative dep_pt_stencil_extent value in subgrid coeffs call", LOG_LEVEL_ERROR )
    endif

    ! Note stencil grid types are of the form:
    !                                   |5|
    !                                   |3|
    ! 1DX --> |4|2|1|3|5|  OR  1DY -->  |1|
    !                                   |2|
    !                                   |4|

    ! NOTE: The default looping limits for this type of field would be
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and some halo cells.

    ! Calculate a0_y, a1_y, a2_y associated with rho_y but in the x_direction
    do cell = 1, ncells_to_iterate
      if (nint(cell_orientation_proxy%data(map_w3(1,cell))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1,cell))) == 4) then
        stencil_map => map_y_w3%get_dofmap(cell)
      else
        stencil_map => map_x_w3%get_dofmap(cell)
      end if
      call subgrid_coeffs_code( nlayers,                                    &
                                rho_approximation,                          &
                                undf_w3,                                    &
                                rho_y_halos_corrected_proxy%data,           &
                                cell_orientation_proxy%data,                &
                                ndf_w3,                                     &
                                rho_stencil_size,                           &
                                stencil_map,                                &
                                x_direction,                                &
                                a0_y_proxy%data,                            &
                                a1_y_proxy%data,                            &
                                a2_y_proxy%data                             &
                                )
    end do


    ! Calculate a0_x, a1_x, a2_x associated with rho_x but in the y_direction
    do cell = 1, ncells_to_iterate
      if (nint(cell_orientation_proxy%data(map_w3(1,cell))) == 2 .or. nint(cell_orientation_proxy%data(map_w3(1,cell))) == 4) then
        stencil_map => map_x_w3%get_dofmap(cell)
      else
        stencil_map => map_y_w3%get_dofmap(cell)
      end if
      call subgrid_coeffs_code( nlayers,                                  &
                                rho_approximation,                        &
                                undf_w3,                                  &
                                rho_x_halos_corrected_proxy%data,         &
                                cell_orientation_proxy%data,              &
                                ndf_w3,                                   &
                                rho_stencil_size,                         &
                                stencil_map,                              &
                                y_direction,                              &
                                a0_x_proxy%data,                          &
                                a1_x_proxy%data,                          &
                                a2_x_proxy%data                           &
                                )
    end do
    call a0_x_proxy%set_dirty()
    call a1_x_proxy%set_dirty()
    call a2_x_proxy%set_dirty()
    call a0_y_proxy%set_dirty()
    call a1_y_proxy%set_dirty()
    call a2_y_proxy%set_dirty()

  end subroutine invoke_subgrid_coeffs_conservative


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
  ! in order to function correctly. See ticket #1058.
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

  mesh => u_n%get_mesh()
  !NOTE: The default looping limits for this type of field would be 
  ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
  ! in order to function correctly. See ticket #1058.
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
!> This call lies within psy-lite as the kernel makes use of stencils
!> and requires psyclone to generate the infrastructure correct calls for the 
!> horizontal looping when coloring is used
!> These will be implemented in psyclone issue #127
subroutine invoke_sample_poly_flux( flux, wind, density, stencil_extent )

  use mesh_mod,                    only: mesh_type
  use reference_element_mod,       only: reference_element_type
  use sample_poly_flux_kernel_mod, only: sample_poly_flux_code
  use stencil_dofmap_mod,          only: stencil_dofmap_type, STENCIL_CROSS

  implicit none

  type(field_type), intent(inout)      :: flux
  type(field_type), intent(in)         :: wind
  type(field_type), intent(in)         :: density
  integer, intent(in)                  :: stencil_extent

  type( field_proxy_type )  :: flux_proxy, wind_proxy, density_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()

  integer, pointer :: map_w2(:,:) => null(), map_w3(:,:) => null()
  integer, pointer :: stencil_map(:,:,:) => null()

  integer :: undf_w3, ndf_w3
  integer :: undf_w2, ndf_w2
  integer :: df_w2, df_wind
  integer :: cell
  integer :: nlayers
  integer :: stencil_size 
  integer :: d
  logical :: swap

  type(mesh_type),               pointer :: mesh              => null()
  class(reference_element_type), pointer :: reference_element => null()
  real(r_def),                   pointer :: nodes_w2(:,:)     => null()

  real(r_def), allocatable :: basis_w2(:,:,:)
  real(r_def), allocatable :: out_face_normal(:,:)

  integer :: dim_w2

  integer :: colour, ncolour
  integer, pointer :: cmap(:,:) =>null(), ncp_colour(:) => null()

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
  map_w3 => density_proxy%vspace%get_whole_dofmap()

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
  reference_element => mesh%get_reference_element()
  call reference_element%get_normals_to_out_faces( out_face_normal )

  ! Look-up colour map
  call flux_proxy%vspace%get_colours(ncolour, ncp_colour, cmap)
  do colour=1,ncolour
    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    do cell=1,mesh%get_last_halo_cell_per_colour_any(colour,1)


      call sample_poly_flux_code( nlayers,                      &
                                  flux_proxy%data,              &
                                  wind_proxy%data,              &
                                  density_proxy%data,           &
                                  stencil_size,                 &
                                  stencil_map(:,:,cmap(colour, cell)),       &
                                  ndf_w2,                       &
                                  undf_w2,                      &
                                  map_w2(:,cmap(colour, cell)), &
                                  basis_w2,                     &
                                  ndf_w3,                       &
                                  undf_w3,                      &
                                  map_w2(:,cmap(colour, cell)), &
                                  out_face_normal )
    end do 
    !$omp end do
    !$omp end parallel
  end do

  call flux_proxy%set_dirty()

end subroutine invoke_sample_poly_flux

!-------------------------------------------------------------------------------
!> This call lies within psy-lite as the kernel makes use of stencils
!> and requires psyclone to generate the infrastructure correct calls for the 
!> horizontal looping when coloring is used.
!> These will be implemented in psyclone issue #127
subroutine invoke_sample_poly_adv( adv, tracer, wind, mt_lumped_inv, &
                                   chi, stencil_extent)

  use sample_poly_adv_kernel_mod, only: sample_poly_adv_code
  use stencil_dofmap_mod,         only: stencil_dofmap_type, STENCIL_CROSS
  use mesh_mod,                   only: mesh_type

  implicit none

  type(field_type), intent(inout)      :: adv
  type(field_type), intent(in)         :: tracer
  type(field_type), intent(in)         :: wind
  type(field_type), intent(in)         :: mt_lumped_inv
  type(field_type), intent(in)         :: chi(3)
  integer, intent(in)                  :: stencil_extent

  type( field_proxy_type )  :: adv_proxy, wind_proxy, tracer_proxy, chi_proxy(3)
  type( field_proxy_type )  :: mt_lumped_inv_proxy

  type(stencil_dofmap_type), pointer :: stencil => null()
  type(stencil_dofmap_type), pointer :: stencil_wx => null()

  integer, pointer :: map_w2(:,:) => null(), map_wt(:,:) => null(), map_wx(:,:) => null()
  integer, pointer :: stencil_map(:,:,:) => null()
  integer, pointer :: stencil_map_wx(:,:,:) => null()

  integer :: undf_wt, ndf_wt
  integer :: undf_w2, ndf_w2
  integer :: undf_wx, ndf_wx
  integer :: ndf_adv
  integer :: df_adv, df_w2, df_wx
  integer :: cell
  integer :: nlayers
  integer :: stencil_size, stencil_size_wx
  integer :: d
  logical :: swap
  type(mesh_type), pointer :: mesh => null()
  real(kind=r_def), pointer :: nodes_adv(:,:) => null()

  real(kind=r_def), allocatable :: basis_w2(:,:,:)
  real(kind=r_def), allocatable :: basis_wx(:,:,:)
  real(kind=r_def), allocatable :: diff_basis_wx(:,:,:)
  integer :: dim_w2
  integer :: dim_wx, diff_dim_wx

  integer :: colour, ncolour
  integer, pointer :: cmap(:,:) =>null(), ncp_colour(:) => null()

  adv_proxy    = adv%get_proxy()
  wind_proxy   = wind%get_proxy()
  tracer_proxy = tracer%get_proxy()
  mt_lumped_inv_proxy = mt_lumped_inv%get_proxy()
  chi_proxy(1) = chi(1)%get_proxy()
  chi_proxy(2) = chi(2)%get_proxy()
  chi_proxy(3) = chi(3)%get_proxy()

  ndf_adv = adv_proxy%vspace%get_ndf()
  nodes_adv => adv_proxy%vspace%get_nodes()

  ndf_wt  = tracer_proxy%vspace%get_ndf()
  undf_wt = tracer_proxy%vspace%get_undf()

  ndf_w2  = wind_proxy%vspace%get_ndf()
  undf_w2 = wind_proxy%vspace%get_undf()
  dim_w2  = wind_proxy%vspace%get_dim_space()

  ndf_wx  = chi_proxy(1)%vspace%get_ndf()
  undf_wx = chi_proxy(1)%vspace%get_undf()
  dim_wx  = chi_proxy(1)%vspace%get_dim_space()
  diff_dim_wx = chi_proxy(1)%vspace%get_dim_space_diff()

  ! Evaluate the basis function
  allocate (basis_w2(dim_w2, ndf_w2, ndf_adv))
  allocate (basis_wx(dim_wx, ndf_wx, ndf_adv))
  allocate (diff_basis_wx(diff_dim_wx, ndf_wx, ndf_adv))
  do df_adv = 1, ndf_adv
    do df_w2 = 1, ndf_w2
      basis_w2(:,df_w2,df_adv) = wind_proxy%vspace%call_function(BASIS,df_w2,nodes_adv(:,df_adv))
    end do
    do df_wx = 1, ndf_wx
      basis_wx(:,df_wx,df_adv) = chi_proxy(1)%vspace%call_function(BASIS,df_wx,nodes_adv(:,df_adv))
      diff_basis_wx(:,df_wx,df_adv) = chi_proxy(1)%vspace%call_function(DIFF_BASIS,df_wx,nodes_adv(:,df_adv))
    end do
  end do

  nlayers = adv_proxy%vspace%get_nlayers()

  stencil => tracer_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                    stencil_extent)
  stencil_size = stencil%get_size()
  stencil_map => stencil%get_whole_dofmap()

  stencil_wx => chi_proxy(1)%vspace%get_stencil_dofmap(STENCIL_CROSS, &
                                                       stencil_extent)
  stencil_size_wx = stencil_wx%get_size()
  stencil_map_wx => stencil_wx%get_whole_dofmap()

  if(wind_proxy%is_dirty(depth=1) ) then
    call wind_proxy%halo_exchange(depth=1)
  end if
  swap = .false.
  do d = 1,stencil_extent+1
    if(tracer_proxy%is_dirty(depth=d) ) swap = .true.
  end do
  if ( swap ) call tracer_proxy%halo_exchange(depth=stencil_extent+1)
  swap = .false.
  do d = 1,stencil_extent+1
    if(chi_proxy(1)%is_dirty(depth=d) ) swap = .true.
  end do
  if ( swap ) then
    call chi_proxy(1)%halo_exchange(depth=stencil_extent+1)
    call chi_proxy(2)%halo_exchange(depth=stencil_extent+1)
    call chi_proxy(3)%halo_exchange(depth=stencil_extent+1)
  end if

  map_w2 => wind_proxy%vspace%get_whole_dofmap()
  map_wt => tracer_proxy%vspace%get_whole_dofmap()
  map_wx => chi_proxy(1)%vspace%get_whole_dofmap()

  mesh => adv%get_mesh()

  ! Look-up colour map
  call adv_proxy%vspace%get_colours(ncolour, ncp_colour, cmap)
  do colour=1,ncolour
    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    do cell=1,mesh%get_last_edge_cell_per_colour(colour)

      call sample_poly_adv_code( nlayers,                     &
                                 adv_proxy%data,              &
                                 tracer_proxy%data,           &
                                 stencil_size,                &
                                 stencil_map(:,:,cmap(colour, cell)),       &
                                 wind_proxy%data,             &
                                 mt_lumped_inv_proxy%data,    &
                                 chi_proxy(1)%data,           &
                                 stencil_size_wx,             &
                                 stencil_map_wx(:,:,cmap(colour, cell)),    &
                                 chi_proxy(2)%data,           &
                                 chi_proxy(3)%data,           &
                                 ndf_wt,                      &
                                 undf_wt,                     &
                                 map_wt(:,cmap(colour, cell)),&
                                 ndf_w2,                      &
                                 undf_w2,                     &
                                 map_w2(:,cmap(colour, cell)),&
                                 basis_w2,                    &
                                 ndf_wx,                      &
                                 undf_wx,                     &
                                 map_wx(:,cmap(colour, cell)),&
                                 basis_wx,                    &
                                 diff_basis_wx                &
                                 )
    end do
    !$omp end do
    !$omp end parallel
  end do

  call adv_proxy%set_dirty()

  deallocate( diff_basis_wx, basis_wx, basis_w2 )

end subroutine invoke_sample_poly_adv

  !------------------------------------------------------------------------------
  ! Needs correct loop limits in the presence of colouring for psyclone
  ! implementation
  subroutine invoke_tracer_viscosity( theta_inc, theta_n ,t_stencil_size, chi )

    use tracer_viscosity_kernel_mod, only : tracer_viscosity_code
    use mesh_mod,                    only : mesh_type 
    use stencil_dofmap_mod,          only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( inout ) :: theta_inc, theta_n, chi(3)
    integer,            intent( in )    :: t_stencil_size
    type( mesh_type ), pointer       :: mesh => null()
    integer                 :: cell, nlayers
    integer                 :: ndf_wt, ndf_chi
    integer                 :: undf_wt, undf_chi

    type( field_proxy_type )  :: theta_inc_proxy, theta_n_proxy, &
                                 chi_proxy(3)

    integer, pointer        :: map_chi(:,:) => null(),map_t(:,:) => null()

    type(stencil_dofmap_type), pointer :: cross_stencil_wt => null()


    integer, pointer        :: cross_stencil_wt_map(:,:,:) => null()
    integer                 :: cross_stencil_wt_size

    theta_inc_proxy = theta_inc%get_proxy()
    theta_n_proxy   = theta_n%get_proxy()

    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    nlayers = theta_inc_proxy%vspace%get_nlayers()
    ndf_wt  = theta_inc_proxy%vspace%get_ndf( )
    undf_wt = theta_inc_proxy%vspace%get_undf()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf()

    mesh => theta_inc%get_mesh()

    if (theta_inc_proxy%is_dirty(depth=1)) CALL theta_inc_proxy%halo_exchange(depth=1)
    if (   chi_proxy(1)%is_dirty(depth=1)) CALL    chi_proxy(1)%halo_exchange(depth=1)
    if (   chi_proxy(2)%is_dirty(depth=1)) CALL    chi_proxy(2)%halo_exchange(depth=1)
    if (   chi_proxy(3)%is_dirty(depth=1)) CALL    chi_proxy(3)%halo_exchange(depth=1)

    if (theta_n_proxy%is_dirty(depth=1+1)) CALL theta_n_proxy%halo_exchange(depth=1+1)

    cross_stencil_wt => theta_n_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, t_stencil_size)
    cross_stencil_wt_map => cross_stencil_wt%get_whole_dofmap()
    cross_stencil_wt_size = cross_stencil_wt%get_size()
    
    map_chi => chi_proxy(1)%vspace%get_whole_dofmap( )
    map_t => theta_inc_proxy%vspace%get_whole_dofmap( )

    do cell = 1, mesh%get_last_edge_cell()

       call tracer_viscosity_code(nlayers,                   &
                                  theta_inc_proxy%data,      &
                                  theta_n_proxy%data,        &
                                  cross_stencil_wt_size,          &
                                  cross_stencil_wt_map(:,:,cell), &
                                  chi_proxy(1)%data,         &
                                  chi_proxy(2)%data,         &
                                  chi_proxy(3)%data,         &
                                  ndf_wt, undf_wt, map_t(:,cell),    &
                                  ndf_chi, undf_chi, map_chi(:,cell) &
                                 )
    end do
     
   call theta_inc_proxy%set_dirty()

  end subroutine invoke_tracer_viscosity

  ! Needs correct loop limits in the presence of colouring for psyclone
  ! implementation
  subroutine invoke_momentum_viscosity( u_inc, u_n, u_stencil_size, chi )
   
    use momentum_viscosity_kernel_mod, only : momentum_viscosity_code
    use mesh_mod,                      only : mesh_type 
    use stencil_dofmap_mod,            only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( inout ) :: u_inc, u_n, chi(3)
    integer,            intent( in )    :: u_stencil_size
    type( mesh_type ), pointer          :: mesh => null()
    integer                             :: cell, nlayers
    integer                             :: ndf_w2, ndf_chi
    integer                             :: undf_w2, undf_chi

    type( field_proxy_type )            :: u_inc_proxy, u_n_proxy, &
                                           chi_proxy(3)

    integer, pointer        :: map_chi(:,:) => null(), map_u(:,:) => null()

    type(stencil_dofmap_type), pointer :: cross_stencil_w2 => null()

    integer, pointer        :: cross_stencil_w2_map(:,:,:) => null()
    integer                 :: cross_stencil_w2_size


    u_inc_proxy = u_inc%get_proxy()
    u_n_proxy   = u_n%get_proxy()

    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    nlayers = u_inc_proxy%vspace%get_nlayers()
    ndf_w2  = u_inc_proxy%vspace%get_ndf( )
    undf_w2 = u_inc_proxy%vspace%get_undf()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf()

    mesh => u_inc%get_mesh()

    if (    u_inc_proxy%is_dirty(depth=1)) CALL     u_inc_proxy%halo_exchange(depth=1)
    if (   chi_proxy(1)%is_dirty(depth=1)) CALL    chi_proxy(1)%halo_exchange(depth=1)
    if (   chi_proxy(2)%is_dirty(depth=1)) CALL    chi_proxy(2)%halo_exchange(depth=1)
    if (   chi_proxy(3)%is_dirty(depth=1)) CALL    chi_proxy(3)%halo_exchange(depth=1)

    if (    u_n_proxy%is_dirty(depth=1+1)) CALL     u_n_proxy%halo_exchange(depth=1+1)

    cross_stencil_w2 => u_n_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, u_stencil_size)
    cross_stencil_w2_map => cross_stencil_w2%get_whole_dofmap()
    cross_stencil_w2_size = cross_stencil_w2%get_size()

    map_chi => chi_proxy(1)%vspace%get_whole_dofmap( )
    map_u => u_inc_proxy%vspace%get_whole_dofmap( )

    do cell = 1, mesh%get_last_halo_cell(1)

       call momentum_viscosity_code(nlayers,                   &
                                    u_inc_proxy%data,          &
                                    u_n_proxy%data,            &
                                    cross_stencil_w2_size,          &
                                    cross_stencil_w2_map(:,:,cell), &
                                    chi_proxy(1)%data,         &
                                    chi_proxy(2)%data,         &
                                    chi_proxy(3)%data,         &
                                    ndf_w2, undf_w2, map_u(:,cell),    &
                                    ndf_chi, undf_chi, map_chi(:,cell) &
                                   )
    end do
     
   call u_inc_proxy%set_dirty()

  end subroutine invoke_momentum_viscosity

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
  !> PSyclone support is required - This routine will be removed by #1261 
  subroutine invoke_mpi_detj_at_w2(detj_at_w2, chi)
    use calc_detj_at_w2_kernel_mod, only : calc_detj_at_w2_code
    use mesh_mod, only                   : mesh_type

    implicit none

    type(field_type), intent(inout) :: detj_at_w2
    type(field_type), intent(in)    :: chi(3)

    type(field_proxy_type)          :: detj_at_w2_p, chi_p(3)
    type(mesh_type), pointer        :: mesh => null()

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

    mesh => detj_at_w2%get_mesh()

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
 
  !-------------------------------------------------------------------------------
  !> This routine is called from psykal_lite due to the variable cell_orientation
  !> being passed into the kernel.
  !> The cell_orientation field should not be halo exchanged across panels as the
  !> orientation of cells is local to its own panel on the cubed-sphere.
  subroutine invoke_fv_divergence( mass_divergence,          &
                                   mass_flux_x,              &
                                   mass_flux_y,              &
                                   cell_orientation,         &
                                   direction )

    use mesh_mod,                         only : mesh_type
    use fv_divergence_kernel_mod,         only : fv_divergence_code
    use flux_direction_mod,               only : x_direction, y_direction

    implicit none

    type(field_type), intent(inout) :: mass_divergence
    type(field_type), intent(in)    :: mass_flux_x
    type(field_type), intent(in)    :: mass_flux_y
    type(field_type), intent(in)    :: cell_orientation
    integer, intent(in)             :: direction

    type(mesh_type), pointer        :: mesh => null()

    integer, pointer                :: map_w3(:,:) => null()
    integer, pointer                :: map_w2(:,:) => null()

    type(field_proxy_type)          :: cell_orientation_proxy
    type(field_proxy_type)          :: mass_flux_x_proxy, mass_flux_y_proxy
    type(field_proxy_type)          :: mass_divergence_proxy

    integer                         :: cell, nlayers
    integer                         :: ndf_w3, undf_w3
    integer                         :: ndf_w2, undf_w2

    cell_orientation_proxy           = cell_orientation%get_proxy()
    mass_divergence_proxy            = mass_divergence%get_proxy()
    mass_flux_x_proxy                = mass_flux_x%get_proxy()
    mass_flux_y_proxy                = mass_flux_y%get_proxy()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = mass_flux_x_proxy%vspace%get_ndf( )
    undf_w2 = mass_flux_x_proxy%vspace%get_undf()

    mesh   => mass_flux_x%get_mesh()
    map_w3 => cell_orientation_proxy%vspace%get_whole_dofmap()
    map_w2 => mass_flux_x_proxy%vspace%get_whole_dofmap()

    ! There is no automatic halo exchange on purpose at the moment. Since if there
    ! was then the x and y directional components would not be respected due to
    ! different panel orientations on the cubed-sphere.
    ! A ticket, #1147, has been created which addresses this issue of dealing
    ! with panel orientation when halo exchanging W2 fields.
    ! A similar implementation was made for invoke_subgrid_coeffs_conservative
    ! which dealt with W3 fields and ticket #1087 implemented this.
    if (direction == x_direction) then

      do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only

        call  fv_divergence_code(  nlayers,                             &
                                   mass_divergence_proxy%data,          &
                                   cell_orientation_proxy%data,         &
                                   mass_flux_x_proxy%data,              &
                                   direction,                           &
                                   ndf_w3,                              &
                                   undf_w3,                             &
                                   map_w3(:,cell),                      &
                                   ndf_w2,                              &
                                   undf_w2,                             &
                                   map_w2(:,cell) )

      end do

    elseif (direction == y_direction) then

      do cell=1,mesh%get_last_edge_cell() ! Loop over core cells only

        call  fv_divergence_code(  nlayers,                             &
                                   mass_divergence_proxy%data,          &
                                   cell_orientation_proxy%data,         &
                                   mass_flux_y_proxy%data,              &
                                   direction,                           &
                                   ndf_w3,                              &
                                   undf_w3,                             &
                                   map_w3(:,cell),                      &
                                   ndf_w2,                              &
                                   undf_w2,                             &
                                   map_w2(:,cell) )

      end do

    end if

  end subroutine invoke_fv_divergence

  !-------------------------------------------------------------------------------
  !> This kernel routine is in psykal_lite due to the use of cell_orientation
  !> variable which cannot be halo exchanged and also the depth of the halos to
  !> loop over includes all halo depths, i.e. all cells in the core and halo.
  subroutine invoke_extract_xy(x_field_out,y_field_out,w2_field_in,cell_orientation)

    use extract_x_kernel_mod,        only : extract_x_code
    use extract_y_kernel_mod,        only : extract_y_code
    use log_mod,                     only : log_event, log_scratch_space,     &
                                            LOG_LEVEL_INFO
    use mesh_mod,                    only : mesh_type

    implicit none

    type(field_type), intent(inout)    :: x_field_out
    type(field_type), intent(inout)    :: y_field_out
    type(field_type), intent(in)       :: w2_field_in
    type(field_type), intent(in)       :: cell_orientation

    type(field_proxy_type)             :: cell_orientation_proxy
    type(field_proxy_type)             :: x_field_out_proxy
    type(field_proxy_type)             :: y_field_out_proxy
    type(field_proxy_type)             :: w2_field_in_proxy

    integer                            :: cell, nlayers
    integer                            :: ndf_w3, undf_w3
    integer                            :: ndf_w2, undf_w2
    integer, pointer                   :: map_w3(:,:) => null()
    integer, pointer                   :: map_w2(:,:) => null()
    integer                            :: halo_depth

    type(mesh_type), pointer           :: mesh => null()

    cell_orientation_proxy = cell_orientation%get_proxy()
    x_field_out_proxy      = x_field_out%get_proxy()
    y_field_out_proxy      = y_field_out%get_proxy()
    w2_field_in_proxy      = w2_field_in%get_proxy()

    nlayers = cell_orientation_proxy%vspace%get_nlayers()
    ndf_w3  = cell_orientation_proxy%vspace%get_ndf( )
    undf_w3 = cell_orientation_proxy%vspace%get_undf()
    ndf_w2  = x_field_out_proxy%vspace%get_ndf( )
    undf_w2 = x_field_out_proxy%vspace%get_undf()

    mesh   => x_field_out%get_mesh()
    map_w2 => x_field_out_proxy%vspace%get_whole_dofmap()
    map_w3 => cell_orientation_proxy%vspace%get_whole_dofmap()

    halo_depth = mesh%get_halo_depth()
    call w2_field_in_proxy%halo_exchange(depth=halo_depth)

    ! Extract the x-component of the W2 field
    do cell = 1, mesh%get_ncells_2d() ! Loop over core and halo cells

      call extract_x_code(  nlayers,                             &
                            cell_orientation_proxy%data,         &
                            w2_field_in_proxy%data,              &
                            x_field_out_proxy%data,              &
                            undf_w3,                             &
                            ndf_w3,                              &
                            map_w3(:,cell),                      &
                            undf_w2,                             &
                            ndf_w2,                              &
                            map_w2(:,cell) )
    end do

    ! Extract the y-component of the W2 field
    do cell = 1, mesh%get_ncells_2d() ! Loop over core and halo cells

      call extract_y_code(  nlayers,                             &
                            cell_orientation_proxy%data,         &
                            w2_field_in_proxy%data,              &
                            y_field_out_proxy%data,              &
                            undf_w3,                             &
                            ndf_w3,                              &
                            map_w3(:,cell),                      &
                            undf_w2,                             &
                            ndf_w2,                              &
                            map_w2(:,cell) )
    end do

  end subroutine invoke_extract_xy

  !-------------------------------------------------------------------------------
  ! Ticket #1156. Stephen Pring
  ! This code is implemented in psykal-lite because the cells to
  ! iterate over include all core cells and all halo cells. At present, the default
  ! iteration is over core cells and a halo depth of 1. The cosmic transport scheme
  ! uses a larger halo depth and this routine requires iteration over all values
  ! in the halo as well.
  subroutine invoke_cosmic_departure_wind(dep_wind_x,dep_wind_y,u_piola_x,u_piola_y,detj_at_w2,direction)
    use cosmic_departure_wind_kernel_mod, only: cosmic_departure_wind_code
    use mesh_mod,                         only: mesh_type
    use flux_direction_mod,               only: x_direction, y_direction
    use log_mod,                          only: log_event, LOG_LEVEL_ERROR

    implicit none

    type(field_type), intent(inout)      :: dep_wind_x, dep_wind_y
    type(field_type), intent(in)         :: u_piola_x, u_piola_y
    integer, intent(in)                  :: direction
    type(field_type), intent(in)         :: detj_at_w2

    type(field_proxy_type) :: dep_wind_x_p, dep_wind_y_p
    type(field_proxy_type) :: u_piola_x_p, u_piola_y_p
    type(field_proxy_type) :: detj_at_w2_p

    integer                 :: cell, nlayers
    integer                 :: ndf_w2
    integer                 :: undf_w2
    integer, pointer        :: map(:) => null()

    type(mesh_type), pointer :: mesh => null()
    integer :: halo_depth

    u_piola_x_p = u_piola_x%get_proxy()
    u_piola_y_p = u_piola_y%get_proxy()
    dep_wind_x_p = dep_wind_x%get_proxy()
    dep_wind_y_p = dep_wind_y%get_proxy()
    detj_at_w2_p = detj_at_w2%get_proxy()

    mesh => u_piola_x%get_mesh()
    halo_depth = mesh%get_halo_depth()
    call detj_at_w2_p%halo_exchange(depth=halo_depth)

    nlayers = u_piola_x_p%vspace%get_nlayers()

    ndf_w2  = u_piola_x_p%vspace%get_ndf()
    undf_w2 = u_piola_x_p%vspace%get_undf()


    ! NOTE: The default looping limits for this type of field would be 
    ! mesh%get_last_halo_cell(1) but this kernel requires a modified loop limit
    ! in order to function correctly. See ticket #1058.
    ! The kernel loops over all core and all halo cells.
    if (direction == x_direction) then
      do cell = 1,mesh%get_ncells_2d()
         map     => u_piola_x_p%vspace%get_cell_dofmap( cell )
         call cosmic_departure_wind_code( nlayers,                                  &
                                          dep_wind_x_p%data,                        &
                                          u_piola_x_p%data,                         &
                                          detj_at_w2_p%data,                        &
                                          ndf_w2, undf_w2, map,                     &
                                          direction                                 &
                                           )
      end do
    elseif (direction == y_direction) then
      do cell = 1,mesh%get_ncells_2d()
         map     => u_piola_y_p%vspace%get_cell_dofmap( cell )
         call cosmic_departure_wind_code( nlayers,                                  &
                                          dep_wind_y_p%data,                        &
                                          u_piola_y_p%data,                         &
                                          detj_at_w2_p%data,                        &
                                          ndf_w2, undf_w2, map,                     &
                                          direction                                 &
                                           )
      end do
    else
      call log_event("Direction incorrectly specified in invoke_cosmic_departure_wind",LOG_LEVEL_ERROR)
    end if

    call dep_wind_x_p%set_dirty()
    call dep_wind_y_p%set_dirty()

  end subroutine invoke_cosmic_departure_wind


  !-------------------------------------------------------------------------------
  ! Ticket #1156. Stephen Pring
  ! This code is implemented in psykal-lite because the cells to
  ! iterate over include all core cells and all halo cells. At present, the default
  ! iteration is over core cells and a halo depth of 1. The cosmic transport scheme
  ! uses a larger halo depth and this routine requires iteration over all values
  ! in the halo as well.
  subroutine invoke_correct_cosmic_wind(wind_x_out,                   &
                                        wind_y_out,                   &
                                        departure_wind_x_in,          &
                                        departure_wind_y_in,          &
                                        orientation_of_cells,         &
                                        direction)

    use correct_cosmic_wind_kernel_mod, only: correct_cosmic_wind_code
    use flux_direction_mod,             only: x_direction, y_direction
    use mesh_mod,                       only: mesh_type
    use log_mod,                        only: log_event, LOG_LEVEL_ERROR

    implicit none

    type(field_type), intent(inout) :: wind_x_out, wind_y_out
    type(field_type), intent(in)    :: departure_wind_x_in,departure_wind_y_in,orientation_of_cells
    integer,          intent(in)    :: direction

    type(field_proxy_type) :: wind_x_in_proxy, wind_y_in_proxy
    type(field_proxy_type) :: wind_x_out_proxy, wind_y_out_proxy
    type(field_proxy_type) :: chi_p(3), orientation_proxy

    integer                 :: cell, nlayers
    integer                 :: ndf_w2, ndf_w3
    integer                 :: undf_w2, undf_w3
    integer, pointer        :: map_w3(:) => null()
    integer, pointer        :: map_w2(:) => null()

    type(mesh_type), pointer :: mesh => null()

    wind_x_in_proxy = departure_wind_x_in%get_proxy()
    wind_y_in_proxy = departure_wind_y_in%get_proxy()
    wind_x_out_proxy = wind_x_out%get_proxy()
    wind_y_out_proxy = wind_y_out%get_proxy()
    orientation_proxy = orientation_of_cells%get_proxy()

    nlayers = orientation_proxy%vspace%get_nlayers()
    ndf_w3  = orientation_proxy%vspace%get_ndf()
    undf_w3 = orientation_proxy%vspace%get_undf()

    ndf_w2  = wind_x_in_proxy%vspace%get_ndf( )
    undf_w2 = wind_x_in_proxy%vspace%get_undf()

    mesh => orientation_of_cells%get_mesh()


    if (direction == x_direction) then
      do cell = 1, mesh%get_ncells_2d()
        map_w3 => orientation_proxy%vspace%get_cell_dofmap(cell)
        map_w2 => wind_x_in_proxy%vspace%get_cell_dofmap(cell)

        call correct_cosmic_wind_code(  nlayers,                        &
                                        wind_x_out_proxy%data,          &
                                        wind_x_in_proxy%data,           &
                                        orientation_proxy%data,         &
                                        undf_w2,                        &
                                        ndf_w2,                         &
                                        map_w2,                         &
                                        undf_w3,                        &
                                        ndf_w3,                         &
                                        map_w3,                         &
                                        direction )

      end do
    elseif (direction == y_direction) then
      do cell = 1, mesh%get_ncells_2d()
        map_w3 => orientation_proxy%vspace%get_cell_dofmap(cell)
        map_w2 => wind_x_in_proxy%vspace%get_cell_dofmap(cell)

        call correct_cosmic_wind_code(  nlayers,                        &
                                        wind_y_out_proxy%data,          &
                                        wind_y_in_proxy%data,           &
                                        orientation_proxy%data,         &
                                        undf_w2,                        &
                                        ndf_w2,                         &
                                        map_w2,                         &
                                        undf_w3,                        &
                                        ndf_w3,                         &
                                        map_w3,                         &
                                        direction )

      end do
    else
      call log_event("Direction incorrectly specified in invoke_correct_cosmic_wind",LOG_LEVEL_ERROR)
    end if

  end subroutine invoke_correct_cosmic_wind

  !------------------------------------------------------------------------------
  !> invoke_initial_exner_sample_kernel: invoke the density initialization for a generic space
  !> Computation of nodal basis function for coordinates chi not currently supported PSyClone,
  !> will be introduced in modification of quadrature strategy, see ticket #723.

  subroutine invoke_initial_exner_sample_kernel( exner, chi, time )

    use initial_exner_sample_kernel_mod, only : initial_exner_sample_code
    use mesh_mod,                      only : mesh_type
    implicit none

    type( field_type ), intent( inout )  :: exner
    type( field_type ), intent( in )     :: chi(3)
    real(kind=r_def),   intent( in )     :: time

    integer          :: cell
    integer          :: ndf_w3, undf_w3, &
                        ndf_chi, undf_chi, dim_chi
    integer, pointer :: map_w3(:)  => null()
    integer, pointer :: map_chi(:) => null()

    type( field_proxy_type ) :: exner_proxy
    type( field_proxy_type ) :: chi_proxy(3)

    real(kind=r_def), allocatable :: basis_chi(:,:,:)

    type(mesh_type), pointer :: mesh => null()

    integer          :: df_w3, df_chi
    real(kind=r_def), pointer :: nodes_w3(:,:) => null()

    exner_proxy    = exner%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ndf_w3       = exner_proxy%vspace%get_ndf( )
    undf_w3      = exner_proxy%vspace%get_undf( )

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi  = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( )

    allocate( basis_chi(dim_chi, ndf_chi, ndf_w3) )

    nodes_w3 => exner_proxy%vspace%get_nodes()
    do df_w3 = 1, ndf_w3
      do df_chi = 1, ndf_chi
        basis_chi(:,df_chi,df_w3) = &
                  chi_proxy(1)%vspace%call_function(BASIS,df_chi,nodes_w3(:,df_w3))
      end do
    end do

    if (chi_proxy(1)%is_dirty(depth=1)) call chi_proxy(1)%halo_exchange(depth=1)
    if (chi_proxy(2)%is_dirty(depth=1)) call chi_proxy(2)%halo_exchange(depth=1)
    if (chi_proxy(3)%is_dirty(depth=1)) call chi_proxy(3)%halo_exchange(depth=1)

    mesh => exner%get_mesh()
    do cell = 1, mesh%get_last_halo_cell(1)

      map_w3 => exner_proxy%vspace%get_cell_dofmap( cell )
      map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      call initial_exner_sample_code(                                  &
                                     exner_proxy%vspace%get_nlayers(), &
                                     ndf_w3,                           &
                                     undf_w3,                          &
                                     map_w3,                           &
                                     exner_proxy%data,                 &
                                     ndf_chi,                          &
                                     undf_chi,                         &
                                     map_chi,                          &
                                     basis_chi,                        &
                                     chi_proxy(1)%data,                &
                                     chi_proxy(2)%data,                &
                                     chi_proxy(3)%data,                &
                                     time                              &
                                     )
    end do

    call exner_proxy%set_dirty()

    deallocate( basis_chi )
  end subroutine invoke_initial_exner_sample_kernel

  !------------------------------------------------------------------------------
  ! Ticket #1470. Kirsty Hanley
  ! Needs correct loop limits in the presence of colouring for psyclone
  ! implementation
  subroutine invoke_tracer_smagorinsky_diff( theta_inc, theta_n, t_stencil_size, delta )

    use tracer_smagorinsky_diff_kernel_mod, only : tracer_smagorinsky_diff_code
    use mesh_mod,                    only : mesh_type 
    use stencil_dofmap_mod,          only : stencil_dofmap_type, STENCIL_CROSS

    implicit none

    type( field_type ), intent( inout ) :: theta_inc, theta_n, delta
    integer(kind=i_def),intent( in )    :: t_stencil_size
    type( mesh_type ), pointer          :: mesh => null()
    integer(kind=i_def)                 :: cell, nlayers
    integer(kind=i_def)                 :: ndf_wt
    integer(kind=i_def)                 :: undf_wt

    type( field_proxy_type )  :: theta_inc_proxy, theta_n_proxy, &
                                 delta_proxy

    integer(kind=i_def), pointer        :: map_t(:,:) => null()

    type(stencil_dofmap_type), pointer  :: cross_stencil_wt => null()


    integer(kind=i_def), pointer        :: cross_stencil_wt_map(:,:,:) => null()
    integer(kind=i_def)                 :: cross_stencil_wt_size

    theta_inc_proxy = theta_inc%get_proxy()
    theta_n_proxy   = theta_n%get_proxy()

    delta_proxy = delta%get_proxy()

    nlayers = theta_inc_proxy%vspace%get_nlayers()
    ndf_wt  = theta_inc_proxy%vspace%get_ndf( )
    undf_wt = theta_inc_proxy%vspace%get_undf()

    mesh => theta_inc%get_mesh()

    if (theta_inc_proxy%is_dirty(depth=1)) CALL theta_inc_proxy%halo_exchange(depth=1)
    if (    delta_proxy%is_dirty(depth=1)) CALL     delta_proxy%halo_exchange(depth=1)

    if (theta_n_proxy%is_dirty(depth=1+1)) CALL theta_n_proxy%halo_exchange(depth=1+1)

    cross_stencil_wt => theta_n_proxy%vspace%get_stencil_dofmap(STENCIL_CROSS, t_stencil_size)
    cross_stencil_wt_map => cross_stencil_wt%get_whole_dofmap()
    cross_stencil_wt_size = cross_stencil_wt%get_size()
    
    map_t => theta_inc_proxy%vspace%get_whole_dofmap( )

    do cell = 1, mesh%get_last_edge_cell()

       call tracer_smagorinsky_diff_code(nlayers,                 &
                                  theta_inc_proxy%data,           &
                                  theta_n_proxy%data,             &
                                  cross_stencil_wt_size,          &
                                  cross_stencil_wt_map(:,:,cell), &
                                  delta_proxy%data,               &
                                  ndf_wt, undf_wt, map_t(:,cell) )
    end do
     
   call theta_inc_proxy%set_dirty()

  end subroutine invoke_tracer_smagorinsky_diff



end module psykal_lite_mod
