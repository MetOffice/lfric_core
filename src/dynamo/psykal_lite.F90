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

module psy

  use field_mod,      only : field_type, field_proxy_type 
  use operator_mod,   only : operator_type, operator_proxy_type
  use quadrature_mod, only : quadrature_type
  use constants_mod,  only : r_def

  implicit none
  public

contains

  !-------------------------------------------------------------------------------  
!> Invoke the solver kernel for a w3 field kernel
  subroutine invoke_w3_solver_kernel( lhs, rhs, chi, qr )

    use w3_solver_kernel_mod, only : solver_w3_code

    type( field_type ), intent( in ) :: lhs
    type( field_type ), intent( in ) :: rhs
    type( field_type ), intent( in ) :: chi(3)
    type( quadrature_type), intent(in) :: qr

    integer                    :: cell, nlayers, nqp_h, nqp_v
    integer                    :: ndf_w3, undf_w3, ndf_w0, undf_w0
    integer                    :: dim_w3, diff_dim_w0
    integer, pointer           :: map_w3(:), map_w0(:) => null() 
    real(kind=r_def), allocatable  :: w3_basis(:,:,:,:),               &
                                      w0_diff_basis(:,:,:,:)

    type( field_proxy_type )        :: lhs_proxy
    type( field_proxy_type )        :: rhs_proxy 
    type( field_proxy_type )        :: chi_proxy(3)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:)   => null(), &
                                 wv(:)   => null()

    lhs_proxy = lhs%get_proxy()
    rhs_proxy = rhs%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    nlayers = rhs_proxy%vspace%get_nlayers()
    ndf_w3  = rhs_proxy%vspace%get_ndf()
    undf_w3 = rhs_proxy%vspace%get_undf()
    dim_w3  = rhs_proxy%vspace%get_dim_space()
    allocate(w3_basis(dim_w3, ndf_w3, nqp_h, nqp_v) )

    ndf_w0 = chi_proxy(1)%vspace%get_ndf( )
    undf_w0 = chi_proxy(1)%vspace%get_undf( )
    diff_dim_w0 = chi_proxy(1)%vspace%get_dim_space_diff( )
    allocate(w0_diff_basis(diff_dim_w0, ndf_w0, nqp_h, nqp_v) )
    
    call rhs_proxy%vspace%compute_basis_function(              &
         w3_basis, ndf_w3, nqp_h, nqp_v, xp, zp )

    call chi_proxy(1)%vspace%compute_diff_basis_function(      &
         w0_diff_basis, ndf_w0, nqp_h, nqp_v, xp, zp )

     do cell = 1, lhs_proxy%vspace%get_ncell()
       map_w3 => lhs_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => chi_proxy(1)%vspace%get_cell_dofmap( cell )

       call solver_w3_code( nlayers,                           &
                            lhs_proxy%data,                    &
                            rhs_proxy%data,                    &
                            chi_proxy(1)%data,                 &
                            chi_proxy(2)%data,                 &
                            chi_proxy(3)%data,                 &
                            ndf_w3, undf_w3,                   &
                            map_w3,                            &
                            w3_basis,                          &
                            ndf_w0, undf_w0,                   &
                            map_w0,                            &
                            w0_diff_basis,                     &
                            nqp_h, nqp_v,                      &
                            wh, wv)
    end do

    deallocate(w3_basis, w0_diff_basis)
  end subroutine invoke_w3_solver_kernel

!-------------------------------------------------------------------------------  
  subroutine invoke_matrix_vector_mm(Ax,x, mm)
    use matrix_vector_mm_mod, only : matrix_vector_mm_code
    use enforce_bc_mod,          only: enforce_bc_w2                            
    use function_space_mod, only : W2 
    implicit none
    type(field_type), intent(in)    :: x
    type(field_type), intent(inout) :: Ax
    type(operator_type), intent(in) :: mm

    integer                 :: cell, nlayers, ndf, undf, fs
    integer, pointer        :: map(:), boundary_dofs(:,:) => null()

    type( field_proxy_type )        :: x_p
    type( field_proxy_type )        :: Ax_p
    type( operator_proxy_type )     :: mm_p

    x_p = x%get_proxy()
    Ax_p = Ax%get_proxy()
    mm_p = mm%get_proxy()

    nlayers = x_p%vspace%get_nlayers()
    ndf = x_p%vspace%get_ndf( )
    undf = x_p%vspace%get_undf( )
    fs = x%which_function_space()  
    if(fs  ==  W2) then 
       boundary_dofs => x_p%vspace%get_boundary_dofs()
    end if
    
    do cell = 1, x_p%vspace%get_ncell()
       map=>x_p%vspace%get_cell_dofmap(cell)
       call matrix_vector_mm_code( cell,               &
                                   nlayers,            &
                                   Ax_p%data,          &
                                   x_p%data,           &
                                   mm_p%ncell_3d,      &
                                   mm_p%local_stencil, &
                                   ndf,                &
                                   undf,               &
                                   map                 &
                                   )
       
       if(fs == W2) then ! this is yuk but haven't done others yet              
          call enforce_bc_w2(nlayers,ndf,undf,map,boundary_dofs,Ax_p%data)           
       end if    
    end do

    

  end subroutine invoke_matrix_vector_mm
  
!-------------------------------------------------------------------------------  
!> Invoke_initial_theta_kernel: Invoke the theta initialisation
  subroutine invoke_initial_theta_kernel( theta, chi )

    use initial_theta_kernel_mod, only : initial_theta_code

    type( field_type ), intent( in ) :: theta
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w0(:) => null()

    type( field_proxy_type )        :: theta_proxy
    type( field_proxy_type )        :: chi_proxy(3)

    theta_proxy  = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
     do cell = 1, theta_proxy%vspace%get_ncell()
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )
       call initial_theta_code( theta_proxy%vspace%get_nlayers(), &
                                theta_proxy%data,   &
                                chi_proxy(1)%data,  &
                                chi_proxy(2)%data,  &
                                chi_proxy(3)%data,  &
                                theta_proxy%vspace%get_ndf( ), &
                                theta_proxy%vspace%get_undf( ), &
                                map_w0 &
                               )
    end do 
  end subroutine invoke_initial_theta_kernel
  
   
!-------------------------------------------------------------------------------  
!> Invoke_rtheta_kernel: Invoke the RHS of the theta equation
  subroutine invoke_rtheta_kernel( r_theta, theta, f, rho, qr )

    use rtheta_kernel_mod, only : rtheta_code

    implicit none

    type( field_type ), intent( in ) :: r_theta, theta, f, rho
    type( quadrature_type), intent( in ) :: qr

    integer                 :: cell, nlayers, nqp_h, nqp_v
    integer                 :: ndf_w0, undf_w0, dim_w0, dim_diff_w0
    integer                 :: ndf_w2, undf_w2, dim_w2
    integer                 :: ndf_w3, undf_w3, dim_w3
    integer, pointer        :: map_w2(:), map_w0(:), orientation_w2(:), map_w3(:) => null()

    type( field_proxy_type )        :: r_theta_proxy, theta_proxy, &
                                       f_proxy, rho_proxy
    
    real(kind=r_def), allocatable  :: basis_w2(:,:,:,:), &
                                      basis_w0(:,:,:,:), &
                                      diff_basis_w0(:,:,:,:), &
                                      basis_w3(:,:,:,:)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:)   => null(), &
                                 wv(:)   => null()

    r_theta_proxy   = r_theta%get_proxy()
    theta_proxy     = theta%get_proxy()
    f_proxy         = f%get_proxy()
    rho_proxy       = rho%get_proxy()

    nlayers = r_theta_proxy%vspace%get_nlayers()

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_w2      = f_proxy%vspace%get_ndf( )
    dim_w2      = f_proxy%vspace%get_dim_space( )
    undf_w2     = f_proxy%vspace%get_undf()
    allocate(basis_w2(dim_w2,ndf_w2,nqp_h,nqp_v))

    ndf_w0      = r_theta_proxy%vspace%get_ndf( )
    dim_w0      = r_theta_proxy%vspace%get_dim_space( )
    dim_diff_w0 = r_theta_proxy%vspace%get_dim_space_diff( )
    undf_w0     = r_theta_proxy%vspace%get_undf()
    allocate(basis_w0(dim_w0,ndf_w0,nqp_h,nqp_v))
    allocate(diff_basis_w0(dim_diff_w0,ndf_w0,nqp_h,nqp_v))

    ndf_w3      = rho_proxy%vspace%get_ndf( )
    dim_w3      = rho_proxy%vspace%get_dim_space( )
    undf_w3     = rho_proxy%vspace%get_undf()
    allocate(basis_w3(dim_w3,ndf_w3,nqp_h,nqp_v))

    call f_proxy%vspace%compute_basis_function(                     &
        basis_w2, ndf_w2, nqp_h, nqp_v, xp, zp)    

    call r_theta_proxy%vspace%compute_basis_function(basis_w0, ndf_w0,      & 
         nqp_h, nqp_v, xp, zp)

    call r_theta_proxy%vspace%compute_diff_basis_function(                  &
         diff_basis_w0, ndf_w0, nqp_h, nqp_v, xp, zp)

    call rho_proxy%vspace%compute_basis_function(                     &
         basis_w3, ndf_w3, nqp_h, nqp_v, xp, zp)    

    do cell = 1, r_theta_proxy%vspace%get_ncell()
       map_w0 => r_theta_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => f_proxy%vspace%get_cell_dofmap( cell )
       map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )

       orientation_w2 => f_proxy%vspace%get_cell_orientation ( cell )

       call rtheta_code( nlayers, &
                         r_theta_proxy%data, &
                         theta_proxy%data, &
                         f_proxy%data, &
                         rho_proxy%data, &
                         ndf_w0, undf_w0, &
                         map_w0, &
                         basis_w0, &
                         diff_basis_w0, &
                         ndf_w2, undf_w2, &
                         map_w2, &
                         basis_w2, &
                         orientation_w2, &
                         ndf_w3, undf_w3, &
                         map_w3, &
                         basis_w3, &
                         nqp_h, nqp_v, wh, wv &
                         )
    end do 

    deallocate(basis_w2, basis_w0, diff_basis_w0, basis_w3)
  end subroutine invoke_rtheta_kernel 
  
!-------------------------------------------------------------------------------  
!> Invoke_ru_kernel: Invoke the RHS of the u equation
  subroutine invoke_ru_kernel( r_u )

    use ru_kernel_mod, only : ru_code

    type( field_type ), intent( in ) :: r_u

    integer                 :: cell, nlayers
    integer                 :: ndf_w2
    integer                 :: undf_w2
    integer, pointer        :: map_w2(:) => null()
    integer, pointer        :: boundary_dofs(:,:) => null()

    type( field_proxy_type )        :: r_u_proxy
    
    r_u_proxy  = r_u%get_proxy()
    
    boundary_dofs => r_u_proxy%vspace%get_boundary_dofs()
    
    nlayers = r_u_proxy%vspace%get_nlayers()

    ndf_w2      = r_u_proxy%vspace%get_ndf( )
    undf_w2     = r_u_proxy%vspace%get_undf()

    do cell = 1, r_u_proxy%vspace%get_ncell()

       map_w2 => r_u_proxy%vspace%get_cell_dofmap( cell )

       call ru_code( nlayers,                                              &
                     r_u_proxy%data,                                       &
                     ndf_w2, undf_w2,                                      &
                     map_w2,                                               &
                     boundary_dofs                                         &
                     )
    end do

  end subroutine invoke_ru_kernel
  
!-------------------------------------------------------------------------------  
!> Invoke_rrho_kernel: Invoke the RHS of the rho equation
  subroutine invoke_rrho_kernel( r_rho, u, qr )

    use rrho_kernel_mod, only : rrho_code

    type( field_type ), intent( in ) :: r_rho, u
    type( quadrature_type), intent( in ) :: qr

    integer                 :: cell, nlayers, nqp_v, nqp_h
    integer                 :: ndf_w2, ndf_w3
    integer                 :: undf_w2, undf_w3
    integer                 :: diff_dim_w2, dim_w3
    integer, pointer  :: map_w3(:), map_w2(:), orientation_w2(:) => null()

    type( field_proxy_type )        :: r_rho_proxy, u_proxy
    
    real(kind=r_def), allocatable  :: basis_w3(:,:,:,:), &
                                      diff_basis_w2(:,:,:,:) 

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()


    r_rho_proxy  = r_rho%get_proxy()
    u_proxy      = u%get_proxy()
    
    nlayers = r_rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_w3  = r_rho_proxy%vspace%get_ndf( )
    dim_w3  = r_rho_proxy%vspace%get_dim_space( )
    undf_w3 = r_rho_proxy%vspace%get_undf()
    allocate(basis_w3(dim_w3,ndf_w3,nqp_h,nqp_v))

    ndf_w2      = u_proxy%vspace%get_ndf( )
    diff_dim_w2 = u_proxy%vspace%get_dim_space_diff( )
    undf_w2     = u_proxy%vspace%get_undf()
    allocate(diff_basis_w2(diff_dim_w2,ndf_w2,nqp_h,nqp_v))

    call r_rho_proxy%vspace%compute_basis_function(basis_w3, ndf_w3,     & 
                                                   nqp_h, nqp_v, xp, zp)

    call u_proxy%vspace%compute_diff_basis_function(                     &
         diff_basis_w2, ndf_w2, nqp_h, nqp_v, xp, zp)

    do cell = 1, r_rho_proxy%vspace%get_ncell()

       map_w3 => r_rho_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => u_proxy%vspace%get_cell_dofmap( cell )
       orientation_w2 => u_proxy%vspace%get_cell_orientation( cell )

       call rrho_code( nlayers,                          &
                       r_rho_proxy%data,                 &
                       u_proxy%data,                     &
                       ndf_w3,                           &
                       undf_w3,                          &
                       map_w3,                           &
                       basis_w3,                         &
                       ndf_w2,                           &
                       undf_w2,                          &
                       map_w2,                           &
                       diff_basis_w2,                    &
                       orientation_w2,                   &
                       nqp_h,                            &
                       nqp_v,                            &
                       wh,                               &
                       wv                                &
                       )
    end do 
  
    deallocate(basis_w3, diff_basis_w2 )

  end subroutine invoke_rrho_kernel

!-------------------------------------------------------------------------------    
!> invoke_compute_mass_matrix_w0: Calculate mass matrix for W0 space
  subroutine invoke_compute_mass_matrix_w0(mm, chi, qr)
    use compute_mass_matrix_kernel_w0_mod, only :  compute_mass_matrix_w0_code
    implicit none

    type( operator_type ),  intent(inout) :: mm
    type( field_type ),     intent(inout) :: chi(3)
    type( quadrature_type), intent(in)    :: qr

    integer :: cell, nlayers, ndf, undf, ndf_chi, dim, diff_dim, dim_diff_chi
    integer :: nqp_h, nqp_v

    integer, pointer        :: map(:)  => null()
    real ( kind=r_def ), allocatable :: diff_basis(:,:,:,:),   &
                                        basis(:,:,:,:),        &
                                        diff_basis_chi(:,:,:,:)

    type( operator_proxy_type) :: mm_proxy
    type( field_proxy_type)    :: chi_proxy(3)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    chi_proxy(1) = chi(1)%get_proxy()                                           
    chi_proxy(2) = chi(2)%get_proxy()                                           
    chi_proxy(3) = chi(3)%get_proxy()                                           
    mm_proxy = mm%get_proxy()

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    dim_diff_chi = chi_proxy(1)%vspace%get_dim_space_diff( )

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    nlayers  = mm_proxy%fs_from%get_nlayers()

    ndf      = chi_proxy(1)%vspace%get_ndf()
    undf     = chi_proxy(1)%vspace%get_undf()
    diff_dim = chi_proxy(1)%vspace%get_dim_space_diff( )
    allocate( diff_basis(diff_dim,ndf,nqp_h,nqp_v),               &
              diff_basis_chi(dim_diff_chi, ndf_chi, nqp_h, nqp_v) )

    dim = mm_proxy%fs_from%get_dim_space_diff( )
    allocate(basis(dim,ndf,nqp_h,nqp_v) )

    call chi_proxy(1)%vspace%compute_diff_basis_function(                   &
         diff_basis,ndf, nqp_h, nqp_v, xp, zp )

    call mm_proxy%fs_from%compute_basis_function(                                 &
         basis, ndf, nqp_h, nqp_v, xp, zp ) 

    call chi_proxy(1)%vspace%compute_diff_basis_function( &
           diff_basis_chi, ndf_chi, nqp_h, nqp_v, xp, zp)

    do cell = 1, chi_proxy(1)%vspace%get_ncell()                                
       map => chi_proxy(1)%vspace%get_cell_dofmap(cell)                     

       call compute_mass_matrix_w0_code( cell,                              &
                                         nlayers,                           &
                                         mm_proxy%ncell_3d,                 &
                                         mm_proxy%local_stencil,            &
                                         chi_proxy(1)%data,                 &
                                         chi_proxy(2)%data,                 &
                                         chi_proxy(3)%data,                 &
                                         ndf,                               &
                                         undf,                              &
                                         map,                               &
                                         basis,                             &
                                         diff_basis_chi,                    &
                                         nqp_h,                             &
                                         nqp_v,                             &
                                         wh,                                &
                                         wv                                 &
                                         )
    end do
    
    deallocate(basis, diff_basis)
    
  end subroutine invoke_compute_mass_matrix_w0

  subroutine invoke_compute_mass_matrix_w1(mm, chi, qr)
    use compute_mass_matrix_kernel_w1_mod, only :  compute_mass_matrix_w1_code
    implicit none

    type( operator_type ),  intent(inout) :: mm
    type( field_type ),     intent(inout) :: chi(3)
    type( quadrature_type), intent(in)    :: qr

    integer :: cell, nlayers, ndf, dim, diff_dim
    integer :: ndf_chi, undf_chi
    integer :: nqp_h, nqp_v

    integer, pointer        ::  map_chi(:)  => null()
    integer, pointer        ::  orientation_w1(:) => null()
    real ( kind=r_def ), allocatable :: diff_basis(:,:,:,:), &
                                        basis(:,:,:,:)

    type( operator_proxy_type) :: mm_proxy
    type( field_proxy_type)    :: chi_proxy(3)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    chi_proxy(1) = chi(1)%get_proxy()                                           
    chi_proxy(2) = chi(2)%get_proxy()                                           
    chi_proxy(3) = chi(3)%get_proxy()                                           
    mm_proxy = mm%get_proxy()  
    
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_chi      = chi_proxy(1)%vspace%get_ndf()
    undf_chi     = chi_proxy(1)%vspace%get_undf()
    diff_dim = chi_proxy(1)%vspace%get_dim_space_diff( )
    allocate(diff_basis(diff_dim,ndf_chi,nqp_h,nqp_v) )

    nlayers  = mm_proxy%fs_from%get_nlayers()
    ndf  = mm_proxy%fs_from%get_ndf()
    dim = mm_proxy%fs_from%get_dim_space( )
    allocate(basis(dim,ndf,nqp_h,nqp_v) )

    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         diff_basis,ndf_chi, nqp_h, nqp_v, xp, zp )

    call mm_proxy%fs_from%compute_basis_function( &
         basis, ndf, nqp_h, nqp_v, xp, zp ) 

    do cell = 1, chi_proxy(1)%vspace%get_ncell()
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)
       orientation_w1 => mm_proxy%fs_from%get_cell_orientation ( cell )

       call compute_mass_matrix_w1_code( cell,                              &
                                         nlayers,                           &
                                         mm_proxy%ncell_3d,                 &
                                         mm_proxy%local_stencil,            &
                                         chi_proxy(1)%data,                 &
                                         chi_proxy(2)%data,                 &
                                         chi_proxy(3)%data,                 &
                                         ndf,                               &
                                         basis,                             &
                                         orientation_w1,                    &
                                         ndf_chi,                           &
                                         undf_chi,                          &
                                         map_chi,                           &
                                         diff_basis,                        &
                                         nqp_h,                             &
                                         nqp_v,                             &
                                         wh,                                &
                                         wv                                 &
                                         )
    end do

    deallocate(basis, diff_basis)
    
  end subroutine invoke_compute_mass_matrix_w1

 subroutine invoke_compute_mass_matrix_w2(mm, chi, qr)
    use compute_mass_matrix_kernel_w2_mod, only :  compute_mass_matrix_w2_code
    implicit none

    type( operator_type ),  intent(inout) :: mm
    type( field_type ),     intent(inout) :: chi(3)
    type( quadrature_type), intent(in)    :: qr

    integer :: cell, nlayers, ndf, dim, diff_dim
    integer :: ndf_chi, undf_chi
    integer :: nqp_h, nqp_v

    integer, pointer        ::  map_chi(:)  => null()
    integer, pointer        :: orientation_w2(:) => null()
    real ( kind=r_def ), allocatable :: diff_basis(:,:,:,:), &
                                        basis(:,:,:,:)

    type( operator_proxy_type) :: mm_proxy
    type( field_proxy_type)    :: chi_proxy(3)

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    mm_proxy = mm%get_proxy()

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_chi      = chi_proxy(1)%vspace%get_ndf()
    undf_chi     = chi_proxy(1)%vspace%get_undf()
    diff_dim = chi_proxy(1)%vspace%get_dim_space_diff( )
    allocate(diff_basis(diff_dim,ndf_chi,nqp_h,nqp_v) )

    nlayers  = mm_proxy%fs_from%get_nlayers()
    ndf  = mm_proxy%fs_from%get_ndf()
    dim = mm_proxy%fs_from%get_dim_space( )
    allocate(basis(dim,ndf,nqp_h,nqp_v) )

    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         diff_basis,ndf_chi, nqp_h, nqp_v, xp, zp )

    call mm_proxy%fs_from%compute_basis_function( &
         basis, ndf, nqp_h, nqp_v, xp, zp ) 

    do cell = 1, chi_proxy(1)%vspace%get_ncell()
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap(cell)

       orientation_w2 => mm_proxy%fs_from%get_cell_orientation( cell )

       call compute_mass_matrix_w2_code( cell,                              &
                                         nlayers,                           &
                                         mm_proxy%ncell_3d,                 &
                                         mm_proxy%local_stencil,            &
                                         chi_proxy(1)%data,                 &
                                         chi_proxy(2)%data,                 &
                                         chi_proxy(3)%data,                 &
                                         ndf,                               &
                                         basis,                             &
                                         orientation_w2,                    &
                                         ndf_chi,                           &
                                         undf_chi,                          &
                                         map_chi,                           &
                                         diff_basis,                        &
                                         nqp_h,                             &
                                         nqp_v,                             &
                                         wh,                                &
                                         wv                                 &
                                         )
    end do

    deallocate(basis, diff_basis)
    
  end subroutine invoke_compute_mass_matrix_w2

!-------------------------------------------------------------------------------    
!> invoke_inner_prod: Calculate inner product of x and y
  subroutine invoke_inner_prod(x,y,inner_prod)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x,y
    real(kind=r_def),    intent(out) :: inner_prod

    type( field_proxy_type)          ::  x_p,y_p
    integer                          :: i,undf

    x_p = x%get_proxy()
    y_p = y%get_proxy()

    undf = x_p%vspace%get_undf()
    !sanity check
    if(undf /= y_p%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inner_prod:x and y live on different w-spaces",LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    inner_prod = 0.0_r_def
    do i = 1,undf
      inner_prod = inner_prod + ( x_p%data(i) * y_p%data(i) )
    end do

  end subroutine invoke_inner_prod
  
!-------------------------------------------------------------------------------   
!> invoke_axpy:  (a * x + y) ; a-scalar, x,y-vector     
  subroutine invoke_axpy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) + field2_proxy%data(i)
    end do
  end subroutine invoke_axpy
  
!-------------------------------------------------------------------------------   
!> invoke_axmy:  (a * x - y) ; a-scalar, x,y-vector
  subroutine invoke_axmy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) - field2_proxy%data(i)
    end do
  end subroutine invoke_axmy
  
!-------------------------------------------------------------------------------   
!> invoke_copy_field_data: copy the data from one field to another ( a = b )
  subroutine invoke_copy_field_data(field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)
    end do
  end subroutine invoke_copy_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_minus_field_data: Subtract values of field2 from values of field 
!> ( c = a - b )
  subroutine invoke_minus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) - field2_proxy%data(i)
    end do
  end subroutine invoke_minus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_plus_field_data:  Add values of field2 to values of field1
!> ( c = a + b )
  subroutine invoke_plus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) + field2_proxy%data(i)
    end do
  end subroutine invoke_plus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_set_field_scalar: set all values in a field to a single value
  subroutine invoke_set_field_scalar(scalar, field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field_res_proxy
    integer                            :: i,undf

    field_res_proxy = field_res%get_proxy()

    undf = field_res_proxy%vspace%get_undf()

    do i = 1,undf
      field_res_proxy%data(i) = scalar
    end do
  end subroutine invoke_set_field_scalar
!-------------------------------------------------------------------------------   
!> invoke_divide_field: divide the values of field1 by field2 and put result in
!>field_res
!> c = a/b
  subroutine invoke_divide_field(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)/field2_proxy%data(i)
    end do
  end subroutine invoke_divide_field

!-------------------------------------------------------------------------------   
!> Invoke_gp_rhs: Invoke the scalar RHS for a Galerkin projection
  subroutine invoke_gp_rhs( rhs, field, chi, qr )

    use gp_rhs_kernel_mod, only : gp_rhs_code

    implicit none

    type( field_type ),  intent( in ) :: rhs, field
    type( field_type ),  intent( in ) :: chi(3)
    type( quadrature_type), intent( in ) :: qr

    type( field_proxy_type)           :: rhs_proxy, field_proxy
    type( field_proxy_type)           :: chi_proxy(3)

    integer          :: cell, nlayers, ndf, undf, dim
    integer          :: ndf_f, undf_f, dim_f,ndf_chi, undf_chi, dim_chi
    integer          :: nqp_h, nqp_v
    integer, pointer :: map(:), map_chi(:), map_f(:) => null()

    real(kind=r_def), allocatable  :: basis(:,:,:,:), f_basis(:,:,:,:), &
                                  chi_diff_basis(:,:,:,:) 

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    rhs_proxy    = rhs%get_proxy()
    field_proxy  = field%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ! Unpack data - get the integers 
    nlayers = rhs_proxy%vspace%get_nlayers( )
    ndf     = rhs_proxy%vspace%get_ndf( )
    undf    = rhs_proxy%vspace%get_undf( )
    dim     = rhs_proxy%vspace%get_dim_space( )        

    ndf_f   = field_proxy%vspace%get_ndf( )
    undf_f  = field_proxy%vspace%get_undf( )
    dim_f   = field_proxy%vspace%get_dim_space( )        

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space_diff( )        

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()    
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    allocate( basis(dim, ndf, nqp_h, nqp_v),       &
              f_basis(dim_f, ndf_f, nqp_h, nqp_v), &
              chi_diff_basis(dim_chi, ndf_chi, nqp_h, nqp_v) )

    call rhs_proxy%vspace%compute_basis_function( &
         basis, ndf, nqp_h, nqp_v, xp, zp)

    call field_proxy%vspace%compute_basis_function( &
         f_basis, ndf_f, nqp_h, nqp_v, xp, zp)

    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         chi_diff_basis, ndf_chi, nqp_h, nqp_v, xp, zp)

    do cell = 1, rhs_proxy%vspace%get_ncell()
       map     => rhs_proxy   %vspace%get_cell_dofmap( cell )
       map_f   => field_proxy %vspace%get_cell_dofmap( cell ) 
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

       call gp_rhs_code( nlayers,                     &
                         rhs_proxy%data, &
                         field_proxy%data, &
                         chi_proxy(1)%data, &
                         chi_proxy(2)%data, &
                         chi_proxy(3)%data , &
                         ndf, undf,  &
                         map, basis, &
                         ndf_f, undf_f, &
                         map_f, f_basis, &
                         ndf_chi, undf_chi, &
                         map_chi, chi_diff_basis, &
                         nqp_h, nqp_v, wh, wv )
    end do

    deallocate(basis, f_basis, chi_diff_basis)

  end subroutine invoke_gp_rhs
  
!-------------------------------------------------------------------------------  
!> Invoke_gp_vector_rhs: Invoke the vector RHS for a Galerkin projection
  subroutine invoke_gp_vector_rhs( rhs, field, chi, qr )

    use gp_vector_rhs_kernel_mod, only : gp_vector_rhs_code

    implicit none

    type( field_type ),  intent( in ) :: rhs(3), field
    type( field_type ),  intent( in ) :: chi(3)
    type( quadrature_type), intent( in ) :: qr

    type( field_proxy_type)           :: rhs_proxy(3), field_proxy
    type( field_proxy_type)           :: chi_proxy(3)

    integer          :: cell, dir, nlayers, nqp_h, nqp_v
    integer          :: ndf, undf, dim
    integer          :: ndf_f, undf_f, dim_f
    integer          :: ndf_chi, undf_chi, dim_chi, dim_diff_chi
    integer, pointer :: map(:), map_chi(:), map_f(:) => null()
    integer, pointer :: orientation(:) => null()

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    real(kind=r_def), allocatable  :: basis(:,:,:,:), f_basis(:,:,:,:), &
                                  chi_basis(:,:,:,:), chi_diff_basis(:,:,:,:)

    do dir = 1,3
      rhs_proxy(dir) = rhs(dir)%get_proxy()
      chi_proxy(dir) = chi(dir)%get_proxy()
    end do
    field_proxy  = field%get_proxy()

    ! Unpack data 
    nlayers = rhs_proxy(1)%vspace%get_nlayers( )
    ndf     = rhs_proxy(1)%vspace%get_ndf( )
    undf    = rhs_proxy(1)%vspace%get_undf( )
    dim     = rhs_proxy(1)%vspace%get_dim_space( ) 

    ndf_f   = field_proxy%vspace%get_ndf( )
    undf_f  = field_proxy%vspace%get_undf( )
    dim_f   = field_proxy%vspace%get_dim_space( )        

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf( )
    dim_chi  = chi_proxy(1)%vspace%get_dim_space( ) 
    dim_diff_chi = chi_proxy(1)%vspace%get_dim_space_diff( ) 

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    allocate( basis(dim, ndf, nqp_h, nqp_v),           &
              f_basis(dim_f, ndf_f, nqp_h, nqp_v), &
              chi_basis(dim_chi, ndf_chi, nqp_h, nqp_v), &
              chi_diff_basis(dim_diff_chi, ndf_chi, nqp_h, nqp_v) )         
    
    call rhs_proxy(1)%vspace%compute_basis_function( &
         basis, ndf, nqp_h, nqp_v, xp, zp)

    call field_proxy%vspace%compute_basis_function( &
         f_basis, ndf_f, nqp_h, nqp_v, xp, zp)

    call chi_proxy(1)%vspace%compute_basis_function( &
         chi_basis, ndf_chi, nqp_h, nqp_v, xp, zp)    

    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         chi_diff_basis, ndf_chi, nqp_h, nqp_v, xp, zp)    

    do cell = 1, rhs_proxy(1)%vspace%get_ncell()
       map     => rhs_proxy(1)%vspace%get_cell_dofmap( cell )
       map_f   => field_proxy %vspace%get_cell_dofmap( cell ) 
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       orientation => field_proxy%vspace%get_cell_orientation ( cell )

       call gp_vector_rhs_code( nlayers, &
                                rhs_proxy(1)%data, &
                                rhs_proxy(2)%data, &
                                rhs_proxy(3)%data, &
                                field_proxy%data, &
                                chi_proxy(1)%data, &
                                chi_proxy(2)%data, &
                                chi_proxy(3)%data, &
                                ndf, undf, &
                                map, basis, &
                                ndf_f, undf_f, &
                                map_f, f_basis, &
                                ndf_chi, undf_chi, &
                                map_chi, chi_basis, &
                                chi_diff_basis, &
                                orientation, &
                                nqp_h, nqp_v, wh, wv )
    end do

  end subroutine invoke_gp_vector_rhs

!-------------------------------------------------------------------------------   
!> invoke_copy_scaled_field_data: copy the scaled data from one field to another ( a = scaler*b )
  subroutine invoke_copy_scaled_field_data(scaler,field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    real( kind=r_def ), intent(in)     :: scaler
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      field_res_proxy%data(i) = scaler*field1_proxy%data(i)
    end do
  end subroutine invoke_copy_scaled_field_data

!-------------------------------------------------------------------------------  
!> Invoke_flux_rhs_kernel: Invoke the RHS of the flux equation Flux = u*f
  subroutine invoke_flux_rhs( rhs, u, f, chi, qr )

    use flux_rhs_kernel_mod, only : flux_rhs_code

    type( field_type ),     intent( in ) :: rhs, f, u
    type( field_type ),     intent( in ) :: chi(3) 
    type( quadrature_type), intent( in ) :: qr

    integer          :: cell
    integer          :: ndf_f, undf_f, dim_f, &
                        ndf_u, undf_u, dim_u, &
                        ndf_chi, undf_chi, dim_diff_chi, &
                        nqp_h, nqp_v

    integer, pointer :: map_f(:) => null(), &
                        map_u(:) => null(), & 
                        map_chi(:) => null(), &
                        boundary_dofs(:,:) => null(), &
                        orientation_u(:) => null()

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    type( field_proxy_type )        :: rhs_proxy, f_proxy, u_proxy
    type( field_proxy_type )        :: chi_proxy(3) 
    
    real(kind=r_def), dimension(:,:,:,:), allocatable :: &
                                  basis_u, &
                                  basis_f, &
                                  diff_basis_chi

    rhs_proxy    = rhs%get_proxy()
    f_proxy      = f%get_proxy()
    u_proxy      = u%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    boundary_dofs => rhs_proxy%vspace%get_boundary_dofs()

    ndf_u  = rhs_proxy%vspace%get_ndf( )
    undf_u = rhs_proxy%vspace%get_undf( )
    dim_u  = rhs_proxy%vspace%get_dim_space( ) 

    ndf_f  = f_proxy%vspace%get_ndf( )
    undf_f = f_proxy%vspace%get_undf( )
    dim_f  = f_proxy%vspace%get_dim_space( ) 

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf( )
    dim_diff_chi = chi_proxy(1)%vspace%get_dim_space_diff( ) 

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    allocate( basis_u(dim_u, ndf_u, nqp_h, nqp_v),           &
              basis_f(dim_f, ndf_f, nqp_h, nqp_v),           &
              diff_basis_chi(dim_diff_chi, ndf_chi, nqp_h, nqp_v) )         
    
    call rhs_proxy%vspace%compute_basis_function( &
         basis_u, ndf_u, nqp_h, nqp_v, xp, zp)  
    call f_proxy%vspace%compute_basis_function( &
         basis_f, ndf_f, nqp_h, nqp_v, xp, zp)  
    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         diff_basis_chi, ndf_chi, nqp_h, nqp_v, xp, zp)  
    
    do cell = 1, rhs_proxy%vspace%get_ncell()
       map_f   => f_proxy%vspace%get_cell_dofmap( cell )
       map_u   => rhs_proxy%vspace%get_cell_dofmap( cell )
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       orientation_u => u_proxy%vspace%get_cell_orientation ( cell )

      call flux_rhs_code( rhs_proxy%vspace%get_nlayers(), &
                           ndf_u, &
                           undf_u, &
                           map_u, &
                           basis_u, &
                           boundary_dofs, &
                           orientation_u, &
                           rhs_proxy%data, &
                           u_proxy%data, &
                           ndf_f, &
                           undf_f, &
                           map_f, &
                           basis_f, &                             
                           f_proxy%data, &
                           ndf_chi, &
                           undf_chi, &
                           map_chi, &
                           diff_basis_chi, &   
                           chi_proxy(1)%data, &
                           chi_proxy(2)%data, &
                           chi_proxy(3)%data, &
                           nqp_h, &
                           nqp_v, &
                           wh, &
                           wv &
                           )           
    end do 
  end subroutine invoke_flux_rhs
 
!-------------------------------------------------------------------------------  
!> Invoke_ru_kernel: Invoke the RHS of the u equation
  subroutine invoke_linear_ru_kernel( r_u, u, rho, theta, phi, chi, qr )

    use linear_ru_kernel_mod, only : linear_ru_code

    type( field_type ), intent( in ) :: r_u, u, rho, theta, phi
    type( field_type ), intent( in ) :: chi(3)
    type( quadrature_type), intent( in ) :: qr

    integer                 :: cell, nlayers, nqp_h, nqp_v
    integer                 :: ndf_w0, ndf_w2, ndf_w3
    integer                 :: undf_w0, undf_w2, undf_w3
    integer                 :: dim_w0, diff_dim_w0, dim_w2, diff_dim_w2,dim_w3
    integer, pointer        :: map_w3(:), map_w2(:), map_w0(:), orientation_w2(:) => null()
    integer, pointer        :: boundary_dofs(:,:) => null()

    type( field_proxy_type )        :: r_u_proxy, u_proxy, rho_proxy, theta_proxy, phi_proxy
    type( field_proxy_type )        :: chi_proxy(3)
    
    real(kind=r_def), allocatable  :: basis_w3(:,:,:,:), &
                                      basis_w2(:,:,:,:), &
                                      basis_w0(:,:,:,:), &
                                      diff_basis_w0(:,:,:,:), &
                                      diff_basis_w2(:,:,:,:) 

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    r_u_proxy    = r_u%get_proxy()
    u_proxy      = u%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()
    phi_proxy    = phi%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    boundary_dofs => r_u_proxy%vspace%get_boundary_dofs()

    nlayers = rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_w3  = rho_proxy%vspace%get_ndf( )
    dim_w3  = rho_proxy%vspace%get_dim_space( )
    undf_w3 = rho_proxy%vspace%get_undf()
    allocate(basis_w3(dim_w3,ndf_w3,nqp_h,nqp_v))

    ndf_w2      = r_u_proxy%vspace%get_ndf( )
    dim_w2      = r_u_proxy%vspace%get_dim_space( )
    diff_dim_w2 = r_u_proxy%vspace%get_dim_space_diff( )
    undf_w2     = r_u_proxy%vspace%get_undf()
    allocate(basis_w2(dim_w2,ndf_w2,nqp_h,nqp_v))
    allocate(diff_basis_w2(diff_dim_w2,ndf_w2,nqp_h,nqp_v))

    ndf_w0      = theta_proxy%vspace%get_ndf( )
    dim_w0      = theta_proxy%vspace%get_dim_space( )
    diff_dim_w0 = theta_proxy%vspace%get_dim_space_diff( )
    undf_w0     = theta_proxy%vspace%get_undf()
    allocate(basis_w0(dim_w0,ndf_w0,nqp_h,nqp_v))
    allocate(diff_basis_w0(diff_dim_w0,ndf_w0,nqp_h,nqp_v))

    call rho_proxy%vspace%compute_basis_function(basis_w3, ndf_w3,         & 
                                                   nqp_h, nqp_v, xp, zp)    

    call r_u_proxy%vspace%compute_basis_function(basis_w2, ndf_w2,         & 
                                                   nqp_h, nqp_v, xp, zp)    

    call r_u_proxy%vspace%compute_diff_basis_function(                     &
         diff_basis_w2, ndf_w2, nqp_h, nqp_v, xp, zp)

    call theta_proxy%vspace%compute_basis_function(basis_w0, ndf_w0,      & 
                                                   nqp_h, nqp_v, xp, zp)    

    call theta_proxy%vspace%compute_diff_basis_function(                  &
         diff_basis_w0, ndf_w0, nqp_h, nqp_v, xp, zp)


    
    do cell = 1, r_u_proxy%vspace%get_ncell()

       map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => r_u_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )

       orientation_w2 => r_u_proxy%vspace%get_cell_orientation ( cell )

       call linear_ru_code( nlayers,                                      &
                            ndf_w2, undf_w2,                              &
                            map_w2, basis_w2, diff_basis_w2,              &
                            boundary_dofs,                                &
                            orientation_w2,                               &
                            r_u_proxy%data,                               &
                            u_proxy%data,                                 &
                            ndf_w3, undf_w3,                              &
                            map_w3, basis_w3,                             &
                            rho_proxy%data,                               &
                            ndf_w0, undf_w0,                              &
                            map_w0, basis_w0, diff_basis_w0,              &   
                            theta_proxy%data,                             &
                            phi_proxy%data,                               &
                            chi_proxy(1)%data,                            &
                            chi_proxy(2)%data,                            &
                            chi_proxy(3)%data,                            &
                            nqp_h, nqp_v, wh, wv                          &
                            )           
    end do

    deallocate(basis_w3, basis_w2, diff_basis_w2, basis_w0, diff_basis_w0)
    
  end subroutine invoke_linear_ru_kernel
  
!> Invoke_compute_geopotential_kernel: Invoke the computation of the
!! geopotential
  subroutine invoke_compute_geopotential_kernel( phi, chi )

    use compute_geopotential_kernel_mod, only : compute_geopotential_code

    type( field_type ), intent( in ) :: phi 
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w0(:) => null()

    type( field_proxy_type )        :: phi_proxy
    type( field_proxy_type )        :: chi_proxy(3)

    phi_proxy    = phi%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
     do cell = 1, phi_proxy%vspace%get_ncell()
       map_w0 => phi_proxy%vspace%get_cell_dofmap( cell )
       call compute_geopotential_code( phi_proxy%vspace%get_nlayers(), &
                                       phi_proxy%data, &
                                       chi_proxy(1)%data, &
                                       chi_proxy(2)%data, &
                                       chi_proxy(3)%data, &
                                       phi_proxy%vspace%get_ndf( ), &
                                       phi_proxy%vspace%get_undf( ), &
                                       map_w0                        &
                                     )
    end do 
  end subroutine invoke_compute_geopotential_kernel
  
!-------------------------------------------------------------------------------    
!> invoke_sum_field: Sum all values of a field x
  subroutine invoke_sum_field( x, field_sum )
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x
    real(kind=r_def),    intent(out) :: field_sum

    type( field_proxy_type)          :: x_p
    integer                          :: df, undf

    x_p = x%get_proxy()   

    undf = x_p%vspace%get_undf()
    
    field_sum = 0.0_r_def
    do df = 1,undf
      field_sum = field_sum + x_p%data(df)
    end do

  end subroutine invoke_sum_field
!-------------------------------------------------------------------------------   
!> invoke_axpby:  z = (a * x + b * y) ; a,b-scalar, x,y-vector     
  subroutine invoke_axpby(a,x,b,y,z)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: x, y
    type( field_type ), intent(inout ) :: z
    real(kind=r_def),   intent(in )    :: a, b
    type( field_proxy_type)            :: x_proxy,y_proxy      &
                                        , z_proxy
    integer                            :: i,undf

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

    do i = 1,undf
      z_proxy%data(i) = (a * x_proxy%data(i)) + (b * y_proxy%data(i))
    end do
  end subroutine invoke_axpby
!-------------------------------------------------------------------------------   
!> invoke_multiply_field: compute y = a*x for scalar a and fields y and x
  subroutine invoke_multiply_field(a, x, y)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in)    :: x
    type( field_type ), intent(inout) :: y
    real(kind=r_def),   intent(in)    :: a
    type( field_proxy_type)           :: x_proxy, y_proxy
    integer                           :: i,undf

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

    do i = 1,undf
      y_proxy%data(i) = a*x_proxy%data(i)
    end do
  end subroutine invoke_multiply_field

  end module psy
