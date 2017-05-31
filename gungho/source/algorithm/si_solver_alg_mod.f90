!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!>@brief Routines for solving the semi-implicit equation set
module si_solver_alg_mod

  use, intrinsic :: ieee_arithmetic

  use constants_mod,             only: r_def, str_def, i_def
  use field_bundle_mod,          only: clone_bundle, &
                                       set_bundle_scalar, &
                                       bundle_axpy, &
                                       copy_bundle, &
                                       minus_bundle, &
                                       bundle_ax, &
                                       bundle_divide, &
                                       bundle_minmax, &
                                       bundle_inner_product
  use finite_element_config_mod, only: wtheta_on
  use runtime_constants_mod,     only: get_mass_matrix_diagonal, &
                                       w2_id, theta_space_id, eye_id
  use field_mod,                 only: field_type
  use formulation_config_mod,    only: eliminate_p
  use lhs_alg_mod,               only: lhs_alg
  use log_mod,                   only: log_event,         &
                                       log_scratch_space, &
                                       LOG_LEVEL_ERROR,   &
                                       LOG_LEVEL_DEBUG,   &
                                       LOG_LEVEL_TRACE,   &
                                       LOG_LEVEL_INFO
  use operator_mod,              only: operator_type
  use solver_config_mod,         only: maximum_iterations, &
                                       si_tolerance, &
                                       si_preconditioner, &
                                       si_postconditioner, &
                                       solver_si_preconditioner_none, &
                                       solver_si_preconditioner_diagonal, &
                                       solver_si_preconditioner_pressure, &
                                       solver_si_postconditioner_none, &
                                       solver_si_postconditioner_diagonal, &
                                       solver_si_postconditioner_pressure, &
                                       gcrk

  use timestepping_config_mod,   only: dt
  use derived_config_mod,        only: si_bundle_size, bundle_size
  use field_indices_mod,         only: igh_u, igh_t, igh_d, igh_p
  use output_config_mod,         only: subroutine_timers 
  use timer_mod,                 only: timer

  implicit none
  type(field_type), allocatable, private :: mm_diagonal(:)
  type(field_type), allocatable, private :: dx(:), Ax(:), residual(:), s(:), &
                                            w(:)
  type(field_type), allocatable, private :: v(:,:), Pv(:,:)
  type(field_type), allocatable, private :: x0_ext(:), rhs0_ext(:)
  real(r_def),      parameter,   private :: sc_err_min = 1.0e-16_r_def
private
  public  :: si_solver_alg
  public  :: si_solver_init
  private :: mixed_gmres_alg
 
contains
!>@brief Setup for the semi-implicit solver, extracts mass matrix diagonals and 
!!         sets up terms for the Newton-Krylov method if needed
!>@details Control routine to solve the system L*dx = rhs using a Krlyov method
!!         sets up terms needed either for a prescribed L
!>@param[inout] x0 The state array to solve for the increment of 
!>@param[in]    rhs0 Fixed rhs forcing for the solver
!>@param[in]    x_ref A reference state used for computing a proscribed L
  subroutine si_solver_alg(x0, rhs0, x_ref)
    use psykal_lite_mod,           only: invoke_set_field_scalar
    use psykal_lite_mod,           only: invoke_copy_field_data
    implicit none

    type(field_type), intent(inout)          :: x0(bundle_size)
    type(field_type), intent(in)             :: rhs0(bundle_size), &
                                                x_ref(bundle_size)
 
    real(kind=r_def)                         :: tau_dt ! tau_dt would eventually be set globally 
                                                       ! (probably the same place as alpha)
    integer(kind=i_def)                      :: i

    if ( subroutine_timers ) call timer('si_solver_alg')
    ! Set up tau_dt: to be used here and in subsequent algorithms
    tau_dt = 0.5_r_def*dt

! PSyclone built-ins support (v.1.3.1) fails for changing index in a loop, so the
! calls below are placeholders for when it becomes available
    if ( eliminate_p ) then
      call mixed_gmres_alg(x0, rhs0, x_ref, tau_dt)      
    else
      do i = 1,bundle_size
        call invoke_copy_field_data(rhs0(i), rhs0_ext(i))
        call invoke_copy_field_data(x0(i), x0_ext(i))
!         call invoke( copy_field(rhs0(i), rhs0_ext(i)), &
!                      copy_field(x0(i), x0_ext(i)) )
      end do
      ! Set initial guess to exner' and r_exner = 0
      call invoke_set_field_scalar(0.0_r_def, x0_ext(si_bundle_size))      
      call invoke_set_field_scalar(0.0_r_def, rhs0_ext(si_bundle_size)) 
!       call invoke( set_field_scalar(0.0_r_def, x0_ext(si_bundle_size)), &     
!                    set_field_scalar(0.0_r_def, rhs0_ext(si_bundle_size)) )     
      call mixed_gmres_alg(x0_ext, rhs0_ext, x_ref, tau_dt)
      do i = 1,bundle_size
        call invoke_copy_field_data(x0_ext(i), x0(i))
!         call invoke( copy_field(x0_ext(i), x0(i)) )
      end do
    end if

    if ( subroutine_timers ) call timer('si_solver_alg')

  end subroutine si_solver_alg
!=============================================================================!

  subroutine si_solver_init(x0)
    use lhs_alg_mod,                     only: lhs_init

    implicit none

    type(field_type),             intent(in) :: x0(bundle_size)
    integer                                  :: iter


    allocate( dx         (si_bundle_size), &
              Ax         (si_bundle_size), &
              residual   (si_bundle_size), &
              s          (si_bundle_size), &
              w          (si_bundle_size), &
              x0_ext     (si_bundle_size), &
              rhs0_ext   (si_bundle_size), &
              mm_diagonal(si_bundle_size), &
              v          (si_bundle_size,gcrk), &
              Pv         (si_bundle_size,gcrk) )
 

    mm_diagonal(igh_u) = get_mass_matrix_diagonal(w2_id)
    mm_diagonal(igh_t) = get_mass_matrix_diagonal(theta_space_id)
    mm_diagonal(igh_d) = get_mass_matrix_diagonal(eye_id)
    if ( .not. eliminate_p ) &
      mm_diagonal(igh_p) = get_mass_matrix_diagonal(eye_id)


   
    if ( eliminate_p ) then
      call clone_bundle(x0, x0_ext, si_bundle_size)
    else
      call clone_bundle(x0, x0_ext(1:bundle_size), bundle_size)
      x0_ext(si_bundle_size) = field_type( vector_space = x0(bundle_size)%get_function_space() )
    end if
    call clone_bundle(x0_ext, rhs0_ext, si_bundle_size)
    call clone_bundle(x0_ext, dx,       si_bundle_size)
    call clone_bundle(x0_ext, Ax,       si_bundle_size)
    call clone_bundle(x0_ext, s,        si_bundle_size)
    call clone_bundle(x0_ext, w,        si_bundle_size)
    call clone_bundle(x0_ext, residual, si_bundle_size)   
    do iter = 1,gcrk
      call clone_bundle(x0_ext,  v(:,iter), si_bundle_size)
      call clone_bundle(x0_ext, Pv(:,iter), si_bundle_size)
    end do
    ! Intitialise lhs fields
    call lhs_init(x0_ext)

  end subroutine si_solver_init

!=============================================================================!
!>@brief GMRES solver adapted for solving the semi-implicit equations
!>@details Standard GMRES algortihm from "Iterative methods for sparse linear
!! systems" by Y Saad, SIAM 2003
!>@param[inout] x0 State to increment 
!>@param[in]    rhs0 Fixed rhs so solve for
!>@param[in]    x_ref Reference state
!>@param[in]    tau_dt The offcentering parameter times the timestep

  subroutine mixed_gmres_alg(x0, rhs0, x_ref, tau_dt)

    implicit none

    type(field_type),             intent(inout) :: x0(si_bundle_size)
    type(field_type),             intent(in)    :: rhs0(si_bundle_size), &
                                                   x_ref(bundle_size)
    real(kind=r_def),             intent(in)    :: tau_dt

    ! The scalars
    real(kind=r_def)         :: h(gcrk+1, gcrk), u(gcrk), g(gcrk+1)
    real(kind=r_def)         :: beta,si, ci, nrm, h1, h2, p, q
    ! Others
    real(kind=r_def)               :: err, sc_err, init_err
    integer(kind=i_def)            :: iter, i, j, k, m
    integer(kind=i_def)            :: max_gmres_iter


    max_gmres_iter = maximum_iterations
    call rhs0(igh_u)%log_minmax(LOG_LEVEL_DEBUG,'max/min r_u = ')
    call rhs0(igh_t)%log_minmax(LOG_LEVEL_DEBUG,'max/min r_t = ')
    call rhs0(igh_d)%log_minmax(LOG_LEVEL_DEBUG,'max/min r_r = ')
    if (.not. eliminate_p) &
      call rhs0(igh_p)%log_minmax(LOG_LEVEL_DEBUG,'max/min r_p = ')

    err = bundle_inner_product(rhs0, rhs0, si_bundle_size)
    sc_err = max( sqrt(err), sc_err_min )
    init_err = sc_err

    if (err < si_tolerance) then
      write( log_scratch_space, '(A,I3,A,E12.4,A,E15.8)' ) &
           "gmres solver_algorithm:converged in ", 0,      &
           " iters, init=", init_err,                      &
           " final=", err
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    else
      write( log_scratch_space, '(A,I3,A,2E15.8)' )        &
             "solver_algorithm[", 0, "]: residual = ", init_err
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

    ! Initial guess
    call set_bundle_scalar(0.0_r_def, dx, si_bundle_size)

    call set_bundle_scalar(0.0_r_def, Ax, si_bundle_size)

    call minus_bundle( rhs0, Ax, residual, si_bundle_size )

    call bundle_preconditioner(s, residual, si_preconditioner, mm_diagonal, si_bundle_size )

    beta = sqrt(bundle_inner_product(s, s, si_bundle_size)) 

    call bundle_ax( 1.0_r_def/beta, s, v(:,1), si_bundle_size )

    h(:,:) = 0.0_r_def
    g(:)   = 0.0_r_def
    g(1)   = beta

    do iter = 1, max_gmres_iter

      do j = 1, gcrk

        call bundle_preconditioner(Pv(:,j), v(:,j), si_postconditioner, mm_diagonal, si_bundle_size)

        call lhs_alg(s, Pv(:,j), x_ref, tau_dt)
        call bundle_preconditioner(w, s, si_preconditioner, mm_diagonal, si_bundle_size )
        do k = 1, j
          h(k,j) =  bundle_inner_product( v(:,k), w, si_bundle_size )
          call bundle_axpy( -h(k,j), v(:,k), w, w, si_bundle_size )
        end do        
        h(j+1,j) = sqrt( bundle_inner_product( w, w, si_bundle_size ))
        if( j < gcrk ) then
          call bundle_ax(1.0_r_def/h(j+1,j), w, v(:,j+1), si_bundle_size)
        end if
      end do

      ! Solve (7.23) of Wesseling (see Saad's book)
      do m = 1, gcrk
        nrm    = sqrt( h(m,m)*h(m,m) + h(m+1,m)*h(m+1,m) )
        si     = h(m+1,m)/nrm
        ci     = h(m,m)/nrm
        p      = ci*g(m) + si*g(m+1)
        q      = -si*g(m) + ci*g(m+1)
        g(m)   = p
        g(m+1) = q
        do j = m, gcrk
          h1       = ci*h(m,j)   + si*h(m+1,j)
          h2       =-si*h(m,j)   + ci*h(m+1,j)
          h(m,j)   = h1
          h(m+1,j) = h2
        end do
      end do

      u(gcrk) = g(gcrk)/h(gcrk,gcrk)
      do i = gcrk-1, 1, -1
        u(i) = g(i)
        do j = i+1, gcrk
          u(i) = u(i) - h(i,j)*u(j)
        end do
        u(i) = u(i)/h(i,i)
      end do

      do i = 1, gcrk
        call bundle_axpy( u(i), Pv(:,i), dx, dx, si_bundle_size )
      end do

      ! Check for convergence
      call lhs_alg(Ax, dx, x_ref, tau_dt)

      call minus_bundle( rhs0, Ax, residual, si_bundle_size )

      beta = sqrt(bundle_inner_product(residual, residual, si_bundle_size))

      err = beta/sc_err
      write( log_scratch_space, '(A,I3,A,E15.8)' ) &
             "solver_algorithm[", iter, "]: residual = ", err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)

      if( err <  si_tolerance ) then
        write( log_scratch_space, '(A,I3,A,E12.4,A,E15.8)' )      &
             "GMRES solver_algorithm:converged in ", iter,        &
             " iters, init=", init_err, " final=", err
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        exit
      end if

      call bundle_preconditioner(s, residual, si_preconditioner, mm_diagonal, si_bundle_size)
      call bundle_ax(1.0_r_def/beta, s, v(:,1), si_bundle_size)

      g(:) = 0.0_r_def
      g(1) = beta

    end do

    if ((iter >= max_gmres_iter .and. err >  si_tolerance) &
        .or. ieee_is_nan(err)) then
      write( log_scratch_space, '(A, I3, A, E15.8)') &
           "GMRES solver_algorithm: NOT converged in ", max_gmres_iter, &
           " iters, Res=", err
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    ! Add increments to field
    call bundle_axpy(1.0_r_def, dx, x0, x0, si_bundle_size)  

  end subroutine mixed_gmres_alg
!=============================================================================!
!>@brief Applies a choosen preconditioner to a state x to produce state y
!>@param[inout] y Preconditioned state
!>@param[in]    x Original state
!>@param[in]    option choice of which preconditioner to use
!>@param[in]    mm Arrays containing diagonal approximation to mass matrices
!>@param[in]    si_bundle_size Number of fields the state arrays
  subroutine bundle_preconditioner(y, x, option, mm, si_bundle_size)
    use psykal_lite_mod,          only: invoke_copy_field_data
    use helmholtz_solver_alg_mod, only: helmholtz_solver_alg

    implicit none
    integer(kind=i_def), intent(in)    :: si_bundle_size
    type(field_type),    intent(inout) :: y(si_bundle_size)
    type(field_type),    intent(in)    :: x(si_bundle_size)
    type(field_type),    optional      :: mm(si_bundle_size)
    integer(kind=i_def), intent(in)    :: option
    integer(kind=i_def)                :: i

! PSyclone built-ins support (v.1.3.1) fails for changing index in a loop, so the
! callS below are placeholders for when it becomes available
    select case ( option)
      case( solver_si_preconditioner_none, solver_si_postconditioner_none )
        do i = 1,si_bundle_size
          call invoke_copy_field_data( x(i), y(i) )
!           call invoke( copy_field(x(i), y(i)) )
        end do
      case( solver_si_preconditioner_diagonal, solver_si_postconditioner_diagonal )
        do i = 1,si_bundle_size
          call invoke_copy_field_data( x(i), y(i) )
!           call invoke( copy_field(x(i), y(i)) )
        end do
        call bundle_divide(y, mm, si_bundle_size)
      case ( solver_si_preconditioner_pressure, solver_si_postconditioner_pressure )
        call set_bundle_scalar(0.0_r_def, y, bundle_size)
        call helmholtz_solver_alg(y, x)
      case default
        call log_event( 'Invalid si pre/postconditioner', LOG_LEVEL_ERROR )
    end select
  end subroutine bundle_preconditioner
!=============================================================================!
end module si_solver_alg_mod
