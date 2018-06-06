!-------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------

!> @brief Abstract base class for iterative solver
!>
!> @detail This class can be used as a base class for iterative solvers
!> which solve the system \f$Ax=b\f$. It contains a linear operator and a
!> preconditioner object which are required by each Krylov-subspace solver.
!> The type also defines an interface for the solver application. Common
!> data such as the relative and absolute solver tolerance are stored as
!> data members.

module iterative_solver_mod
  use constants_mod,        only : r_def, i_def
  use vector_mod,           only : abstract_vector_type
  use linear_operator_mod,  only : abstract_linear_operator_type
  use preconditioner_mod,   only : abstract_preconditioner_type
  use log_mod,              only : log_event, LOG_LEVEL_INFO, &
                                   LOG_LEVEL_ERROR, &
                                   log_scratch_space
  
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

  implicit none
  private

  !> @brief abstract solver type for the solver API
  type, public, abstract :: abstract_iterative_solver_type
     private
     !> Linear operator
     class(abstract_preconditioner_type),  pointer :: prec => null()
     class(abstract_linear_operator_type), pointer :: lin_op => null()
      !> preconditioner

      ! relative tolerance
     real(kind=r_def)                              :: r_tol
     ! absolute tolerance
     real(kind=r_def)                              :: a_tol
     ! maximal number of iterations
     integer(kind=i_def)                           :: max_iter
   contains
     procedure (apply_interface), deferred :: apply
  end type abstract_iterative_solver_type

  abstract interface
!> @brief solve linear system \f$Ax=b\f$ for \f$x\f$
!> @detailled apply the iterative solver for a given right hand side
!> \f$\b\f$
!!
!> @param [inout] x resulting solution \f$x\f$
!> @param [in] b right hand side vector \f$b\f$
     subroutine apply_interface(self, x, b)
       import :: abstract_vector_type
       import :: abstract_iterative_solver_type
       class(abstract_iterative_solver_type), intent(inout) :: self
       class(abstract_vector_type),           intent(inout) :: x
       class(abstract_vector_type),           intent(inout) :: b
     end subroutine apply_interface
  end interface

  !! ---------- End of the abstract type ------------ !!
  !! --  Now what follows are the declarations and -- !!
  !! --  interfaces for of the procedures of the   -- !!
  !! --  extended types. Actual procedures are in  -- !!
  !! --  submodules.                               -- !!
  !! ------------------------------------------------ !!

  !! ---  Conjugate Gradient type declaration and interfaces --- !!
  
  type, public, extends(abstract_iterative_solver_type) :: conjugate_gradient_type
     private
   contains
     procedure :: apply => cg_solve
     procedure :: cg_solve
  end type conjugate_gradient_type

  ! overload the default structure constructor
  interface conjugate_gradient_type
     module procedure cg_constructor
  end interface

  ! the constructor will be in a submodule
  interface 
     module function cg_constructor( lin_op, prec, r_tol, a_tol, max_iter) & 
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       type(conjugate_gradient_type) :: self
     end function
  end interface
  
  interface 
     module subroutine cg_solve(self, x, b)
       class(conjugate_gradient_type), intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine
  end interface

  !! --- BiCGStab type declarations and interfaces --- !!
  type, public, extends(abstract_iterative_solver_type) :: BiCGstab_type
     private
   contains
     procedure :: apply => bicgstab_solve
     procedure :: bicgstab_solve
  end type BiCGstab_type

  ! overload the default structure constructor
  interface BiCGstab_type
     module procedure bicgstab_constructor
  end interface

  interface 
     module function bicgstab_constructor( lin_op, prec, r_tol, a_tol, max_iter) & 
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       type(bicgstab_type) :: self
     end function
  end interface
  
  interface 
     module subroutine bicgstab_solve(self, x, b)
       class(bicgstab_type),           intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine
  end interface  


  !! --- GMRES type declarations and interfaces --- !!
  type, public, extends(abstract_iterative_solver_type) :: gmres_type
     private
     integer(kind=i_def) :: gcrk
   contains
     procedure :: apply => gmres_solve
     procedure :: gmres_solve
  end type gmres_type

  ! overload the default structure constructor
  interface gmres_type
     module procedure gmres_constructor
  end interface

  interface 
     module function gmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter) & 
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       integer(kind=i_def),                          intent(in) :: gcrk
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       type(gmres_type) :: self
     end function gmres_constructor
  end interface
  
  interface 
     module subroutine gmres_solve(self, x, b)
       class(gmres_type),           intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine gmres_solve
  end interface  
  
end module iterative_solver_mod

!  Submodule with procedures for conjugate gradient !!
submodule(iterative_solver_mod) conjugate_gradient_smod
contains
  !> constructs a <code>conjugate_gradient</code> solver
  !! sets the values for the solver such as the residual (r_tol) and
  !! points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @return the constructed conjugate gradient solver
  module function cg_constructor(lin_op, prec, r_tol, a_tol, max_iter) result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    type(conjugate_gradient_type) :: self
    
    self%lin_op => lin_op
    self%prec   => prec
    self%r_tol  = r_tol
    self%a_tol  = a_tol
    self%max_iter    = max_iter
  end function

  !> CG solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !! This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine cg_solve(self, x, b)
    implicit none
    class(conjugate_gradient_type), intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
    class(abstract_vector_type),    intent(inout) :: b

    ! written in terms of abstract types
    integer(kind=i_def) :: iter
    real(kind=r_def)    :: alpha, beta
    real(kind=r_def)    :: r_nrm, r_nrm_0, r_nrm_old, rz, rz_new
    logical             :: converged
    integer :: astat

    ! temporary vectors 
    class(abstract_vector_type), allocatable :: r
    class(abstract_vector_type), allocatable :: p
    class(abstract_vector_type), allocatable :: z

    call x%duplicate(r)
    call r%set_scalar(0.0_r_def)
    call x%duplicate(p)
    call p%set_scalar(0.0_r_def)
    call x%duplicate(z)

    converged=.false.

    !set up the algorithm
    call self%lin_op%apply(x,r) ! r = A.x
    call r%scale(-1.0_r_def)     ! r = -A.x
    call r%axpy(1.0_r_def, b)   ! r = b - A.x
    r_nrm_0 = r%norm()                   ! r_0 = ||r||_2
    write(log_scratch_space,'(A,E15.8)')  &
         "iterative_solver:cg starting ||b|| = ", r_nrm_0
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    call z%set_scalar(0.0_r_def)
    call self%prec%apply(r,z)         ! z = P^{-1}.r    
    rz = r%dot(z)                        ! rz = <r,z>
    r_nrm_old = r_nrm_0
    call p%copy(z)

    write(log_scratch_space,'("iter      ||r_i||        ||r_i||/||r_0||  ||r_i/r_{i-1}||")')
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    ! iterate until maximal number of iterations is reached
    do iter=1, self%max_iter
       call self%lin_op%apply(p,z)       ! z = A.p
       alpha = rz / p%dot(z)             ! alpha = <r,z> / <p,A.p>
       call x%axpy(alpha,p)              ! x -> x + alpha*p
       call r%axpy(-alpha,z)             ! r -> r - alpha*A.p
       r_nrm = r%norm()                  ! r = ||r||_2
       write(log_scratch_space,'(I6, "    ",E12.5,"   ",E12.5,"   ",F8.4)')&
            iter, r_nrm, r_nrm/r_nrm_0, r_nrm/r_nrm_old
    call log_event(log_scratch_space,LOG_LEVEL_INFO)       
       ! exit if either absolute or relative tolerance is reached
       if (      ( r_nrm/r_nrm_0 <= self%r_tol ) &
            .or. ( r_nrm <= self%a_tol ) ) then
          converged=.true.
          exit
       end if
       call self%prec%apply(r,z)         ! z = P^{-1}.r
       rz_new = r%dot(z)                 ! rz_new = <r,z>
       beta = rz_new/rz                  ! beta = <r_{new},z_{new}> / <r,z>
       call p%aypx(beta,z)               ! p -> z + beta*p
       rz = rz_new
       r_nrm_old = r_nrm
    end do
    if (converged) then
       write(log_scratch_space, &
            '("Solver converged after ",I6," iterations")') iter
    else
       write(log_scratch_space, &
            '("Solver failed to converge after ",I6," iterations")') iter
    end if
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

  end subroutine cg_solve

end submodule conjugate_gradient_smod

! submodule for the bicgstab procedures -- !!
submodule(iterative_solver_mod) bicgstab_smod
contains
  !> constructs a <code>bicgstab_type</code> solver
  !! sets the values for the solver such as the residual (r_tol) and
  !! points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @return the constructed conjugate gradient solver  
  module function bicgstab_constructor( lin_op, prec, r_tol, a_tol, max_iter) & 
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    type(bicgstab_type) :: self
    
    write(log_scratch_space,'(A)') "bicgstab_constructor:Hello, World"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    self%lin_op => lin_op
    self%prec   => prec
    self%r_tol  = r_tol
    self%a_tol  = a_tol
    self%max_iter    = max_iter
  end function bicgstab_constructor

  !> bicgstab solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !! This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine bicgstab_solve(self, x, b)
    implicit none
    class(bicgstab_type),           intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
     class(abstract_vector_type),    intent(inout) :: b

    ! tempory vectors
    class(abstract_vector_type), allocatable :: r
    class(abstract_vector_type), allocatable :: p
    class(abstract_vector_type), allocatable :: v
    class(abstract_vector_type), allocatable :: t
    class(abstract_vector_type), allocatable :: s
    class(abstract_vector_type), allocatable :: y
    class(abstract_vector_type), allocatable :: z
    
    !temporary scalars
    real(kind=r_def) :: alpha, beta, rho, omega, norm
    real(kind=r_def) :: err, sc_err, tt, ts

    integer(kind=r_def) :: iter

    ! compute the starting residual
    call x%duplicate(r)
    call r%copy(b)
    sc_err = r%norm()
    sc_err = max(sc_err,0.001_r_def)    
    write( log_scratch_space, '(A,E15.8)' ) &
         "iterative solver: bicgstab starting ... ||b|| = ", sc_err
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    alpha = 1.0_r_def
    omega = 1.0_r_def
    norm  = 1.0_r_def

    call x%duplicate(v)
    call x%duplicate(p)
    call x%duplicate(s)
    call x%duplicate(t)
    call v%set_scalar(0.0_r_def)
    call p%set_scalar(0.0_r_def)

    call x%duplicate(y)
    call x%duplicate(z)

    do iter = 1, self%max_iter
       rho = r%dot(b)
       beta = (rho/norm) * (alpha/omega)
       ! p = r + beta*(p - omega*v)
       call p%aypx(beta, r)
       call p%axpy(-beta*omega, v)
       ! precondition
       call self%prec%apply(p,y)
       ! apply the matrix
       call self%lin_op%apply(y,v)

       norm = b%dot(v)
       alpha = rho/norm

       ! s = r - alpha * v
       call s%copy(r)
       call s%axpy(-alpha,v)
       ! apply the preconditioner
       call self%prec%apply(s,z)
       ! apply the operator
       call self%lin_op%apply(z,t)

       ! final scalars
       tt = t%dot(t)
       ts = t%dot(s)
       omega = ts/tt
       ! final updates
       ! x = x + omega*z + alpha * y
       call x%axpy(omega,z)
       call x%axpy(alpha,y)
       ! compute the residual vector
       call r%copy(s)
       call r%axpy(-omega,t)

       !update the scalars
       norm = rho

       ! check for convergence
       err = r%norm()/sc_err
       write( log_scratch_space, '(A,I2,A, E15.8)' ) "iterative_solver:bicgstab[", &
            iter, "]: res = ", err
       call log_event(log_scratch_space, LOG_LEVEL_INFO)

       if (err < self%r_tol) then
          write( log_scratch_space, '(A, I2, A, E15.8)' ) &
               "iterative_solver:bicgstab:converged in ", iter,    &
             " iters, final=", err
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        exit
      end if
    end do
    
    if(iter >= self%max_iter) then
       write(log_scratch_space, '(A, I3, A, E15.8)') &
           "iterative_solver:bicgstab: NOT converged in", iter, " iters, Res=", err
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
   
  end subroutine bicgstab_solve
  
end submodule bicgstab_smod 
   
submodule(iterative_solver_mod) gmres_smod
contains
  !> constructs a <code>gmres_type</code> solver
  !! sets the values for the solver such as the residual (r_tol) and
  !! points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] gcrk integer, number of internal vectors to use known as the restart value
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_iter, integer the maximum number of iterations
  !> @return the constructed GMRES solver    
  module function gmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter) & 
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    integer(kind=i_def),                          intent(in) :: gcrk    
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    type(gmres_type) :: self
    
    self%lin_op => lin_op
    self%prec   => prec
    self%gcrk   = gcrk
    self%r_tol  = r_tol
    self%a_tol  = a_tol
    self%max_iter    = max_iter
  end function gmres_constructor

  !> gmres_solve. Over-rides the abstract interface to do the actual solve.
  !> @detail The solver implements left-preconditioning, i.e. is solving M{-1}.A.x = M{-1}.b
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !! This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner  
  module subroutine gmres_solve(self, x, b)
    implicit none
    class(gmres_type),              intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
    class(abstract_vector_type),    intent(inout) :: b

    ! temporary vectors
    class(abstract_vector_type), allocatable :: s
    class(abstract_vector_type), allocatable :: w
    class(abstract_vector_type), allocatable :: dx
    class(abstract_vector_type), allocatable :: Ax    
    class(abstract_vector_type), allocatable :: res
    class(abstract_vector_type), allocatable, dimension(:) :: v

    ! temporary scalars
    real(kind=r_def), allocatable, dimension(:)   :: u, g
    real(kind=r_def), allocatable, dimension(:,:) :: h
    real(kind=r_def)                              :: beta, si, ci, nrm, h1, h2, p, q
    real(kind=r_def)                              :: err, sc_err, init_err

    ! iterators
    integer(kind=i_def) :: iv, ivj, iter

    call x%duplicate(dx)
    call dx%set_scalar(0.0_r_def)

    call x%duplicate(Ax)
    call Ax%set_scalar(0.0_r_def)

    call x%duplicate(res) 
    call res%copy(b) ! assumes Ax is zero initial guess
    
    sc_err = res%norm()
    sc_err = max(sc_err,self%a_tol)    
    write( log_scratch_space, '(A,E15.8,":",E15.8)' ) &
         "iterative solver: GMRES starting ... ||b|| = ", b%norm(),sc_err
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    init_err = sc_err

    !initial guess
    call x%duplicate(s)
    call x%duplicate(w)
    call s%set_scalar(0.0_r_def)
    ! pre-cond
    call self%prec%apply(res,s)

    beta = s%norm()
    allocate(v(self%gcrk), source=x)
    call v(1)%copy(s)
    call v(1)%scale(1.0_r_def/beta)
    
    allocate( h(self%gcrk+1, self%gcrk) )
    allocate( g(self%gcrk+1) )
    allocate( u(self%gcrk) )
    
    g(:)   = 0.0_r_def
    g(1)   = beta

    ! initialisation complete, lets go to work.
    do iter = 1, self%max_iter
       do iv = 1, self%gcrk
          
          call w%copy(v(iv))
          ! apply the operator
          call self%lin_op%apply( w, s )
          ! apply the preconditioner
          call self%prec%apply( s, w )

          ! compute the h values
          do ivj = 1, iv
             h(ivj,iv) = v(ivj)%dot(w)
             ! a x y z, z = ax +y
             call w%axpy(-h(ivj,iv), v(ivj))
          end do
          h(iv+1, iv) = w%norm()
          if(iv < self%gcrk) then
             call v(iv+1)%copy(w)
             call v(iv+1)%scale(1.0_r_def/h(iv+1,iv) )
          end if
       end do


       ! Solve (7.23) of Wesseling (see Saad's book)
       do iv = 1, self%gcrk
          nrm    = sqrt( h(iv,iv)*h(iv,iv) + h(iv+1,iv)*h(iv+1,iv) )
          si     = h(iv+1,iv)/nrm
          ci     = h(iv,iv)/nrm
          p      = ci*g(iv) + si*g(iv+1)
          q      = -si*g(iv) + ci*g(iv+1)
          g(iv)   = p
          g(iv+1) = q
          do ivj = iv, self%gcrk
             h1       = ci*h(iv,ivj)   + si*h(iv+1,ivj)
             h2       =-si*h(iv,ivj)   + ci*h(iv+1,ivj)
             h(iv,ivj)   = h1
             h(iv+1,ivj) = h2
          end do
       end do

       u(self%gcrk) = g(self%gcrk)/h(self%gcrk,self%gcrk)
       do iv = self%gcrk-1, 1, -1
          u(iv) = g(iv)
          do ivj = iv+1, self%gcrk
             u(iv) = u(iv) - h(iv,ivj)*u(ivj)
          end do
          u(iv) = u(iv)/h(iv,iv)
       end do

       ! compute the increments
       do iv = 1, self%gcrk
          ! y, x : y = Px
          call s%copy(v(iv))
          call dx%axpy(u(iv), s)
       end do
       
       ! check for convergence
       call self%lin_op%apply(dx, Ax)
       call res%copy(Ax)
       call res%aypx(-1.0_r_def, b)

       beta = res%norm()
       err = beta/sc_err
       if (err < self%r_tol ) then
          write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
               "GMRES solver_algorithm: converged in ", iter,        &
               " iters, init=", init_err,                           &
               " final=", err
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          exit ! break out of loop
       end if

       call self%prec%apply(res,s)
       call v(1)%set_scalar(0.0_r_def)
       call v(1)%axpy(1.0_r_def/beta, s)

       g(:) = 0.0_r_def
       g(1) = beta
    end do

    if( (iter >= self%max_iter .and. err > self%r_tol ) &
         .or. ieee_is_nan(err) ) then
       write( log_scratch_space, '(A, I3, A, E15.8)')    &
            "GMRES solver_algorithm: NOT converged in ", &
            self%max_iter, " iters, Res=", err
       call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    call x%axpy(1.0_r_def, dx)
       
    deallocate(h, g, u)
    deallocate(v)
    
  end subroutine gmres_solve
end submodule gmres_smod
