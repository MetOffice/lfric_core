!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------

!> @brief Abstract base class for iterative solver.
!>
!> This class can be used as a base class for iterative solvers which solve the
!> system \f$Ax=b\f$. It contains a linear operator and a preconditioner object
!> which are required by each Krylov-subspace solver. The type also defines an
!> interface for the solver application. Common data such as the relative and
!> absolute solver tolerance are stored as data members.

module sci_iterative_solver_mod

  use constants_mod,        only : r_def, i_def, EPS, l_def
  use vector_mod,           only : abstract_vector_type
  use sci_linear_operator_mod, &
                            only : abstract_linear_operator_type
  use sci_preconditioner_mod, &
                            only : abstract_preconditioner_type
  use log_mod,              only : log_event, LOG_LEVEL_INFO, &
                                   LOG_LEVEL_DEBUG,           &
                                   LOG_LEVEL_ERROR,           &
                                   log_scratch_space,         &
                                   log_at_level
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

  implicit none
! Removing the following "private" statement is a workaround for a bug that
! appeared in Intel v19. Every item in the module has an explicit access set,
! so not setting the default has no effect. See ticket #3326 for details
!  private

  !> @brief Abstract solver type for the solver API.
  !>
  type, public, abstract :: abstract_iterative_solver_type
     private
     !> Preconditioner
     class(abstract_preconditioner_type),  pointer :: prec => null()
     !> Linear operator
     class(abstract_linear_operator_type), pointer :: lin_op => null()

     ! relative tolerance
     real(kind=r_def)                              :: r_tol
     ! absolute tolerance
     real(kind=r_def)                              :: a_tol
     ! maximal number of iterations
     integer(kind=i_def)                           :: max_iter
     ! monitor the error
     logical(kind=l_def)                           :: monitor_convergence
     ! fail if solver does not converge
     logical(kind=l_def)                           :: fail_on_non_converged
   contains
     procedure (apply_interface), deferred :: apply
  end type abstract_iterative_solver_type

  abstract interface

    !> @brief Solve linear system \f$Ax=b\f$ for \f$x\f$.
    !>
    !> @details Apply the iterative solver for a given right hand side
    !> \f$b\f$.
    !>
    !> @param[inout] x  Resulting solution \f$x\f$
    !> @param[in] b  Right hand side vector \f$b\f$
    !>
    subroutine apply_interface(self, x, b)
      import :: abstract_vector_type
      import :: abstract_iterative_solver_type
      class(abstract_iterative_solver_type), intent(inout) :: self
      class(abstract_vector_type),           intent(inout) :: x
      class(abstract_vector_type),           intent(inout) :: b
    end subroutine apply_interface
  end interface

  ! ---------- End of the abstract type ------------ !
  ! --  Now what follows are the declarations and -- !
  ! --  interfaces for of the procedures of the   -- !
  ! --  extended types. Actual procedures are in  -- !
  ! --  submodules.                               -- !
  ! ------------------------------------------------ !

  ! ---  Conjugate Gradient type declaration and interfaces --- !

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
     module function cg_constructor( lin_op, prec, r_tol, a_tol, max_iter, &
                                     monitor_convergence, fail_on_non_converged) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
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

  ! --- BiCGStab type declarations and interfaces --- !
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
     module function bicgstab_constructor( lin_op, prec, r_tol, a_tol, max_iter, &
                                           monitor_convergence, fail_on_non_converged) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
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


  ! --- GMRES type declarations and interfaces --- !
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
     module function gmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                        monitor_convergence, fail_on_non_converged) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       integer(kind=i_def),                          intent(in) :: gcrk
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
       type(gmres_type) :: self
     end function gmres_constructor
  end interface

  interface
     module subroutine gmres_solve(self, x, b)
       class(gmres_type),              intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine gmres_solve
  end interface

  ! --- FGMRES type declarations and interfaces --- !
  type, public, extends(abstract_iterative_solver_type) :: fgmres_type
     private
     integer(kind=i_def) :: gcrk
   contains
     procedure :: apply => fgmres_solve
     procedure :: fgmres_solve
  end type fgmres_type

  ! overload the default structure constructor
  interface fgmres_type
     module procedure fgmres_constructor
  end interface

  interface
     module function fgmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                         monitor_convergence, fail_on_non_converged) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       integer(kind=i_def),                          intent(in) :: gcrk
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
       type(fgmres_type) :: self
     end function fgmres_constructor
  end interface

  interface
     module subroutine fgmres_solve(self, x, b)
       class(fgmres_type),              intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine fgmres_solve
  end interface

  ! --- GCR type declarations and interfaces --- !
  type, public, extends(abstract_iterative_solver_type) :: gcr_type
     private
     integer(kind=i_def) :: gcrk
   contains
     procedure :: apply => gcr_solve
     procedure :: gcr_solve
  end type gcr_type

  ! overload the default structure constructor
  interface gcr_type
     module procedure gcr_constructor
  end interface

  interface
     module function gcr_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                      monitor_convergence, fail_on_non_converged) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       integer(kind=i_def),                          intent(in) :: gcrk
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
       type(gcr_type) :: self
     end function gcr_constructor
  end interface

  interface
     module subroutine gcr_solve(self, x, b)
       class(gcr_type),              intent(inout) :: self
       class(abstract_vector_type),  intent(inout) :: x
       class(abstract_vector_type),  intent(inout) :: b
     end subroutine gcr_solve
  end interface

  ! --- BLOCK_GCR type declarations and interfaces --- !
  ! BLOCK_GCR breaks the API and will be removed in ticket #???
  type, public, extends(abstract_iterative_solver_type) :: block_gcr_type
     private
     integer(kind=i_def) :: gcrk
   contains
     procedure :: apply => block_gcr_solve
     procedure :: block_gcr_solve
  end type block_gcr_type

  ! overload the default structure constructor
  interface block_gcr_type
     module procedure block_gcr_constructor
  end interface

  interface
     module function block_gcr_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                            monitor_convergence, fail_on_non_converged ) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       integer(kind=i_def),                          intent(in) :: gcrk
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged

       type(block_gcr_type) :: self
     end function block_gcr_constructor
  end interface

  interface
     module subroutine block_gcr_solve(self, x, b)
       class(block_gcr_type),        intent(inout) :: self
       class(abstract_vector_type),  intent(inout) :: x
       class(abstract_vector_type),  intent(inout) :: b
     end subroutine block_gcr_solve
  end interface


  ! ---  Precondition only type declaration and interfaces --- !

  type, public, extends(abstract_iterative_solver_type) :: precondition_only_type
     private
   contains
     procedure :: apply => precondition_only_solve
     procedure :: precondition_only_solve
  end type precondition_only_type

  ! overload the default structure constructor
  interface precondition_only_type
     module procedure precondition_only_constructor
  end interface

  ! the constructor will be in a submodule
  interface
     module function precondition_only_constructor( lin_op, prec, &
                                                    monitor_convergence ) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       type(precondition_only_type) :: self
     end function
  end interface

  interface
     module subroutine precondition_only_solve(self, x, b)
       class(precondition_only_type), intent(inout) :: self
       class(abstract_vector_type),   intent(inout) :: x
       class(abstract_vector_type),   intent(inout) :: b
     end subroutine
  end interface

  ! ---  Jacobi type declaration and interfaces --- !

  type, public, extends(abstract_iterative_solver_type) :: jacobi_type
     private
     real(kind=r_def)    :: rho_relax ! Overrelaxation factor
   contains
     procedure :: apply => jacobi_solve
     procedure :: jacobi_solve
  end type jacobi_type

  ! overload the default structure constructor
  interface jacobi_type
     module procedure jacobi_constructor
  end interface

  ! the constructor will be in a submodule
  interface
     module function jacobi_constructor( lin_op, prec, r_tol, a_tol, max_iter,       &
                                         monitor_convergence, fail_on_non_converged, &
                                         rho_relax) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
       real(kind=r_def),                             intent(in) :: rho_relax
       type(jacobi_type) :: self
     end function
  end interface

  interface
     module subroutine jacobi_solve(self, x, b)
       class(jacobi_type), intent(inout) :: self
       class(abstract_vector_type),    intent(inout) :: x
       class(abstract_vector_type),    intent(inout) :: b
     end subroutine
  end interface

  ! ---  Chebyshev type declaration and interfaces --- !

  type, public, extends(abstract_iterative_solver_type) :: chebyshev_type
     private
     real(kind=r_def) :: lmin, lmax
   contains
     procedure :: apply => chebyshev_solve
     procedure :: chebyshev_solve
  end type chebyshev_type

  ! overload the default structure constructor
  interface chebyshev_type
     module procedure chebyshev_constructor
  end interface

  ! the constructor will be in a submodule
  interface
     module function chebyshev_constructor( lin_op, prec, r_tol, a_tol, max_iter, &
                                            monitor_convergence, fail_on_non_converged, lmin, lmax) &
          result(self)
       class(abstract_linear_operator_type), target, intent(in) :: lin_op
       class(abstract_preconditioner_type),  target, intent(in) :: prec
       real(kind=r_def),                             intent(in) :: r_tol
       real(kind=r_def),                             intent(in) :: a_tol
       integer(kind=i_def),                          intent(in) :: max_iter
       logical(kind=l_def),                          intent(in) :: monitor_convergence
       logical(kind=l_def),                          intent(in) :: fail_on_non_converged
       real(kind=r_def),                             intent(in) :: lmin
       real(kind=r_def),                             intent(in) :: lmax
       type(chebyshev_type) :: self
     end function
  end interface

  interface
     module subroutine chebyshev_solve(self, x, b)
       class(chebyshev_type),       intent(inout) :: self
       class(abstract_vector_type), intent(inout) :: x
       class(abstract_vector_type), intent(inout) :: b
     end subroutine
  end interface

  ! ---------- End of the extended types ------------!

  ! Derived data type that contains an allocatable
  ! abstract_vector_type. This DDT is used to avoid a CCE bug (between
  ! version 15.0 and 18.01) that cause issue when an array of abstract
  ! type is allocated. The workaround is to have an array of this DDT
  ! and then allocate one at the time the abstract type. See ticket
  ! #4451 for details.
  type, private :: array_abstract_vector_type
    class(abstract_vector_type), allocatable :: vt
  end type array_abstract_vector_type

end module sci_iterative_solver_mod

!  Submodule with procedures for conjugate gradient !
submodule(sci_iterative_solver_mod) conjugate_gradient_smod
contains

  !> Constructs a <code>conjugate_gradient</code> solver.
  !>
  !> Sets the values for the solver such as the residual (<code>r_tol</code>)
  !> and points the linear operator and preconditioner at those passed in.
  !>
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @return the constructed conjugate gradient solver
  !>
  module function cg_constructor(lin_op, prec, r_tol, a_tol, max_iter, &
                                 monitor_convergence, fail_on_non_converged) result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged

    type(conjugate_gradient_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged

  end function

  !> CG solve.
  !>
  !> Over-rides the abstract interface to do the actual solve.
  !>
  !> @param[inout] b  "RHS" or boundary conditions.
  !> @param[inout] x  Solution.
  !>
  module subroutine cg_solve(self, x, b)
    implicit none
    class(conjugate_gradient_type), intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
    class(abstract_vector_type),    intent(inout) :: b

    ! written in terms of abstract types
    integer(kind=i_def) :: iter
    real(kind=r_def)    :: alpha, beta
    real(kind=r_def)    :: r_nrm, r_nrm_0, r_nrm_old, rz, rz_new
    logical(kind=l_def) :: converged

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
    call r%scale(-1.0_r_def)    ! r = -A.x
    call r%axpy(1.0_r_def, b)   ! r = b - A.x
    r_nrm_0 = r%norm()          ! r_0 = ||r||_2

    ! Check if r == 0 (to avoid divide by zero problems)
    if ( r_nrm_0 < EPS ) then
      write( log_scratch_space, '(A,E15.8)' ) &
           "cg converged in 0 iterations... ||b|| = ", r_nrm_0
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      return
    end if

    if ( self%monitor_convergence ) then
      write(log_scratch_space,'(A,E15.8)')  &
           "cg starting ||r|| = ||b - A.x|| = ", r_nrm_0
      call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
    end if

    call z%set_scalar(0.0_r_def)
    call self%prec%apply(r,z)         ! z = P^{-1}.r
    rz = r%dot(z)                        ! rz = <r,z>
    r_nrm_old = r_nrm_0
    call p%copy(z)

    if ( self%monitor_convergence ) then
      write(log_scratch_space,'("iter      ||r_i||        ||r_i||/||r_0||  ||r_i/r_{i-1}||")')
      call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
    end if
    ! iterate until maximal number of iterations is reached
    do iter=1, self%max_iter
       call self%lin_op%apply(p,z)       ! z = A.p
       alpha = rz / p%dot(z)             ! alpha = <r,z> / <p,A.p>
       call x%axpy(alpha,p)              ! x -> x + alpha*p
       call r%axpy(-alpha,z)             ! r -> r - alpha*A.p

       if ( self%monitor_convergence ) then
         r_nrm = r%norm()                  ! r = ||r||_2
         write(log_scratch_space,'(I6, "    ",E12.5,"   ",E12.5,"   ",F8.4)')&
              iter, r_nrm, r_nrm/r_nrm_0, r_nrm/r_nrm_old
         call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
         ! exit if either absolute or relative tolerance is reached
         if (      ( r_nrm/r_nrm_0 <= self%r_tol ) &
              .or. ( r_nrm <= self%a_tol ) ) then
            converged=.true.
            exit
         end if
       end if
       call self%prec%apply(r,z)         ! z = P^{-1}.r
       rz_new = r%dot(z)                 ! rz_new = <r,z>
       beta = rz_new/rz                  ! beta = <r_{new},z_{new}> / <r,z>
       call p%aypx(beta,z)               ! p -> z + beta*p
       rz = rz_new
       r_nrm_old = r_nrm
    end do
    if ( self%monitor_convergence ) then
      if ( converged ) then
         write(log_scratch_space, &
              '("cg converged after ",I6," iterations")') iter
         call log_event(log_scratch_space,LOG_LEVEL_INFO)
      else
         write(log_scratch_space, &
              '("cg failed to converge after ",I6," iterations")') iter
         if ( self%fail_on_non_converged ) then
           call log_event(log_scratch_space,LOG_LEVEL_ERROR)
         else
           call log_event(log_scratch_space,LOG_LEVEL_INFO)
         end if
      end if
    end if

  end subroutine cg_solve

end submodule conjugate_gradient_smod

! submodule for the bicgstab procedures -- !
submodule(sci_iterative_solver_mod) bicgstab_smod
contains

  !> Constructs a <code>bicgstab_type</code> solver.
  !>
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !>
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !>
  module function bicgstab_constructor( lin_op, prec, r_tol, a_tol, max_iter, &
                                        monitor_convergence, fail_on_non_converged) &
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    type(bicgstab_type) :: self

    write(log_scratch_space,'(A)') "bicgstab_constructor:"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    self%lin_op                => lin_op
    self%prec                  => prec
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged

  end function bicgstab_constructor

  !> bicgstab solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine bicgstab_solve(self, x, b)
    implicit none
    class(bicgstab_type),           intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
    class(abstract_vector_type),    intent(inout) :: b

    ! tempory vectors
    class(abstract_vector_type), allocatable :: r
    class(abstract_vector_type), allocatable :: r0
    class(abstract_vector_type), allocatable :: p
    class(abstract_vector_type), allocatable :: v
    class(abstract_vector_type), allocatable :: t
    class(abstract_vector_type), allocatable :: s
    class(abstract_vector_type), allocatable :: y
    class(abstract_vector_type), allocatable :: z

    !temporary scalars
    real(kind=r_def) :: alpha, beta, rho, omega, rho_old
    real(kind=r_def) :: err, sc_err, tt, ts

    integer(kind=i_def) :: iter

    ! compute the starting residual
    ! r = b
    call x%duplicate(r)
    call r%copy(b)

    ! v = Ax
    call x%duplicate(v)
    call v%set_scalar(0.0_r_def)
    call self%lin_op%apply(x,v)
    ! r = b - Ax
    call r%axpy(-1.0_r_def,v)
    ! store initial residual
    call r%duplicate(r0)
    call r0%copy(r)

    sc_err = r%norm()

    ! Check if r == 0 (to avoid divide by zero problems)
    if ( sc_err < EPS ) then
       write( log_scratch_space, '(A,E15.8)' ) &
            "bicgstab converged in 0 iterations... ||b|| = ", sc_err
       call log_event(log_scratch_space, LOG_LEVEL_INFO)
       return
    end if

    sc_err = max(sc_err,self%a_tol)
    if ( self%monitor_convergence ) then
      write( log_scratch_space, '(A,E15.8)' ) &
           " bicgstab starting ... ||r|| = ||b - A.x|| = ", sc_err
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

    alpha = 1.0_r_def
    omega = 1.0_r_def
    rho_old  = 1.0_r_def


    call x%duplicate(p)
    call x%duplicate(s)
    call x%duplicate(t)
    call v%set_scalar(0.0_r_def)
    call p%set_scalar(0.0_r_def)

    call x%duplicate(y)
    call x%duplicate(z)

    do iter = 1, self%max_iter
       rho = r%dot(r0)
       beta = (rho/rho_old) * (alpha/omega)
       ! p =C( r - beta*omega*v) + beta P This is post-conditioning
       ! in two stages
       ! stage 1.
       call t%copy(r)
       call t%axpy(-beta*omega, v)
       ! stage 2 post-condition
       call self%prec%apply(t,y)
       ! now add on beta P
       call p%aypx(beta, y)
       ! apply the matrix
       call self%lin_op%apply(p,v)

       alpha = rho/r0%dot(v)
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
       ! x = x + omega*z + alpha * p
       call x%axpy(omega,z)
       call x%axpy(alpha,p)
       ! compute the residual vector
       call r%copy(s)
       call r%axpy(-omega,t)

       !update the scalars
       rho_old = rho

       ! check for convergence
       if ( self%monitor_convergence ) then
         err = r%norm()/sc_err
         write( log_scratch_space, '(A,I4,A, E15.8)' ) "bicgstab[", &
              iter, "]: res = ", err
         call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

         if (err < self%r_tol) then
            write( log_scratch_space, '(A, I4, A, E15.8)' ) &
               "bicgstab:converged in ", iter,              &
               " iters, final=", err
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          exit
        end if
      end if
    end do

    if(iter >= self%max_iter .and. self%monitor_convergence) then
       write(log_scratch_space, '(A, I3, A, E15.8)') &
           "bicgstab: NOT converged in", iter, " iters, Res=", err
      if ( self%fail_on_non_converged ) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      else
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
      end if
    end if

  end subroutine bicgstab_solve

end submodule bicgstab_smod

submodule(sci_iterative_solver_mod) gmres_smod
contains
  !> constructs a <code>gmres_type</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] gcrk integer, number of internal vectors to use known as the restart value
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_iter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @return the constructed GMRES solver
  module function gmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                     monitor_convergence, fail_on_non_converged) &
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    integer(kind=i_def),                          intent(in) :: gcrk
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    type(gmres_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%gcrk                  = gcrk
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged

  end function gmres_constructor

  !> gmres_solve. Over-rides the abstract interface to do the actual solve.
  !> @detail The solver implements left-preconditioning, i.e. is solving M{-1}.A.x = M{-1}.b
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
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
    class(abstract_vector_type), allocatable :: Ax
    class(abstract_vector_type), allocatable :: res

    ! Use a special DDT to avoid a CCE bug.
    type(array_abstract_vector_type), allocatable, dimension(:) :: v

    ! temporary scalars
    real(kind=r_def), allocatable, dimension(:)   :: u, g
    real(kind=r_def), allocatable, dimension(:,:) :: h
    real(kind=r_def)                              :: beta, si, ci, nrm, h1, h2, p, q, res_norm
    real(kind=r_def)                              :: err, sc_err, init_err

    ! iterators
    integer(kind=i_def) :: iv, ivj, iter

    call x%duplicate(Ax)
    call Ax%set_scalar(0.0_r_def)

    call b%duplicate(res)
    ! compute res = b -Ax ... in stages
    call self%lin_op%apply(x,Ax)
    call res%copy(b)
    call res%axpy(-1.0_r_def,Ax)

    sc_err = res%norm()

    ! Check if r == 0 (to avoid divide by zero problems)
    if ( sc_err < EPS ) then
      write( log_scratch_space, '(A,E15.8)' ) &
           "gmres converged in 0 iterations... ||b|| = ", sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      return
    end if

    if ( self%monitor_convergence ) then
      sc_err = max(sc_err,self%a_tol)
      write( log_scratch_space, '(A,E15.8,":",E15.8)' ) &
           "GMRES starting ... ||r|| = ||b - A.x||", res%norm(),sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      init_err = sc_err
    end if

    !initial guess
    call x%duplicate(s)
    call x%duplicate(w)

    ! Use special DDT to avoid CCE bug.
    allocate(v(self%gcrk))

    do iv = 1, self%gcrk
      allocate(v(iv)%vt, source=x)
    end do

    allocate( h(self%gcrk+1, self%gcrk) )
    allocate( g(self%gcrk+1) )
    allocate( u(self%gcrk) )

    ! initialisation complete, lets go to work.
    do iter = 1, self%max_iter
       call s%set_scalar(0.0_r_def)
       call self%prec%apply(res,s)
       beta = s%norm()
       call v(1)%vt%copy(s)
       call v(1)%vt%scale(1.0_r_def/beta)

       g(:) = 0.0_r_def
       g(1) = beta

       do iv = 1, self%gcrk

          call w%copy(v(iv)%vt)
          ! apply the operator
          call self%lin_op%apply( w, s )
          ! apply the preconditioner
          call self%prec%apply( s, w )

          ! compute the h values
          do ivj = 1, iv
             h(ivj,iv) = v(ivj)%vt%dot(w)
             ! a x y z, z = ax +y
             call w%axpy(-h(ivj,iv), v(ivj)%vt)
          end do
          h(iv+1, iv) = w%norm()
          if(iv < self%gcrk) then
             call v(iv+1)%vt%copy(w)
             call v(iv+1)%vt%scale(1.0_r_def/h(iv+1,iv) )
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

       ! compute the increments and update the solution
       do iv = 1, self%gcrk
          ! y, x : y = Px
          call s%copy(v(iv)%vt)
          call x%axpy(u(iv), s)
       end do

       call Ax%set_scalar(0.0_r_def)
       call self%lin_op%apply(x, Ax)
       call res%copy(Ax)
       call res%aypx(-1.0_r_def, b)

       ! check for convergence
       if ( self%monitor_convergence ) then
         res_norm = res%norm()
         err = res_norm/sc_err
         if (err < self%r_tol ) then
            write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
                 "GMRES solver_algorithm: converged in ", &
                 iter, " iters, init=", init_err, " final=", err
            call log_event( log_scratch_space, LOG_LEVEL_INFO )
            exit ! break out of loop
         end if
       end if

    end do

    if ( self%monitor_convergence ) then
      write( log_scratch_space, '(A, I3, A, E15.8)')    &
             "GMRES solver_algorithm: NOT converged in ", &
              self%max_iter, " iters, Res=", err
      if( (iter >= self%max_iter .and. err > self%r_tol .and. self%fail_on_non_converged ) &
           .or. ieee_is_nan(err) ) then
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      else
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
      end if
    end if

    deallocate(h, g, u)

    do iv = 1, self%gcrk
      deallocate(v(iv)%vt)
    end do

    deallocate(v)


  end subroutine gmres_solve
end submodule gmres_smod

submodule(sci_iterative_solver_mod) fgmres_smod
contains
  !> constructs a <code>fgmres_type</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] gcrk integer, number of internal vectors to use known as the restart value
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_iter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @return the constructed fGMRES solver
  module function fgmres_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                      monitor_convergence, fail_on_non_converged) &
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    integer(kind=i_def),                          intent(in) :: gcrk
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    type(fgmres_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%gcrk                  = gcrk
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged

  end function fgmres_constructor

  !> fgmres_solve. Over-rides the abstract interface to do the actual solve.
  !> @detail The solver implements flexible right-preconditioning, i.e. is solving A.M{-1}.M.x = b
  !>         as 1)  A.M{-1}.y = b
  !>            2)  M{-1}y = x
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine fgmres_solve(self, x, b)
    implicit none
    class(fgmres_type),              intent(inout) :: self
    class(abstract_vector_type),    intent(inout) :: x
    class(abstract_vector_type),    intent(inout) :: b

    ! temporary vectors
    class(abstract_vector_type), allocatable :: s
    class(abstract_vector_type), allocatable :: w
    class(abstract_vector_type), allocatable :: dx
    class(abstract_vector_type), allocatable :: Ax
    class(abstract_vector_type), allocatable :: res

    ! Use a special DDT to avoid a CCE bug.
    type(array_abstract_vector_type), allocatable, dimension(:) :: v, Pv

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

    ! Check if r == 0 (to avoid divide by zero problems)
    if ( sc_err < EPS ) then
      write( log_scratch_space, '(A,E15.8)' ) &
           "fGMRES converged in 0 iterations... ||b|| = ", sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      return
    end if

    if ( self%monitor_convergence ) then
      sc_err = max(sc_err,self%a_tol)
      write( log_scratch_space, '(A,E15.8,":",E15.8)' ) &
           "fGMRES starting ... ||b|| = ", b%norm(),sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      init_err = sc_err
    end if

    ! Initial guess
    call x%duplicate(s)
    call x%duplicate(w)
    call s%copy(res)

    beta = s%norm()

    ! Use special DDT to avoid CCE bug.
    allocate(v(self%gcrk))
    allocate(Pv(self%gcrk))

    do iv = 1, self%gcrk
      allocate(Pv(iv)%vt, source=x)
      allocate(v(iv)%vt, source=x)
    end do

    call v(1)%vt%copy(s)
    call v(1)%vt%scale(1.0_r_def/beta)

    allocate( h(self%gcrk+1, self%gcrk) )
    allocate( g(self%gcrk+1) )
    allocate( u(self%gcrk) )

    g(:)   = 0.0_r_def
    g(1)   = beta

    ! initialisation complete, lets go to work.
    do iter = 1, self%max_iter
       do iv = 1, self%gcrk

          call self%prec%apply(v(iv)%vt, Pv(iv)%vt)
          ! apply the operator
          call self%lin_op%apply( Pv(iv)%vt, s )
          call w%copy(s)
          ! compute the h values
          do ivj = 1, iv
             h(ivj,iv) = v(ivj)%vt%dot(w)
             ! a x y z, z = ax +y
             call w%axpy(-h(ivj,iv), v(ivj)%vt)
          end do
          h(iv+1, iv) = w%norm()
          if(iv < self%gcrk) then
             call v(iv+1)%vt%copy(w)
             call v(iv+1)%vt%scale(1.0_r_def/h(iv+1,iv) )
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
          call dx%axpy(u(iv), Pv(iv)%vt)
       end do

       ! check for convergence
       call self%lin_op%apply(dx, Ax)
       call res%copy(Ax)
       call res%aypx(-1.0_r_def, b)

       beta = res%norm()

       if ( self%monitor_convergence ) then
         err = beta/sc_err
         if (err < self%r_tol ) then
            write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
                 "fGMRES solver_algorithm: converged in ", &
                 iter, " iters, init=", init_err, " final=", err
            call log_event( log_scratch_space, LOG_LEVEL_INFO )
            exit ! break out of loop
         end if
       end if

       call s%copy(res)
       call v(1)%vt%set_scalar(0.0_r_def)
       call v(1)%vt%axpy(1.0_r_def/beta, s)

       g(:) = 0.0_r_def
       g(1) = beta
    end do

    if( (iter >= self%max_iter .and. err > self%r_tol .and. self%monitor_convergence) &
         .or. ieee_is_nan(err) ) then
       write( log_scratch_space, '(A, I3, A, E15.8)')    &
            "fGMRES solver_algorithm: NOT converged in ", &
            self%max_iter, " iters, Res=", err
       if ( self%fail_on_non_converged .or. ieee_is_nan(err) ) then
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       else
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
       end if
    end if

    call x%axpy(1.0_r_def, dx)

    deallocate(h, g, u)

    do iv = 1, self%gcrk
      deallocate(Pv(iv)%vt)
      deallocate(v(iv)%vt)
    end do

    deallocate(v ,Pv)


  end subroutine fgmres_solve
end submodule fgmres_smod

submodule(sci_iterative_solver_mod) gcr_smod
contains
  !> constructs a <code>gcr_type</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] gcrk integer, number of internal vectors to use known as the restart value
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_iter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @return the constructed GCR solver
  module function gcr_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                   monitor_convergence, fail_on_non_converged) &
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    integer(kind=i_def),                          intent(in) :: gcrk
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    type(gcr_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%gcrk                  = gcrk
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged
  end function gcr_constructor

  !> gcr_solve. Over-rides the abstract interface to do the actual solve.
  !> @detail The solver implements right-preconditioning, i.e. is solving A.M{-1}.M.x = b
  !>         as 1)  A.M{-1}.y = b
  !>            2)  M{-1}y = x
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine gcr_solve(self, x, b)
    implicit none
    class(gcr_type),              intent(inout) :: self
    class(abstract_vector_type),  intent(inout) :: x
    class(abstract_vector_type),  intent(inout) :: b

    ! temporary vectors
    class(abstract_vector_type), allocatable :: dx
    class(abstract_vector_type), allocatable :: Ax
    class(abstract_vector_type), allocatable :: res

    ! Use a special DDT to avoid a CCE bug.
    type(array_abstract_vector_type), allocatable, dimension(:) :: v, Pv

    ! temporary scalars
    real(kind=r_def) :: alpha, beta
    real(kind=r_def) :: err, sc_err, init_err

    ! iterators
    integer(kind=i_def) :: iv, ivj, iv_final, iter

    call x%duplicate(dx)
    call x%duplicate(res)
    call x%duplicate(Ax)

    ! initial guess
    call dx%set_scalar(0.0_r_def)
    call self%prec%apply( b, dx )
    call self%lin_op%apply( dx, Ax )
    !res = b - Ax
    call res%copy(b)
    call res%axpy(-1.0_r_def, Ax)

    !sc_err = res%norm()
    sc_err = b%norm()

    ! Check if r == 0 (to avoid divide by zero problems)
    if ( sc_err < EPS ) then
      write( log_scratch_space, '(A,E15.8)' ) &
           "GCR converged in 0 iterations... ||b|| = ", sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      return
    end if

    if ( self%monitor_convergence ) then
      sc_err = max(sc_err,self%a_tol)
      write( log_scratch_space, '(A,E15.8,":",E15.8)' ) &
           "GCR starting ... ||b|| = ", b%norm(),sc_err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
      init_err = sc_err
    end if

    ! Use special DDT to avoid CCE bug.
    allocate(v(self%gcrk))
    allocate(Pv(self%gcrk))

    do iv = 1, self%gcrk
      allocate(Pv(iv)%vt, source=x)
      allocate(v(iv)%vt, source=x)
    end do

    ! initialisation complete, lets go to work.
    do iter = 1, self%max_iter
       do iv = 1, self%gcrk

          ! apply the preconditioner
          call self%prec%apply( res, Pv(iv)%vt )
          ! apply the operator
          call self%lin_op%apply( Pv(iv)%vt, v(iv)%vt )

          do ivj = 1, iv-1
            alpha = v(iv)%vt%dot(v(ivj)%vt)
            call v(iv)%vt%axpy(-alpha, v(ivj)%vt)
            call Pv(iv)%vt%axpy(-alpha, Pv(ivj)%vt)
          end do
          alpha = v(iv)%vt%norm()
          beta  = 1.0_r_def/alpha
          call v(iv)%vt%scale(beta)
          call Pv(iv)%vt%scale(beta)
          alpha = res%dot(v(iv)%vt)
          call dx%axpy(alpha, Pv(iv)%vt)
          call res%axpy(-alpha, v(iv)%vt)

          if ( self%monitor_convergence ) then
            err = res%norm()/sc_err
            iv_final = iv
            if (err < self%r_tol ) exit
          end if
        end do

        if ( self%monitor_convergence ) then
          write( log_scratch_space, '(A, I2, A, I2, A, E12.4, A, E15.8)' ) &
               "GCR solver_algorithm: [",iter,",",iv_final, &
               "], iters, init=", init_err, " final=", err
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          if (err < self%r_tol ) exit
        end if
    end do

    if( (iter >= self%max_iter .and. err > self%r_tol .and. self%monitor_convergence ) &
         .or. ieee_is_nan(err) ) then
       write( log_scratch_space, '(A, I3, A, E15.8)')    &
            "GCR solver_algorithm: NOT converged in ", &
            self%max_iter, " iters, Res=", err
       if ( self%fail_on_non_converged .or. ieee_is_nan(err) ) then
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       else
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
       end if
    end if

    call x%axpy(1.0_r_def, dx)

    do iv = 1, self%gcrk
      deallocate(Pv(iv)%vt)
      deallocate(v(iv)%vt)
    end do

    deallocate(v ,Pv)

  end subroutine gcr_solve
end submodule gcr_smod

submodule(sci_iterative_solver_mod) block_gcr_smod
contains
  !> constructs a <code>block_gcr_type</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] gcrk integer, number of internal vectors to use known as the restart value
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_iter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @return the constructed GCR solver
  module function block_gcr_constructor( lin_op, prec, gcrk, r_tol, a_tol, max_iter, &
                                         monitor_convergence, fail_on_non_converged ) &
       result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    integer(kind=i_def),                          intent(in) :: gcrk
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    type(block_gcr_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%gcrk                  = gcrk
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged

  end function

  !> block_gcr_solve. Over-rides the abstract interface to do the actual solve.
  !> @todo NB This breaks the API through use of field_norms.  Ticket #??? will address this.
  !> @detail The solver implements right-preconditioning, i.e. is solving A.M{-1}.M.x = b
  !>         as 1)  A.M{-1}.y = b
  !>            2)  M{-1}y = x
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine block_gcr_solve(self, x, b)
    implicit none
    class(block_gcr_type),              intent(inout) :: self
    class(abstract_vector_type),  intent(inout) :: x
    class(abstract_vector_type),  intent(inout) :: b

    ! temporary vectors
    class(abstract_vector_type), allocatable :: Ax
    class(abstract_vector_type), allocatable :: res

    ! Use a special DDT to avoid a CCE bug.
    type(array_abstract_vector_type), allocatable, dimension(:) :: v, Pv

    ! temporary scalars
    real(kind=r_def) :: alpha, beta
    real(kind=r_def) :: err, aerr, init_err

    ! iterators
    integer(kind=i_def) :: iv, ivj, iv_final, iter

    ! diagnostics
    real(kind=r_def), allocatable, dimension(:) :: initial_error, final_error, relative_error, scaled_error
    integer(kind=i_def)                         :: n_fields, n

    logical(kind=l_def), allocatable :: converged(:)

    n_fields = x%vector_size()
    allocate( initial_error(n_fields), final_error(n_fields), &
              relative_error(n_fields), scaled_error(n_fields), &
              converged(n_fields))


    call x%duplicate(res)
    call x%duplicate(Ax)

    ! initial guess
    call x%set_scalar(0.0_r_def)
    call self%prec%apply( b, x )
    call self%lin_op%apply( x, Ax )

    call res%copy(b)
    call res%axpy(-1.0_r_def, Ax)

    if ( self%monitor_convergence ) then
      ! if the input field has error smaller than the product of absolute and relative
      ! tolerance, then use this for the initial error.  This will ensure convergence
      ! in a single iteration.
      do n = 1,n_fields
        initial_error(n) = max(b%field_norm(n), self%a_tol*self%r_tol)
      end do
      init_err = sum(initial_error)

      if (log_at_level( log_level_info )) then
        ! Only compute the norm if we are going to write it
        write( log_scratch_space, '(A,E15.8,":",E15.8)' ) &
             "BLOCK_GCR starting ... ||b|| = ", b%norm(),init_err
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
      end if
    end if

    ! Use special DDT to avoid CCE bug.
    allocate(v(self%gcrk))
    allocate(Pv(self%gcrk))

    do iv = 1, self%gcrk
      allocate(Pv(iv)%vt, source=x)
      allocate(v(iv)%vt, source=x)
    end do

    converged(:)=.false.
    ! initialisation complete, lets go to work.
    do iter = 1, self%max_iter
       do iv = 1, self%gcrk

          ! apply the preconditioner
          call self%prec%apply( res, Pv(iv)%vt )
          ! apply the operator
          call self%lin_op%apply( Pv(iv)%vt, v(iv)%vt )

          do ivj = 1, iv-1
            alpha = v(iv)%vt%dot(v(ivj)%vt)
            call v(iv)%vt%axpy(-alpha, v(ivj)%vt)
            call Pv(iv)%vt%axpy(-alpha, Pv(ivj)%vt)
          end do
          alpha = v(iv)%vt%norm()
          beta  = 1.0_r_def/alpha
          call v(iv)%vt%scale(beta)
          call Pv(iv)%vt%scale(beta)
          alpha = res%dot(v(iv)%vt)
          call x%axpy(alpha, Pv(iv)%vt)
          call res%axpy(-alpha, v(iv)%vt)

          if ( self%monitor_convergence ) then

            iv_final = iv
            final_error(:) = 0.0_r_def
            ! In gungho, for the semi-implicit solver the fields are organised:
            ! (pressure, horizontal velocity, vertical velocity).
            ! Since the vertical velocity term generally controls the convergence
            ! this is tested first and so the testing loop goes backwards through the field vector
            do n = n_fields,1,-1
              final_error(n) = res%field_norm(n)
              if (final_error(n)/initial_error(n) < self%r_tol &
                 .or. final_error(n) < self%a_tol )then
                ! This field is converged
                converged(n) = .true.
              else
                ! This field is not converged so don't bother checking any others
                exit
              end if
            end do

            aerr = sum(final_error)
            err = sum(final_error/initial_error)

            if (all(converged))then
              exit
            else
              converged(:)=.false.
            end if

            relative_error = final_error/initial_error
            do n = 1,n_fields
              write( log_scratch_space, '(A, I2, 3E16.8)') 'Intermediate BLOCK_GCR errors (I/F/R): ', &
                 n, initial_error(n), final_error(n), relative_error(n)
              call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
            end do

          end if
        end do
        if ( self%monitor_convergence ) then
          write( log_scratch_space, '(A, I2, A, I2, A, E12.4, A, E15.8, A, E15.8)' ) &
               "BLOCK_GCR solver_algorithm: [",iter,",",iv_final, &
               "], iters, init=", init_err, " final=", err, " abs=", aerr
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          if (all(converged)) exit
        end if
    end do

    if ( self%monitor_convergence ) then
      relative_error = final_error/initial_error
      do n = 1,n_fields
        write( log_scratch_space, '(A, I2, 3E16.8)') 'BLOCK_GCR field errors (Init/Final/Rel): ', &
          n, initial_error(n), final_error(n), relative_error(n)
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
      end do

      if( (iter >= self%max_iter .and. err > self%r_tol) &
           .or. ieee_is_nan(err) ) then
         write( log_scratch_space, '(A, I3, A, E15.8)')    &
              "BLOCK_GCR solver_algorithm: NOT converged in ", &
              self%max_iter, " iters, Res=", err
        if ( self%fail_on_non_converged ) then
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        else
           call log_event( log_scratch_space, LOG_LEVEL_INFO )
        end if
      end if
    end if

    do iv = 1, self%gcrk
      deallocate(Pv(iv)%vt)
      deallocate(v(iv)%vt)
    end do

    deallocate(v ,Pv)

  end subroutine block_gcr_solve
end submodule block_gcr_smod

!  Submodule with procedures for Precondition only !
submodule(sci_iterative_solver_mod) precondition_only_smod
contains
  !> constructs a <code>precondition_only</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @return the constructed conjugate gradient solver
  module function precondition_only_constructor(lin_op, prec, monitor_convergence) result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    type(precondition_only_type) :: self

    self%lin_op              => lin_op
    self%prec                => prec
    self%monitor_convergence = monitor_convergence

    if( .not. self%monitor_convergence ) then
       call log_event("Precondition only: No diagnostic norm output", LOG_LEVEL_INFO)
    end if

  end function

  !> Precondition only solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine precondition_only_solve(self, x, b)
    implicit none
    class(precondition_only_type), intent(inout) :: self
    class(abstract_vector_type),   intent(inout) :: x
    class(abstract_vector_type),   intent(inout) :: b

    class(abstract_vector_type), allocatable     :: res
    real(r_def)                                  :: e, e0
    class(abstract_vector_type), allocatable     :: Ax


    call log_event("Precondition only starting", LOG_LEVEL_DEBUG)
    call self%prec%apply(b,x)         ! x = P^{-1}.b

    if( self%monitor_convergence ) then
       ! Compute initial and final error
       call b%duplicate(res)
       call res%copy(b)
       e0 = res%norm()
       call x%duplicate(Ax)
       call self%lin_op%apply(x, Ax)
       call res%axpy(-1.0_r_def, Ax)
       e = res%norm()
       write(log_scratch_space,'(A,3E15.8)')  &
            "Precondition only error,init, relative: = ", e,e0,e/e0
       call log_event(log_scratch_space,LOG_LEVEL_DEBUG)
    else
       call log_event("Precondition only: ... finished", LOG_LEVEL_DEBUG)
    end if

  end subroutine precondition_only_solve

end submodule precondition_only_smod

!  Submodule with procedures for Jacobi solver !
submodule(sci_iterative_solver_mod) jacobi_smod
contains
  !> Constructs a <code>Jacobi</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @param[in] rho_relax, overrelaxation factor
  !> @return the constructed jacobi solver
  module function jacobi_constructor(lin_op, prec, r_tol, a_tol, max_iter,       &
                                     monitor_convergence, fail_on_non_converged, &
                                     rho_relax) result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    real(kind=r_def),                             intent(in) :: rho_relax
    type(jacobi_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged
    self%rho_relax             = rho_relax

  end function

  !> Jacobi solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine jacobi_solve(self, x, b)
    implicit none
    class(jacobi_type),          intent(inout) :: self
    class(abstract_vector_type), intent(inout) :: x
    class(abstract_vector_type), intent(inout) :: b

    ! written in terms of abstract types
    integer(kind=i_def) :: iter
    real(kind=r_def)    :: r_nrm, r_nrm_0
    logical(kind=l_def) :: converged

    ! temporary vectors
    class(abstract_vector_type), allocatable :: r
    class(abstract_vector_type), allocatable :: z

    real(kind=r_def) :: const

    call b%duplicate(r)
    call b%duplicate(z)
    call r%set_scalar(0.0_r_def)
    call z%set_scalar(0.0_r_def)

    converged=.false.
    const = -self%rho_relax

    ! Solve Ax = b
    do iter=1, self%max_iter
      call self%lin_op%apply(x, r) ! r = Ax
      call r%axpy(-1.0_r_def, b) ! r = r - b = Ax - b

      ! Check for convergence
      if ( self%monitor_convergence ) then
        r_nrm = r%norm()
        if ( iter == 1 ) r_nrm_0 = r_nrm
        write(log_scratch_space, &
             '("  jacobi iteration ",I6,": ||r|| = ",E12.4," ||r||/||r_0|| = ",E12.4)') iter, r_nrm, r_nrm/r_nrm_0
         call log_event(log_scratch_space,LOG_LEVEL_INFO)

        if (      ( r_nrm/r_nrm_0 <= self%r_tol ) &
             .or. ( r_nrm <= self%a_tol ) ) then
           converged=.true.
           exit
        end if
      end if

      call self%prec%apply(r, z) ! z = P^{-1}.r = P^{-1}.(Ax - b)
      call x%axpy(const, z) ! x = x + const*z = x + const*P^{-1}(Ax - b)

    end do

    if ( self%monitor_convergence ) then
      if (converged) then
         write(log_scratch_space, &
              '("jacobi converged after ",I6," iterations")') iter
         call log_event(log_scratch_space,LOG_LEVEL_INFO)
      else
         write(log_scratch_space, &
              '("jacobi failed to converge after ",I6," iterations")') iter
         if ( self%fail_on_non_converged ) then
           ! Reached maximum number of iterations so flag failure to converge as
           ! an error
           call log_event(log_scratch_space,LOG_LEVEL_ERROR)
         else
           ! Only performing a fixed number of iterations to so flag failure to converge as
           ! information only
           call log_event(log_scratch_space,LOG_LEVEL_INFO)
         end if
      end if
    end if

  end subroutine jacobi_solve
end submodule jacobi_smod


submodule(sci_iterative_solver_mod) chebyshev_smod
contains
  !> constructs a <code>chebyshev</code> solver
  !> sets the values for the solver such as the residual (r_tol) and
  !> points the linear operator and preconditioner at those passed in.
  !> @param[in] lin_op The linear operator the solver will use
  !> @param[in] prec The preconditioner the solver will use
  !> @param[in] r_tol real, the relative tolerance halting condition
  !> @param[in] a_tol real, the absolute tolerance halting condition
  !> @param[in] max_inter, integer the maximum number of iterations
  !> @param[in] monitor_convergence Monitor the convergence and error in the
  !>                                solver
  !> @param[in] fail_on_non_converged Exit with error if the solver does not
  !>                                  converge
  !> @param[in] lmin Lower bound on the eigenvalues of the matrix
  !> @param[in] lmax Upper bound on the eigenvalues of the matrix
  !> @return the constructed chebyshev solver
  module function chebyshev_constructor(lin_op, prec, r_tol, a_tol, max_iter,       &
                                        monitor_convergence, fail_on_non_converged, &
                                        lmin, lmax) result(self)
    implicit none
    class(abstract_linear_operator_type), target, intent(in) :: lin_op
    class(abstract_preconditioner_type),  target, intent(in) :: prec
    real(kind=r_def),                             intent(in) :: r_tol
    real(kind=r_def),                             intent(in) :: a_tol
    integer(kind=i_def),                          intent(in) :: max_iter
    logical(kind=l_def),                          intent(in) :: monitor_convergence
    logical(kind=l_def),                          intent(in) :: fail_on_non_converged
    real(kind=r_def),                             intent(in) :: lmin
    real(kind=r_def),                             intent(in) :: lmax
    type(chebyshev_type) :: self

    self%lin_op                => lin_op
    self%prec                  => prec
    self%r_tol                 = r_tol
    self%a_tol                 = a_tol
    self%max_iter              = max_iter
    self%monitor_convergence   = monitor_convergence
    self%fail_on_non_converged = fail_on_non_converged
    self%lmin                  = lmin
    self%lmax                  = lmax

  end function

  !> chebyshev solve. Over-rides the abstract interface to do the actual solve.
  !> @param[inout] b an abstract vector which will be an actual vector of unkown extended type
  !> This the "RHS" or boundary conditions,
  !> @param[inout] x an abstract vector which is the solution
  !> @param[self] The solver which has pointers to the lin_op and preconditioner
  module subroutine chebyshev_solve(self, x, b)
    implicit none
    class(chebyshev_type),       intent(inout) :: self
    class(abstract_vector_type), intent(inout) :: x
    class(abstract_vector_type), intent(inout) :: b

    integer(i_def) :: iter
    real(r_def) :: init_norm, final_norm
    real(r_def) :: a1, a2, w, a_over_b, wa_over_b

    class(abstract_vector_type), allocatable :: z
    class(abstract_vector_type), allocatable :: r
    class(abstract_vector_type), allocatable :: xp
    class(abstract_vector_type), allocatable :: xo

    ! Chebyshev iteration for solving M y = f with preconditioner D

    ! Set initial guess
    call x%set_scalar(0.0_r_def)
    call x%duplicate(xp)
    call xp%set_scalar(0.0_r_def)
    call x%duplicate(xo)
    call xo%set_scalar(0.0_r_def)
    call x%duplicate(z)
    call x%duplicate(r)

    if ( self%monitor_convergence ) init_norm = max(1.0_r_def, b%norm())

    ! Set up scalars
    a1 = 2.0_r_def/(self%lmax - self%lmin)
    a2 = (self%lmax + self%lmin)/(self%lmax - self%lmin)
    a_over_b = a1/a2
    w = 1.0_r_def

    do iter = 1, self%max_iter

      ! r = b-M*xo
      call self%lin_op%apply(xo,z)
      call r%axpby(1.0_r_def, b, -1.0_r_def, z)

      ! z = D^{-1}.r
      call self%prec%apply(r,z)

      ! x = w*(a/b*z+xo) + (1-w)*xp
      w = 1.0_r_def/(1.0_r_def - w/(4.0_r_def*a2**2))
      wa_over_b = w * a_over_b
      call x%axpby(wa_over_b, z, w, xo)
      call x%axpy(1.0_r_def-w,xp)

      if ( iter < self%max_iter ) then
        ! xp = xo
        call xp%copy(xo)

        ! xo = x
        call xo%copy(x)
      end if

    end do
    ! residiual = norm(b - M*x)
    if ( self%monitor_convergence ) then
      call self%lin_op%apply(x,z) ! z = M.x
      call z%axpy(-1.0_r_def, b)  ! z = M.x-b
      final_norm = z%norm()
      write(log_scratch_space, &
          '("chebyshev[",I4,"], redidual = ",E16.8)') self%max_iter, final_norm/init_norm
      if ( self%fail_on_non_converged ) then
        call log_event(log_scratch_space,LOG_LEVEL_ERROR)
      else
        call log_event(log_scratch_space,LOG_LEVEL_INFO)
      end if
    end if

  end subroutine chebyshev_solve

end submodule chebyshev_smod
