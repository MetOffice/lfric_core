! Preprocessor is used to template test code.

#define TMPL_PASTE(a) a
#define TMPL_MODULE_NAME TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_mod_test
#define TMPL_TEST_CLASS TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_test_type
#define TMPL_SCALAR_MOD TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_mod
#define TMPL_SCALAR_TYPE TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_type
#define TMPL_LITERAL(v) TMPL_PASTE(TMPL_PASTE(v)_)TMPL_KIND
#define integer 1
#define real 2
#if TMPL_TYPE == 1
#  define TMPL_ASSERT(e, f, t) @assertEqual(e, f)
#  define TMPL_CAST nint
#elif TMPL_TYPE == 2
#  define TMPL_ASSERT(e, f, t) @assertEqual(e, f, tolerance=TMPL_LITERAL(t))
#  define TMPL_CAST real
#else
#  define TMPL_ASSERT(e, f, t) @fail("Unrecognised type TMPL_TYPE (Is case correct?)")
#  define TMPL_CAST real
#endif
#undef real
#undef integer

module TMPL_MODULE_NAME

  use, intrinsic :: iso_fortran_env, only : TMPL_KIND
  use constants_mod, only: i_def
  use mpi_mod, only : global_mpi, &
                      lfric_comm_type
  use pfunit

  implicit none

  private
  public :: set_up, tear_down, test_all
  !
  ! pFUnit relies on symbol spillage
  !
  public :: MPITestMethod, MPITestParameter

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  @before
  subroutine set_up(this)

    implicit none

    class(MPITestMethod), intent(inout) :: this
    type(lfric_comm_type) :: lfric_comm

    call lfric_comm%set_comm_mpi_val(this%getMpiCommunicator())

    !Store the MPI communicator for later use
    call global_mpi%initialise(lfric_comm)

  end subroutine set_up

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  @after
  subroutine tear_down(this)

    implicit none

    class(MPITestMethod), intent(inout) :: this

    ! Clear the stored MPI communicator
    call global_mpi%finalise()

  end subroutine tear_down

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  @Test( npes=[2] )
  subroutine test_all( this )

    use TMPL_SCALAR_MOD, only: TMPL_SCALAR_TYPE

    implicit none

    class(MPITestMethod), intent(inout) :: this

    type(TMPL_SCALAR_TYPE) :: s1, s2

    TMPL_TYPE(TMPL_KIND) :: s1_val, s2_val
    TMPL_TYPE(TMPL_KIND) :: sum_val, min_val, max_val

    integer(i_def) :: local_rank

    local_rank  = global_mpi%get_comm_rank()

    s1_val = TMPL_CAST(10.0) + local_rank

    s1%value = s1_val

    ! Check the scalar type has been correctly constructed and that
    ! we have access to the data
    TMPL_ASSERT( s1_val, s1%value, 1.0e-2 )

    ! Test the global reductions
    ! - a bit trivial until we have parallel tests
    sum_val = s1%get_sum()
    TMPL_ASSERT( TMPL_CAST(21.0), sum_val, 1.0e-2 )

    min_val = s1%get_min()
    TMPL_ASSERT( TMPL_CAST(10.0), min_val, 1.0e-2 )

    max_val = s1%get_max()
    TMPL_ASSERT( TMPL_CAST(11.0), max_val, 1.0e-2 )

    ! Test alternative way to create a scalar - by calling the initialiser
    s2_val = TMPL_CAST(20.0) + local_rank
    call s2%initialise(s2_val, global_mpi)
    ! Make sure we can call one of the reductions on this scalar
    sum_val = s2%get_sum()
    TMPL_ASSERT( TMPL_CAST(41.0), sum_val, 1.0e-2 )

  end subroutine test_all

end module TMPL_MODULE_NAME
