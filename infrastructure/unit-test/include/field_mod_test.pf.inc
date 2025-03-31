! Preprocessor is used to template test code.

#define TMPL_PASTE(a) a
#define TMPL_MODULE_NAME TMPL_PASTE(TMPL_PASTE(field_)TMPL_KIND)_mod_test
#define TMPL_TEST_FIELD_PROC TMPL_PASTE(test_field_)TMPL_KIND
#define TMPL_TEST_FIELD_PROXY_PROC TMPL_PASTE(TMPL_PASTE(test_field_)TMPL_KIND)_proxy
#define TMPL_TEST_EMPTY_FIELD_PROC TMPL_PASTE(test_empty_field_)TMPL_KIND
#define TMPL_TEST_CLASS TMPL_PASTE(TMPL_PASTE(field_)TMPL_KIND)_test_type
#define TMPL_FIELD_MOD TMPL_PASTE(TMPL_PASTE(field_)TMPL_KIND)_mod
#define TMPL_FIELD_TYPE TMPL_PASTE(TMPL_PASTE(field_)TMPL_KIND)_type
#define TMPL_FIELD_PROXY_TYPE TMPL_PASTE(TMPL_PASTE(field_)TMPL_KIND)_proxy_type
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

  use, intrinsic :: iso_fortran_env,  only : TMPL_KIND
  use constants_mod,                  only : i_def, str_def, imdi
  use halo_comms_mod,                 only : initialise_halo_comms, &
                                             finalise_halo_comms
  use halo_routing_collection_mod,    only : halo_routing_collection_type, &
                                             halo_routing_collection
  use halo_comms_mod,                 only : halo_routing_type
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use local_mesh_mod,                 only : local_mesh_type
  use mesh_collection_mod,            only : mesh_collection_type, &
                                             mesh_collection
  use mesh_mod,                       only : mesh_type, PLANE_BI_PERIODIC
  use function_space_mod,             only : function_space_type
  use fs_continuity_mod,              only : W2, W3
  use mpi_mod,                        only : global_mpi, &
                                             lfric_comm_type
  use pfunit

  implicit none

  private
  public :: TMPL_TEST_FIELD_PROC, &
            TMPL_TEST_FIELD_PROXY_PROC, &
            TMPL_TEST_EMPTY_FIELD_PROC
  !
  ! pFUnit depends on symbol spillover
  !
  public :: MPITestCase, MPITestParameter


  @TestCase
  type, extends(MPITestCase), public :: TMPL_TEST_CLASS
    private
    type(function_space_type), pointer :: w2_fs => null()
    type(function_space_type), pointer :: w3_fs => null()
    type(local_mesh_type)              :: unit_test_local_mesh
    type(mesh_type), pointer           :: mesh => null()
  contains
    procedure setUp
    procedure tearDown
    procedure get_local_mesh_ptr
  end type TMPL_TEST_CLASS

  integer(i_def), parameter :: element_order_h = 0
  integer(i_def), parameter :: element_order_v = 0

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setUp( this )

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: this

    type(local_mesh_type), pointer :: unit_test_local_mesh_ptr
    type(mesh_type) :: unit_test_mesh
    integer(i_def)  :: mesh_id
    type(lfric_comm_type) :: lfric_comm

    call lfric_comm%set_comm_mpi_val(this%getMpiCommunicator())

    !Store the MPI communicator for later use
    call global_mpi%initialise(lfric_comm)

    ! Initialise halo functionality
    call initialise_halo_comms(lfric_comm)

    ! Create top level mesh collection
    mesh_collection = mesh_collection_type()
    ! Create top level function space collection
    function_space_collection = function_space_collection_type()
    ! Create top level halo_routing collection
    halo_routing_collection = halo_routing_collection_type()

    call this%unit_test_local_mesh%initialise()
    unit_test_local_mesh_ptr => this%get_local_mesh_ptr()

    ! Dummy mesh mod has 9 cells, 3 layers and is uniform in vertical
    unit_test_mesh = mesh_type( PLANE_BI_PERIODIC, unit_test_local_mesh_ptr )
    mesh_id = mesh_collection%add_new_mesh( unit_test_mesh )
    this%mesh => mesh_collection%get_mesh( mesh_id )

    ! Make a function space to be used to create fields in the tests
    this%w2_fs => function_space_collection%get_fs( this%mesh, &
                                                    element_order_h,            &
                                                    element_order_v, W2 )
    this%w3_fs => function_space_collection%get_fs( this%mesh, &
                                                    element_order_h,            &
                                                    element_order_v, W3 )

  end subroutine setUp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tearDown( this )

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: this

    call function_space_collection%clear()
    call mesh_collection%clear()
    call halo_routing_collection%clear()

    ! Finalise halo functionality
    call finalise_halo_comms()
    ! Clear the stored MPI communicator
    call global_mpi%finalise()

  end subroutine tearDown

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_local_mesh_ptr(this) result(unit_test_local_mesh_ptr)

    implicit none

    class(TMPL_TEST_CLASS), intent(inout), target :: this

    type(local_mesh_type), pointer :: unit_test_local_mesh_ptr

    unit_test_local_mesh_ptr => this%unit_test_local_mesh

  end function get_local_mesh_ptr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  @Test( npes=[1] )
  subroutine TMPL_TEST_FIELD_PROC( context )

    use TMPL_FIELD_MOD, only: TMPL_FIELD_TYPE, TMPL_FIELD_PROXY_TYPE

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: context

    type(TMPL_FIELD_TYPE)          :: f1
    type(TMPL_FIELD_TYPE)          :: f2
    type(TMPL_FIELD_PROXY_TYPE)    :: f1_p
    type(TMPL_FIELD_PROXY_TYPE)    :: f2_p

    character(str_def) :: name
    integer(i_def)     :: halo_depth

    type(function_space_type), pointer :: fs => null()
    type(mesh_type),           pointer :: mesh_out => null()
    type(mesh_type),           pointer :: mesh => null()

    integer(i_def) :: order

    ! Make a field with a function space
    call f1%initialise( vector_space = context%w2_fs, name='f1', halo_depth=3 )

    ! Get  and check the the name of the field
    name = f1%get_name()
    @assertEqual( 'f1', trim(name) )

    ! Check getter for halo depth
    halo_depth=f1%get_field_halo_depth()
    @assertEqual ( 3, halo_depth )

    ! get the proxy
    f1_p = f1%get_proxy()

    ! Check the size of the data array is the same as undf from the fs
    @assertEqual( size( f1_p%data ) , f1_p%vspace%get_undf() )
    ! Check the the function space the field lives on is the one it was
    ! created with
    @assertEqual( f1%which_function_space(), W2 )

    ! Check the get function space returns W2 by comparing undf
    fs => f1%get_function_space()
    @assertEqual( fs%get_undf(), context%w2_fs%get_undf() )

    ! Check the order the field thinks it is
    order=f1%get_element_order_h()
    @assertEqual( element_order_h, order )

    order=f1%get_element_order_v()
    @assertEqual( element_order_v, order )

    ! Check what the field thinks the mesh_id is
    mesh => f1%get_mesh()
    @assertTrue( associated(context%mesh, mesh) )

    ! Copy field operator tests follow...
    call f1%copy_field_serial(f2)
    f2_p = f2%get_proxy()

    ! Check that contained vector_space pointers point to same address
    @assertTrue(associated(f1_p%vspace, f2_p%vspace), message="Assigned field_type vspace pointers differ.")

    ! Test copy_field with an optional name
    call f1%copy_field_serial(f2, name='f2')
    name = f2%get_name()
    @assertEqual( 'f2', trim(name) )

  end subroutine TMPL_TEST_FIELD_PROC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  @Test( npes=[1] )
  subroutine TMPL_TEST_FIELD_PROXY_PROC( context )

    use TMPL_FIELD_MOD, only: TMPL_FIELD_TYPE, TMPL_FIELD_PROXY_TYPE

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: context

    type(TMPL_FIELD_TYPE)       :: f1
    type(TMPL_FIELD_TYPE)       :: f2
    type(TMPL_FIELD_PROXY_TYPE) :: f1_p
    type(TMPL_FIELD_PROXY_TYPE) :: f2_p
    integer(i_def)                :: i
    integer(i_def)                :: halo_depth
    type(halo_routing_type), pointer :: halo_routing
    integer(i_def)                :: order
    TMPL_TYPE(TMPL_KIND)                  :: test_sum, test_min, test_max
    logical                       :: dirty

    ! Make a field with a function space
    call f1%initialise( vector_space = context%w2_fs, name='f1', halo_depth=3 )

    ! get the proxy
    f1_p = f1%get_proxy()

    !Initialise field data
    do i=1,f1_p%vspace%get_undf()
      f1_p%data(i)=real(i, kind=TMPL_KIND)
    end do

    ! Make a second field, identical to the first
    call f2%initialise( vector_space = context%w2_fs )
    f2_p = f2%get_proxy()
    do i=1,f2_p%vspace%get_undf()
      f2_p%data(i)=real(i, kind=TMPL_KIND)
    end do

    ! Perform a halo swap on the second field
    call f2_p%halo_exchange_start(1)
    ! If I had any calculations to overlap with comms - I'd put them here!
    call f2_p%halo_exchange_finish(1)

    ! For serial case, the halo swapped field should be the same as the original
    TMPL_ASSERT(f1_p%data, f2_p%data, 1.0e-2)

    ! Test the global reductions
    test_sum=f1_p%get_sum()
    TMPL_ASSERT(TMPL_CAST(4095.0), test_sum, 1.0e-2)

    test_min=f1_p%get_min()
    TMPL_ASSERT(TMPL_CAST(1.0), test_min, 1.0e-2)

    test_max=f1_p%get_max()
    TMPL_ASSERT(TMPL_CAST(90.0), test_max, 1.0e-2)

    halo_depth=f1_p%get_field_proxy_halo_depth()
    @assertEqual ( 3, halo_depth )

    halo_routing=>f1_p%get_halo_routing()
    order = halo_routing%get_element_order_h()
    @assertEqual( element_order_h, order )

    order = halo_routing%get_element_order_v()
    @assertEqual( element_order_v, order )

    dirty=f1_p%is_dirty(0)
    @assertFalse ( dirty )

    call f1_p%set_dirty()
    dirty=f1_p%is_dirty(0)
    @assertFalse ( dirty )

    call f1_p%set_dirty()
    dirty=f1_p%is_dirty(1)
    @assertTrue ( dirty )

    call f1_p%set_clean(1)
    dirty=f1_p%is_dirty(1)
    @assertFalse ( dirty )

  end subroutine TMPL_TEST_FIELD_PROXY_PROC

  @Test( npes=[1] )
  subroutine TMPL_TEST_EMPTY_FIELD_PROC( context )

    use TMPL_FIELD_MOD, only: TMPL_FIELD_TYPE, TMPL_FIELD_PROXY_TYPE

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: context

    type(TMPL_FIELD_TYPE)       :: f1
    type(TMPL_FIELD_PROXY_TYPE) :: f1_p

    TMPL_TYPE(TMPL_KIND), target :: empty_data(1)

    ! Make a field, but provide override data
    call f1%initialise( vector_space = context%w2_fs, name = 'f1', override_data = empty_data )

    ! get the proxy
    f1_p = f1%get_proxy()

    ! Check that the proxy returns the empty_data array
    @assertTrue( associated( f1_p%data, empty_data ) )

    ! Re-initialise as a normal field
    call f1%initialise( vector_space = context%w2_fs, name = 'f1' )

    ! get the proxy
    f1_p = f1%get_proxy()

    ! Check that the proxy does not return the empty_data array
    @assertFalse( associated( f1_p%data, empty_data ) )

    ! Re-initialise again with the override data
    call f1%initialise( vector_space = context%w2_fs, &
                        name = 'f1', override_data = empty_data )

    ! get the proxy
    f1_p = f1%get_proxy()

    ! Check that the proxy returns the empty_data array once again
    @assertTrue( associated( f1_p%data, empty_data ) )

  end subroutine TMPL_TEST_EMPTY_FIELD_PROC

end module TMPL_MODULE_NAME
