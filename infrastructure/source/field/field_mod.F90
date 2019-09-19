!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------
!
!> @brief A module providing field related classes.
!>
!> @details Both a representation of a field which provides no access to the
!> underlying data (to be used in the algorithm layer) and an accessor class
!> (to be used in the Psy layer) are provided.


module field_mod

  use constants_mod,      only: r_def, r_double, i_def, i_halo_index, l_def, &
                                str_def
  use halo_routing_collection_mod, &
                          only: halo_routing_collection_type, &
                                halo_routing_collection
  use halo_routing_mod,   only: halo_routing_type
  use fs_continuity_mod,  only: WCHI
  use function_space_mod, only: function_space_type
  use mesh_mod,           only: mesh_type

  use yaxt,               only: xt_redist,  xt_request, &
                                xt_redist_s_exchange, &
                                xt_redist_a_exchange, xt_request_wait
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------


   !> Abstract field type that is the parent of any field type in the field
   !> object hierarchy
   type, extends(linked_list_data_type), public, abstract :: field_parent_type
   contains
   end type field_parent_type

  !> Algorithm layer representation of a field.
  !>
  !> Objects of this type hold all the data of the field privately.
  !> Unpacking the data is done via the proxy type accessed by the Psy layer
  !> alone.
  !>
  type, extends(field_parent_type), public :: field_type
    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer         :: vspace => null( )
    !> Each field also holds an integer enaumerated value for the
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), allocatable :: data( : )
    !> Holds metadata information about a field - like dofmap or routing tables
    type(halo_routing_type), pointer :: halo_routing => null()
    !> Flag that holds whether each depth of halo is clean or dirty (dirty=1)
    integer(kind=i_def), allocatable :: halo_dirty(:)
    !> Flag that determines whether the copy constructor should copy the data.
    !! false to start with, true thereafter.
    logical(kind=l_def) :: data_extant = .false.
    !> Name of the field
    character(str_def) :: name

    ! IO interface procedure pointers

    procedure(write_interface), nopass, pointer            :: write_method => null()
    procedure(read_interface), nopass, pointer             :: read_method => null()
    procedure(checkpoint_write_interface), nopass, pointer :: checkpoint_write_method => null()
    procedure(checkpoint_read_interface), nopass, pointer  :: checkpoint_read_method => null()

  contains

    !> Function to get a proxy with public pointers to the data in a
    !! field_type.
    procedure, public :: get_proxy

    ! Routine to return a deep, but empty copy of a field
    procedure, public :: copy_field_properties

    ! Logging procedures
    procedure, public :: log_field
    procedure, public :: log_dofs
    procedure, public :: log_minmax
    procedure, public :: log_absmax

    !> Function returns the enumerated integer for the functions_space on which
    !! the field lives
    procedure, public :: which_function_space

    !> Function returns a pointer to the function space on which
    !! the field lives
    procedure, public :: get_function_space

    !> Setter for the field write method 
    procedure, public :: set_write_behaviour

    !> Getter for the field write method
    procedure, public :: get_write_behaviour

    !> Setter for the read method 
    procedure, public :: set_read_behaviour

    !> Getter for the read method
    procedure, public :: get_read_behaviour

    !> Setter for the checkpoint method 
    procedure, public :: set_checkpoint_write_behaviour

    !> Setter for the restart method 
    procedure, public :: set_checkpoint_read_behaviour

    !> Routine to return whether field can be checkpointed
    procedure, public :: can_checkpoint

    !> Routine to return whether field can be written
    procedure, public :: can_write

    !> Routine to return whether field can be read
    procedure, public :: can_read

    !> Routine to write field
    procedure         :: write_field

    !> Routine to read field
    procedure         :: read_field

    !> Routine to read a checkpoint netCDF file
    procedure         :: read_checkpoint

    !> Routine to write a checkpoint netCDF file
    procedure         :: write_checkpoint

    !> Returns the name of the field
    procedure         :: get_name

    !> Routine to return the mesh used by this field
    procedure         :: get_mesh
    procedure         :: get_mesh_id
    !> Routine to return the order of the FEM elements
    procedure         :: get_element_order

    !> Overloaded assigment operator
    procedure         :: field_type_assign

    !> Routine to destroy field_type
    procedure         :: field_final

    !> Finalizers for scalar and arrays of field_type objects
    final             :: field_destructor_scalar, &
                         field_destructor_array1d, &
                         field_destructor_array2d

    !> Override default assignment for field_type pairs.
    generic           :: assignment(=) => field_type_assign

  end type field_type

  interface field_type
    module procedure field_constructor
  end interface

   !> A pointer to a field
   !>
   !> We want to be able hold pointers to fields but when these are passed
   !> around, Fortran has a habit of automatically dereferencing them for you.
   !> If we store them in an object that contains the pointer - they won't
   !> get dereferenced.
   !>
   type, extends(field_parent_type), public :: field_pointer_type
     type(field_type), pointer, public :: field_ptr
   contains
     final :: field_pointer_destructor
   end type field_pointer_type

  interface field_pointer_type
    module procedure field_pointer_constructor
  end interface

  !> Psy layer representation of a field.
  !>
  !> This is an accessor class that allows access to the actual field
  !> information with each element accessed via a public pointer.
  !>
  type, public :: field_proxy_type

    private

    !> An unused allocatable integer that prevents an intenal compiler error
    !> with the Gnu Fortran compiler. Adding an allocatable forces the compiler
    !> to accept that the object has a finaliser. It gets confused without it.
    !> This is a workaround for GCC bug id 61767 - when this bug is fixed, the
    !> integer can be removed. 
    integer(kind=i_def), allocatable :: dummy_for_gnu
    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer, public :: vspace => null()
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), public, pointer :: data( : ) => null()
    !> Holds metadata information about a field - like dofmap or routing tables
    type(halo_routing_type), pointer :: halo_routing => null()
    !> pointer to array that holds halo dirtiness
    integer(kind=i_def), pointer :: halo_dirty(:) => null()
    !> Unique identifier used to identify a halo exchange, so the start of an
    !> asynchronous halo exchange can be matched with the end
    type(xt_request) :: halo_request

  contains

    !> Performs a blocking halo exchange operation on the field.
    !> @todo This is temporarily required by PSyclone for initial development
    !! Eventually, the PSy layer will call the asynchronous versions of 
    !! halo_exchange and this function should be removed.
    !! @param[in] depth The depth to which the halos should be exchanged
    procedure, public :: halo_exchange
    !> Starts a halo exchange operation on the field. The halo exchange
    !> is non-blocking, so this call only starts the process. On Return
    !> from this call, outbound data will have been transferred, but no
    !> guarantees are made for in-bound data elements at this stage.
    !! @param[in] depth The depth to which the halos should be exchanged
    procedure, public :: halo_exchange_start

    !> Wait (i.e. block) until the transfer of data in a halo exchange
    !> (started by a call to halo_exchange_start) has completed.
    !! @param[in] depth The depth to which the halos have been exchanged
    procedure, public :: halo_exchange_finish

    !> Perform a global sum operation on the field
    !> @return The global sum of the field values over all ranks
    procedure, public :: get_sum

    !> Calculate the global minimum of the field
    !> @return The minimum of the field values over all ranks
    procedure, public :: get_min

    !> Calculate the global maximum of the field
    !> @return The maximum of the field values over all ranks
    procedure, public :: get_max

    !> Wait (i.e. block) until all current non-blocking reductions
    !> (sum, max, min) are complete.
    !>
    !> We presently have only blocking reductions, so this
    !> subroutine currently returns without waiting.
    procedure reduction_finish

    !> Returns whether the halos at the given depth are dirty or clean
    !! @param[in] depth The depth at which to check the halos
    !! @return True if the halos are dirty or false if they are clean
    procedure is_dirty

    !> Flags all halos as being dirty
    procedure set_dirty

    !> Flags all the halos up the given depth as clean
    !! @param[in] depth The depth up to which to set the halo to clean
    procedure set_clean

  end type field_proxy_type

  !----- Interface for the deferred function in the field_parent_type -----
  interface
   type(field_proxy_type)function get_proxy_interface(self)
     import field_parent_type
     import field_proxy_type
     class(field_parent_type), target, intent(in) :: self
   end function
  end interface

  ! Define the IO interfaces

  abstract interface

    subroutine write_interface(field_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      type(field_proxy_type ), intent(in)  :: field_proxy
    end subroutine write_interface

    subroutine read_interface(field_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)    :: field_name
      type(field_proxy_type ), intent(inout) :: field_proxy
    end subroutine read_interface

    subroutine checkpoint_write_interface(field_name, file_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      character(len=*),        intent(in)  :: file_name
      type(field_proxy_type ), intent(in)  :: field_proxy
    end subroutine checkpoint_write_interface

    subroutine checkpoint_read_interface(field_name, file_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      character(len=*),        intent(in)  :: file_name
      type(field_proxy_type ), intent(inout)  :: field_proxy
    end subroutine checkpoint_read_interface

  end interface

 public :: write_interface
 public :: read_interface
 public :: checkpoint_write_interface
 public :: checkpoint_read_interface

contains

  !> Function to create a proxy with access to the data in the field_type.
  !>
  !> @return The proxy type with public pointers to the elements of
  !> field_type
  type(field_proxy_type ) function get_proxy(self)
    implicit none
    class(field_type), target, intent(in)  :: self

    get_proxy % vspace                 => self % vspace
    get_proxy % data                   => self % data
    get_proxy % halo_routing           => self % halo_routing
    get_proxy % halo_dirty             => self % halo_dirty

  end function get_proxy

  !> Constructor for a field pointer
  !>
  !> @param [in] field_ptr A pointer to the field that is to be
  !>                       stored as a reference
  !> @return The field_pointer type
  function field_pointer_constructor(field_ptr) result(self)
    implicit none

    type(field_type), pointer :: field_ptr
    type(field_pointer_type) :: self

    self%field_ptr => field_ptr
  end function field_pointer_constructor

  !> Finaliser for a field pointer
  !>
  ! The following finaliser doesn't do anything. Without it, the Gnu compiler
  ! tries to create its own, but only ends up producing an Internal Compiler
  ! Error, so its included here to prevent that.
  subroutine field_pointer_destructor(self)
    implicit none
    type(field_pointer_type), intent(inout) :: self
  end subroutine field_pointer_destructor

  !> Construct a <code>field_type</code> object.
  !>
  !> @param [in] vector_space the function space that the field lives on
  !> @param [in] name The name of the field
  !> @return self the field
  !>
  function field_constructor(vector_space, name) result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    type(function_space_type), target, intent(in) :: vector_space
    character(*), optional, intent(in)            :: name

    type(field_type), target :: self
    ! only associate the vspace pointer, copy constructor does the rest.
    self%vspace => vector_space
    self%data_extant = .false.
    if ( self%vspace%which() /= WCHI ) then ! chi fields are never halo exchanged - so don't neeed a routing table
      self%halo_routing => &
        halo_routing_collection%get_halo_routing( self%vspace%get_mesh_id(),&
                                              self%vspace%get_element_order(), &
                                              self%vspace%which() )
    end if
    ! Set the name of the field if given, otherwise default to 'none'
    if (present(name)) then
      self%name = name
    else
      self%name = 'none'
    end if 

  end function field_constructor


  !> Destroy a scalar field_type instance.
  subroutine field_final(self)

    use log_mod, only : log_event, LOG_LEVEL_ERROR

    implicit none

    class(field_type), intent(inout)    :: self

    if(allocated(self%data)) then
      deallocate(self%data)
    end if

    if ( allocated(self%halo_dirty) ) deallocate(self%halo_dirty)

    nullify( self%vspace,                   &
             self%write_method,             &
             self%read_method,              &
             self%checkpoint_write_method,  &
             self%checkpoint_read_method )

  end subroutine field_final



  !> Finalizer for a scalar <code>field_type</code> instance.
  subroutine field_destructor_scalar(self)

    implicit none

    type(field_type), intent(inout)    :: self

    call self%field_final()

  end subroutine field_destructor_scalar

  !> Finalizer for a 1d array of <code>field_type</code> instances.
  subroutine field_destructor_array1d(self)

    implicit none

    type(field_type), intent(inout)    :: self(:)
    integer(i_def) :: i

    do i=lbound(self,1), ubound(self,1)
      call self(i)%field_final()
    end do

  end subroutine field_destructor_array1d

  !> Finalizer for a 2d array of <code>field_type</code> instances.
  subroutine field_destructor_array2d(self)

    implicit none

    type(field_type), intent(inout)    :: self(:,:)
    integer(i_def) :: i,j

    do i=lbound(self,1), ubound(self,1)
      do j=lbound(self,2), ubound(self,2)
        call self(i,j)%field_final()
      end do
    end do

  end subroutine field_destructor_array2d

  !> Create new empty field that inherits the properties of source field
  !>
  !> @param[in]  self  field_type 
  !> @param[out] dest  field_type new field
  subroutine copy_field_properties(self, dest)
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none
    class(field_type), target, intent(in)  :: self
    class(field_type), target, intent(out) :: dest

    type (mesh_type), pointer   :: mesh => null()
    real(kind=r_def), pointer   :: data_ptr( : ) => null()

    dest%vspace => self%vspace
    dest%write_method => self%write_method
    dest%read_method => self%read_method
    dest%checkpoint_write_method => self%checkpoint_write_method
    dest%checkpoint_read_method => self%checkpoint_read_method
    dest%name = self%name
    dest%halo_routing => self%halo_routing

    allocate( dest%data(self%vspace%get_last_dof_halo()) )

    ! Set the data_extant to be .true. now that the data array has been allocated.
    dest%data_extant = .true.

    ! Create a flag for holding whether a halo depth is dirty or not
    mesh=>dest%vspace%get_mesh()
    allocate(dest%halo_dirty(mesh%get_halo_depth()))
    dest%halo_dirty(:) = 1

    nullify(data_ptr)
    nullify(mesh)

  end subroutine copy_field_properties

  !> Assignment operator between field_type pairs.
  !>
  !> @param[out] dest   field_type lhs
  !> @param[in]  source field_type rhs
  subroutine field_type_assign(dest, source)

    implicit none
    class(field_type), intent(in)  :: source
    class(field_type), intent(out) :: dest

    call source%copy_field_properties(dest)

    if(source%data_extant) then
       dest%data(:) = source%data(:)
       dest%halo_dirty(:)=source%halo_dirty(:)
    end if
  end subroutine field_type_assign

  !> Setter for field write behaviour
  !> @param[in,out]  self  field_type 
  !> @param [in] write_behaviour - pointer to procedure implementing write method
  subroutine set_write_behaviour(self, write_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(write_interface), pointer, intent(in) :: write_behaviour
    self%write_method => write_behaviour
  end subroutine set_write_behaviour

  !> Getter to get pointer to field write behaviour
  !> @param[in]  self  field_type 
  !> @param [in] write_behaviour - 
  !>             pointer to procedure implementing write method 
  !> @return pointer to procedure for field write behaviour
  subroutine get_write_behaviour(self, write_behaviour)

    implicit none

    class(field_type), intent(in) :: self
    procedure(write_interface), pointer, intent(inout) :: write_behaviour

    write_behaviour => self%write_method

    return
  end subroutine get_write_behaviour

  !> Setter for read behaviour
  !> @param[in,out]  self  field_type 
  !> @param [in] read_behaviour - pointer to procedure implementing read method
  subroutine set_read_behaviour(self, read_behaviour)
    implicit none
    class(field_type), intent(inout)               :: self
    procedure(read_interface), pointer, intent(in) :: read_behaviour
    self%read_method => read_behaviour
  end subroutine set_read_behaviour

  !> Getter to get pointer to read behaviour
  !> @param[in]  self  field_type 
  !> @param [in] read_behaviour - 
  !>             pointer to procedure implementing read method 
  !> @return pointer to procedure for read behaviour
  subroutine get_read_behaviour(self, read_behaviour)

    implicit none

    class(field_type), intent(in) :: self
    procedure(read_interface), pointer, intent(inout) :: read_behaviour

    read_behaviour => self%read_method

    return
  end subroutine get_read_behaviour


  !> Setter for checkpoint write behaviour
  !>
  !> @param [in] checkpoint_write_behaviour - 
  !>             pointer to procedure implementing checkpoint write method 
  subroutine set_checkpoint_write_behaviour(self, checkpoint_write_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(checkpoint_write_interface), pointer, intent(in)   :: checkpoint_write_behaviour
    self%checkpoint_write_method => checkpoint_write_behaviour
  end subroutine set_checkpoint_write_behaviour

  !> Setter for checkpoint read behaviour
  !>
  !> @param [in] checkpoint_read_behaviour - 
  !>             pointer to procedure implementing checkpoint read method 
  subroutine set_checkpoint_read_behaviour(self, checkpoint_read_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(checkpoint_read_interface), pointer, intent(in)   :: checkpoint_read_behaviour
    self%checkpoint_read_method => checkpoint_read_behaviour
  end subroutine set_checkpoint_read_behaviour

  !> Returns whether field can be checkpointed
  !>
  !> @return .true. or .false.
  function can_checkpoint(self) result(checkpointable)

    implicit none

    class(field_type), intent(in) :: self
    logical(l_def) :: checkpointable

    if (associated(self%checkpoint_write_method) .and. &
       associated(self%checkpoint_read_method)) then
      checkpointable = .true.
    else
      checkpointable = .false.
    end if

  end function can_checkpoint

  !> Returns whether field can be written
  !>
  !> @return .true. or .false.
  function can_write(self) result(writeable)

    implicit none

    class(field_type), intent(in) :: self
    logical(l_def) :: writeable

    writeable = associated(self%write_method)

  end function can_write

  !> Returns whether field can be read
  !>
  !> @return .true. or .false.
  function can_read(self) result(readable)

    implicit none

    class(field_type), intent(in) :: self
    logical(l_def) :: readable

    readable = associated(self%read_method)

  end function can_read

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  !> Gets the name of the field.
  !>
  !> @return field name
  function get_name(self) result(name)

    implicit none

    class(field_type), intent(in) :: self
    character(str_def) :: name

    name = self%name

  end function get_name

  !> Function to get mesh information from the field.
  !>
  !> @return Mesh object
  function get_mesh(self) result(mesh)

    implicit none

    class(field_type), intent(in) :: self
    type(mesh_type), pointer :: mesh

    mesh => self%vspace%get_mesh()

  end function get_mesh

  !> Function to get mesh id from the field.
  !>
  !> @return mesh_id
  function get_mesh_id(self) result(mesh_id)
    implicit none

    class (field_type) :: self
    integer(i_def) :: mesh_id

    mesh_id = self%vspace%get_mesh_id()

    return
  end function get_mesh_id

  !> Function to get element order from the field.
  !>
  !> @return Element order of this field
  function get_element_order(self) result(elem)
    implicit none
    
    class (field_type) :: self
    integer(i_def) :: elem
    
    elem = self%vspace%get_element_order()
    
    return
  end function get_element_order
  !> Sends field contents to the log.
  !!
  !! @param[in] dump_level The level to use when sending the dump to the log.
  !! @param[in] label A title added to the log before the data is written out
  !>
  subroutine log_field( self, dump_level, label )

    use constants_mod, only : r_double, i_def
    use log_mod, only : log_event,         &
                        log_scratch_space, &
                        LOG_LEVEL_INFO,    &
                        LOG_LEVEL_TRACE

    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: dump_level
    character( * ),              intent(in) :: label

    integer(i_def)          :: cell
    integer(i_def)          :: layer
    integer(i_def)          :: df
    integer(i_def), pointer :: map(:) => null()

    write( log_scratch_space, '( A, A)' ) trim( label ), " =["
    call log_event( log_scratch_space, dump_level )

    do cell=1,self%vspace%get_ncell()
      map => self%vspace%get_cell_dofmap( cell )
      do df=1,self%vspace%get_ndf()
        do layer=0,self%vspace%get_nlayers()-1
          write( log_scratch_space, '( I6, I6, I6, E16.8 )' ) &
              cell, df, layer+1, self%data( map( df ) + layer )
          call log_event( log_scratch_space, dump_level )
        end do
      end do
    end do

    call log_event( '];', dump_level )

  end subroutine log_field

  !> Sends the field contents to the log
  !!
  !! @param[in] log_level The level to use for logging.
  !! @param[in] label A title added to the log before the data is written out
  !!
  subroutine log_dofs( self, log_level, label )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: label

    integer(i_def) :: df

    call log_event( label, log_level )

    do df=1,self%vspace%get_undf()
      write( log_scratch_space, '( I6, E16.8 )' ) df,self%data( df )
      call log_event( log_scratch_space, log_level )
    end do

  end subroutine log_dofs

  !> Sends the min/max of a field to the log
  !!
  !! @param[in] log_level The level to use for logging.
  !! @param[in] label A title added to the log before the data is written out
  !!
  subroutine log_minmax( self, log_level, label )

    use log_mod,    only : log_event, log_scratch_space,     &
                           application_log_level => log_level
    use scalar_mod, only : scalar_type
    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: label
    integer(i_def)                          :: undf
    type(scalar_type)                       :: fmin, fmax

    ! If we aren't going to log the min and max then we don't need to
    ! do any further work here.
    if ( log_level < application_log_level() ) return

    undf = self%vspace%get_last_dof_owned()
    fmin = scalar_type( minval( self%data(1:undf) ) )
    fmax = scalar_type( maxval( self%data(1:undf) ) )
 
    write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
         "Min/max ", trim( label ),                   &
         " = ", fmin%get_min(), fmax%get_max()
    call log_event( log_scratch_space, log_level )

  end subroutine log_minmax

  !> Sends the max of the absolute value of field to the log
  !!
  !! @param[in] log_level The level to use for logging.
  !! @param[in] label A title added to the log before the data is written out.
  !!
  subroutine log_absmax( self, log_level, label )

    use log_mod,    only : log_event, log_scratch_space,     &
                           application_log_level => log_level
    use scalar_mod, only : scalar_type
    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: label
    integer(i_def)                          :: undf
    type(scalar_type)                       :: fmax

    ! If we aren't going to log the abs max then we don't need to
    ! do any further work here.
    if ( log_level < application_log_level() ) return

    undf = self%vspace%get_last_dof_owned()
    fmax = scalar_type( maxval( abs(self%data(1:undf)) ) )
 
    write( log_scratch_space, '( A, A, E16.8 )' ) &
         trim( label ), " = ", fmax%get_max()
    call log_event( log_scratch_space, log_level )

  end subroutine log_absmax

  !> Function to integer id of the function space from the field
  !>
  !> @return fs
  function which_function_space(self) result(fs)
    implicit none
    class(field_type), intent(in) :: self
    integer(i_def) :: fs

    fs = self%vspace%which()
    return
  end function which_function_space


  !> Function to get pointer to function space from the field.
  !>
  !> @return vspace
  function get_function_space(self) result(vspace)
    implicit none

    class (field_type), target :: self
    type(function_space_type), pointer :: vspace

    vspace => self%vspace

    return
  end function get_function_space

  !> Calls the underlying IO implementation for writing a field
  !> throws an error if this has not been set
  !> @param [in] field_name - field name / id to write
  subroutine write_field(this, field_name)

    use log_mod,           only : log_event, &
                                  LOG_LEVEL_ERROR

    implicit none 

    class(field_type),   intent(in)     :: this
    character(len=*),    intent(in)     :: field_name

    if (associated(this%write_method)) then

      call this%write_method(trim(field_name), this%get_proxy())

    else

      call log_event( 'Error trying to write field '// trim(field_name) // &
                      ', write_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine write_field

  !> Calls the underlying IO implementation for reading into the field
  !> throws an error if this has not been set
  !> @param [in] field_name - field name / id to read
  subroutine read_field( self, field_name)
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR, &
                                LOG_LEVEL_INFO

    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(in)              :: field_name


    type( field_proxy_type )                      :: tmp_proxy

    if (associated(self%read_method)) then

      tmp_proxy = self%get_proxy()

      call self%read_method(trim(field_name), tmp_proxy)

      ! Set halos dirty here as for parallel read we only read in data for owned
      ! dofs and the halos will not be set

      self%halo_dirty(:) = 1

    else

      call log_event( 'Error trying to read into field '// trim(field_name) // &
                      ', read_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine read_field

  !> Reads a checkpoint file into the field
  !> @param [in] field_name - field name / id to read
  !> @param [in] file_name - file name to read from
  subroutine read_checkpoint( self, field_name, file_name)
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR

    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(in)              :: field_name
    character(len=*),     intent(in)              :: file_name


    type( field_proxy_type )                      :: tmp_proxy


    if (associated(self%checkpoint_read_method)) then

      tmp_proxy = self%get_proxy()

      call self%checkpoint_read_method(trim(field_name), trim(file_name), tmp_proxy)

      ! Set halos dirty here as for parallel read we only read in data for owned
      ! dofs and the halos will not be set

      self%halo_dirty(:) = 1

    else

      call log_event( 'Error trying to read checkpoint for field '// trim(field_name) // &
                      ', checkpoint_read_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine read_checkpoint

  !> Writes a checkpoint file
  !> @param [in] field_name - field name / id to write
  !> @param [in] file_name - file name to write to
  subroutine write_checkpoint( self, field_name, file_name )
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR


    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(in)              :: field_name
    character(len=*),     intent(in)              :: file_name

    if (associated(self%checkpoint_write_method)) then

      call self%checkpoint_write_method(trim(field_name), trim(file_name), self%get_proxy())

    else

      call log_event( 'Error trying to write checkpoint for field '// trim(field_name) // &
                      ', checkpoint_write_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine write_checkpoint


  !! Perform a blocking halo exchange operation on the field
  !!
  subroutine halo_exchange( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    use count_mod,       only : halo_calls
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(xt_redist) :: redist
    type(mesh_type), pointer   :: mesh => null()

    if( self%vspace%is_comms_fs() ) then
      mesh=>self%vspace%get_mesh()
      if( depth > mesh%get_halo_depth() ) &
        call log_event( 'Error in field: '// &
                        'attempt to exchange halos with depth out of range.', &
                        LOG_LEVEL_ERROR )
        ! Start a blocking (synchronous) halo exchange
        redist = self%halo_routing%get_redist(depth)
        call xt_redist_s_exchange(redist, self%data, self%data)

      ! Halo exchange is complete so set the halo dirty flag to say it
      ! is clean (or more accurately - not dirty)
      self%halo_dirty(1:depth) = 0
      ! If a halo counter has been set up, increment it
      if (allocated(halo_calls)) call halo_calls%counter_inc()
    end if

    nullify( mesh )

  end subroutine halo_exchange

  !! Start a halo exchange operation on the field
  !!
  subroutine halo_exchange_start( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(xt_redist) :: redist
    type(mesh_type), pointer   :: mesh => null()

    if( self%vspace%is_comms_fs() ) then

      mesh=>self%vspace%get_mesh()
      if( depth > mesh%get_halo_depth() ) &
        call log_event( 'Error in field: '// &
                        'attempt to exchange halos with depth out of range.', &
                        LOG_LEVEL_ERROR )
      ! Start an asynchronous halo exchange
      redist = self%halo_routing%get_redist(depth)
      call xt_redist_a_exchange(redist, self%data, self%data, self%halo_request)

      nullify( mesh )
    end if

  end subroutine halo_exchange_start

  !! Wait for a halo exchange to complete
  !!
  subroutine halo_exchange_finish( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    use count_mod,       only : halo_calls
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(mesh_type), pointer   :: mesh => null()

    if( self%vspace%is_comms_fs() ) then

      mesh=>self%vspace%get_mesh()
      if( depth > mesh%get_halo_depth() ) &
        call log_event( 'Error in field: '// &
                        'attempt to exchange halos with depth out of range.', &
                        LOG_LEVEL_ERROR )
      ! Wait for the asynchronous halo exchange to complete
      call xt_request_wait(self%halo_request)

      ! Halo exchange is complete so set the halo dirty flag to say it
      ! is clean (or more accurately - not dirty)
      self%halo_dirty(1:depth) = 0
      ! If a halo counter has been set up, increment it
      if (allocated(halo_calls)) call halo_calls%counter_inc()

      nullify( mesh )
    end if

  end subroutine halo_exchange_finish

  !! Start performing a global sum operation on the field
  !!
  function get_sum(self) result (answer)

    use mpi_mod, only: global_sum
    implicit none

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: l_sum
    real(r_def) :: answer

    integer(i_def) :: i

    if( self%vspace%is_comms_fs() ) then

      ! Generate local sum
      l_sum = 0.0
      do i = 1, self%vspace%get_last_dof_owned()
        l_sum = l_sum + self%data(i)
      end do

      call global_sum( l_sum, answer )
    end if
  end function get_sum

  !! Start the calculation of the global minimum of the field
  !!
  function get_min(self) result (answer)

    use mpi_mod, only: global_min
    implicit none

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: l_min
    real(r_def) :: answer

    integer(i_def) :: i

    if( self%vspace%is_comms_fs() ) then

      ! Generate local min
      l_min = self%data(1)
      do i = 2, self%vspace%get_last_dof_owned()
        if( self%data(i) < l_min ) l_min = self%data(i)
      end do

      call global_min( l_min, answer )
    end if

  end function get_min

  !! Start the calculation of the global maximum of the field
  !!
  function get_max(self) result (answer)

    use mpi_mod, only: global_max
    implicit none

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: l_max
    real(r_def) :: answer

    integer(i_def) :: i

    if( self%vspace%is_comms_fs() ) then

      ! Generate local max
      l_max = self%data(1)
      do i = 2, self%vspace%get_last_dof_owned()
        if( self%data(i) > l_max ) l_max = self%data(i)
      end do

      call global_max( l_max, answer )
    end if

  end function get_max

  !! Wait for any current (non-blocking) reductions (sum, max, min) to complete
  !!
  !! We presently have only blocking reductions, so there is
  !! no need to ever call this subroutine. It is left in here to complete the
  !! API so when non-blocking reductions are implemented, we can support them
  subroutine reduction_finish(self)

    implicit none

    class(field_proxy_type), intent(in) :: self

    logical(l_def) :: is_dirty_tmp

    if( self%vspace%is_comms_fs() ) then

      is_dirty_tmp=self%is_dirty(1)    ! reduction_finish currently does nothing.
                                    ! The "self" that is passed in automatically
                                    ! to a type-bound subroutine is not used -
                                    ! so the compilers complain -  have to use
                                    ! it for something harmless.


    end if

  end subroutine reduction_finish

  ! Returns true if a halo depth is dirty
  ! @param[in] depth The depth of halo to inquire about
  function is_dirty(self, depth) result(dirtiness)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class(field_proxy_type), intent(in) :: self
    integer(i_def), intent(in) :: depth
    logical(l_def) :: dirtiness
    type(mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'call to is_dirty() with depth out of range.', &
                      LOG_LEVEL_ERROR )    

    dirtiness = .false.
    if(self%halo_dirty(depth) == 1)dirtiness = .true.

    nullify( mesh )
  end function is_dirty

  ! Sets a halo depth to be flagged as dirty
  ! @param[in] depth The depth up to which to make the halo dirty
  subroutine set_dirty( self )

    implicit none

    class(field_proxy_type), intent(inout) :: self

    self%halo_dirty(:) = 1

  end subroutine set_dirty

  ! Sets the halos up to depth to be flagged as clean
  ! @param[in] depth The depth up to which to make the halo clean
  subroutine set_clean(self, depth)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class(field_proxy_type), intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'call to set_clean() with depth out of range.', &
                      LOG_LEVEL_ERROR )    

    self%halo_dirty(1:depth) = 0
    nullify( mesh )
  end subroutine set_clean

end module field_mod
