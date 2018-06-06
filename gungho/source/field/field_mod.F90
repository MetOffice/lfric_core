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

  use constants_mod,      only: r_def, r_double, i_def, i_halo_index, l_def
  use function_space_mod, only: function_space_type
  use mesh_mod,           only: mesh_type

  use yaxt,               only: xt_redist,  xt_request, &
                                xt_redist_s_exchange, &
                                xt_redist_a_exchange, xt_request_wait

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  !> Algorithm layer representation of a field.
  !>
  !> Objects of this type hold all the data of the field privately.
  !> Unpacking the data is done via the proxy type accessed by the Psy layer
  !> alone.
  !>
  type, public :: field_type
    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer         :: vspace => null( )
    !> Each field also holds an integer enaumerated value for the
    !> function space it will be output on
    integer(kind=i_def)                          :: ospace
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), allocatable :: data( : )
    !> Flag that holds whether each depth of halo is clean or dirty (dirty=1)
    integer(kind=i_def), allocatable :: halo_dirty(:)
    !> Flag that determines whether the copy constructor should copy the data.
    !! false to start with, true thereafter.
    logical(kind=l_def) :: data_extant = .false.

    ! IO interface procedure pointers

    procedure(write_interface), nopass, pointer  :: write_field_method => null()
    procedure(checkpoint_interface), nopass, pointer  :: checkpoint_method => null()
    procedure(restart_interface), nopass, pointer  :: restart_method => null()

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

    !> Function returns the enumerated integer for the functions_space on which
    !! the field lives
    procedure, public :: which_function_space

    !> Function returns the enumerated integer for the output function space
    procedure, public :: which_output_function_space

    !> Function returns a pointer to the function space on which
    !! the field lives
    procedure, public :: get_function_space

    !> Setter for the field write method 
    procedure, public :: set_write_field_behaviour

    !> Getter for the field write method
    procedure         :: get_write_field_behaviour

    !> Setter for the checkpoint method 
    procedure, public :: set_checkpoint_behaviour

    !> Setter for the restart method 
    procedure, public :: set_restart_behaviour

    !> Routine to write field
    procedure         :: write_field

    !> Routine to read a restart netCDF file
    procedure         :: read_restart

    !> Routine to write a checkpoint netCDF file
    procedure         :: write_checkpoint

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
    !> Each field also has a pointer to the function space it will be
    !> output on
    type( function_space_type ), pointer         :: ospace => null( )
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), public, pointer         :: data( : ) => null()
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

  ! Define the IO interfaces

  abstract interface

    subroutine write_interface(field_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      type(field_proxy_type ), intent(in)  :: field_proxy
    end subroutine write_interface

    subroutine checkpoint_interface(field_name, file_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      character(len=*),        intent(in)  :: file_name
      type(field_proxy_type ), intent(in)  :: field_proxy
    end subroutine checkpoint_interface

    subroutine restart_interface(field_name, file_name, field_proxy)
      import r_def, field_proxy_type
      character(len=*),        intent(in)  :: field_name
      character(len=*),        intent(in)  :: file_name
      type(field_proxy_type ), intent(inout)  :: field_proxy
    end subroutine restart_interface

  end interface

 public :: write_interface
 public :: checkpoint_interface
 public :: restart_interface

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
    get_proxy % halo_dirty             => self % halo_dirty

  end function get_proxy


  !> Construct a <code>field_type</code> object.
  !>
  !> @param [in] vector_space the function space that the field lives on
  !> @param [in] output_space the function space used for field output
  !> @return self the field
  !>
  function field_constructor(vector_space, output_space) result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    type(function_space_type), target, intent(in) :: vector_space
    integer(i_def), optional, intent(in)               :: output_space

    type(field_type), target :: self
    ! only associate the vspace pointer, copy constructor does the rest.
    self%vspace => vector_space
    self%data_extant = .false.

    ! Set the output function space if given, otherwise default
    ! to native function space
    if (present(output_space)) then
      self%ospace = output_space
    else
      self%ospace = self%vspace%which()
    end if 

  end function field_constructor


  !> Destroy a scalar field_type instance.
  subroutine field_final(self)

    use log_mod, only : log_event, LOG_LEVEL_ERROR

    implicit none

    class(field_type), intent(inout)    :: self
    integer(i_def) :: rc

    if(allocated(self%data)) then
      deallocate(self%data)
    end if

    if ( allocated(self%halo_dirty) ) deallocate(self%halo_dirty)

    nullify( self%vspace,             &
             self%write_field_method, &
             self%checkpoint_method,  &
             self%restart_method )

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

    integer(i_def) :: rc
    integer(i_def) :: halo_start, halo_finish

    type (mesh_type), pointer   :: mesh => null()
    real(kind=r_def), pointer   :: data_ptr( : ) => null()

    dest%vspace => self%vspace
    dest%ospace = self%ospace
    dest%write_field_method => self%write_field_method
    dest%checkpoint_method => self%checkpoint_method
    dest%restart_method => self%restart_method

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
  !>
  !> @param [in] pointer to procedure implementing write method 
  subroutine set_write_field_behaviour(self, write_field_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(write_interface), pointer, intent(in)   :: write_field_behaviour
    self%write_field_method => write_field_behaviour
  end subroutine set_write_field_behaviour

  !> Getter to get pointer to field write behaviour
  !>
  !> @return pointer to procedure for field write
  subroutine get_write_field_behaviour(self, write_field_behaviour)

    implicit none

    class (field_type) :: self
    procedure(write_interface), pointer, intent(inout)   :: write_field_behaviour

    write_field_behaviour => self%write_field_method

    return
  end subroutine get_write_field_behaviour

  !> Setter for checkpoint behaviour
  !>
  !> @param [in] pointer to procedure implementing checkpoint method 
  subroutine set_checkpoint_behaviour(self, checkpoint_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(checkpoint_interface), pointer, intent(in)   :: checkpoint_behaviour
    self%checkpoint_method => checkpoint_behaviour
  end subroutine set_checkpoint_behaviour

  !> Setter for restart behaviour
  !>
  !> @param [in] pointer to procedure implementing restart method 
  subroutine set_restart_behaviour(self, restart_behaviour)
    implicit none
    class(field_type), intent(inout)                  :: self
    procedure(restart_interface), pointer, intent(in)   :: restart_behaviour
    self%restart_method => restart_behaviour
  end subroutine set_restart_behaviour


  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

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

    use log_mod,    only : log_event, log_scratch_space, LOG_LEVEL_DEBUG
    use scalar_mod, only : scalar_type
    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: label
    integer(i_def)                          :: undf
    type(scalar_type)                       :: fmin, fmax

    undf = self%vspace%get_last_dof_owned()
    fmin = scalar_type( minval( self%data(1:undf) ) )
    fmax = scalar_type( maxval( self%data(1:undf) ) )
 
    write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
         "Min/max ", trim( label ),                   &
         " = ", fmin%get_min(), fmax%get_max()
    call log_event( log_scratch_space, log_level )

  end subroutine log_minmax


  function which_function_space(self) result(fs)
    implicit none
    class(field_type), intent(in) :: self
    integer(i_def) :: fs

    fs = self%vspace%which()
    return
  end function which_function_space

  function which_output_function_space(self) result(fs)
    implicit none
    class(field_type), intent(in) :: self
    integer(i_def) :: fs

    fs = self%ospace
    return
  end function which_output_function_space

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
  subroutine write_field(this, field_name)

    use log_mod,           only : log_event, &
                                  LOG_LEVEL_ERROR

    implicit none 

    class(field_type),   intent(in)     :: this
    character(len=*),    intent(in)     :: field_name

    if (associated(this%write_field_method)) then

      call this%write_field_method(trim(field_name), this%get_proxy())

    else

      call log_event( 'Error trying to write field '// trim(field_name) // &
                      ', write_field_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine write_field

  !> Reads a restart file into the field
  !>
  subroutine read_restart( self, field_name, file_name)
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR

    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(in)              :: field_name
    character(len=*),     intent(in)              :: file_name


    type( field_proxy_type )                      :: tmp_proxy


    if (associated(self%restart_method)) then

      tmp_proxy = self%get_proxy()

      call self%restart_method(trim(field_name), trim(file_name), tmp_proxy)

      ! Set halos dirty here as for parallel read we only read in data for owned
      ! dofs and the halos will not be set

      self%halo_dirty(:) = 1

    else

      call log_event( 'Error trying to restart field '// trim(field_name) // &
                      ', restart_method not set up', LOG_LEVEL_ERROR )
    end if

  end subroutine read_restart

  !> Writes a checkpoint file
  !>
  subroutine write_checkpoint( self, field_name, file_name )
    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR


    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(in)              :: field_name
    character(len=*),     intent(in)              :: file_name

    if (associated(self%checkpoint_method)) then

      call self%checkpoint_method(trim(field_name), trim(file_name), self%get_proxy())

    else

      call log_event( 'Error trying to checkpoint field '// trim(field_name) // &
                      ', checkpoint_method not set up', LOG_LEVEL_ERROR )
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
        redist=self%vspace%get_redist(depth)
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
      redist=self%vspace%get_redist(depth)
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
    type(xt_redist) :: redist
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

    integer(i_def) :: rc
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
