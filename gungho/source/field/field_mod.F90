!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
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

  use ESMF

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
    !> Pointer array of type real which holds the values of the field
    real(kind=r_def), pointer             :: data( : ) => null()
    !> The data for each field is held within an ESMF array component
    type(ESMF_Array) :: esmf_array
    !> Flag that holds whether each depth of halo is clean or dirty (dirty=1)
    integer(kind=i_def), allocatable :: halo_dirty(:)
    !> Flag that determines whether the copy constructor should copy the data.
    !! false to start with, true thereafter.
    logical(kind=l_def) :: extant = .false.

    ! IO interface procedure pointers

    procedure(write_interface), nopass, pointer  :: write_field_method => null()

  contains

    !> Function to get a proxy with public pointers to the data in a
    !! field_type.
    procedure, public :: get_proxy

    procedure, public :: log_field
    procedure, public :: log_dofs
    procedure, public :: log_minmax

    !> function returns the enumerated integer for the functions_space on which
    !! the field lives
    procedure, public :: which_function_space

    !> function returns the enumerated integer for the output function space
    procedure, public :: which_output_function_space

    !> function returns a pointer to the function space on which
    !! the field lives
    procedure, public :: get_function_space

    !> Setter for the field write method 
    procedure, public :: set_write_field_behaviour
    !> Getter for the field write method
    procedure         :: get_write_field_behaviour

    !> Routine to write field
    procedure         :: write_field

    !> Routine to read a restart netCDF file
    procedure         :: read_restart

    !> Routine to write a checkpoint netCDF file
    procedure         :: write_checkpoint

    !> Routine to return the mesh used by this field
    procedure         :: get_mesh
    procedure         :: get_mesh_id

    !> Overloaded assigment operator
    procedure         :: field_type_assign

    !> Routine to destroy field_type
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
  !> This is an accessor class that allows access to the actual field information
  !> with each element accessed via a public pointer.
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
    !> pointer to the ESMF array
    type(ESMF_Array), pointer :: esmf_array => null()
    !> pointer to array that holds halo dirtiness
    integer(kind=i_def), pointer :: halo_dirty(:) => null()

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
    !> ESMF have only implemented blocking reductions, so this
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

  end interface

 public :: write_interface

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
    get_proxy % esmf_array             => self % esmf_array

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
    self%extant = .false.

    ! Set the output function space if given, otherwise default
    ! to native function space
    if (present(output_space)) then
      self%ospace = output_space
    else
      self%ospace = self%vspace%which()
    end if 
    
  end function field_constructor

  !> Destroy a scalar <code>field_type</code> instance.
  subroutine field_destructor_scalar(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none
    type(field_type), intent(inout)    :: self
    integer(i_def) :: rc

    nullify(self%vspace)
    if(associated(self%data)) then
      call ESMF_ArrayDestroy(self%esmf_array, noGarbage=.TRUE., rc=rc)
      if (rc /= ESMF_SUCCESS ) &
        call log_event( "ESMF failed to destroy a field.", &
                        LOG_LEVEL_ERROR )
    end if

    if(allocated(self%halo_dirty)) then
      deallocate(self%halo_dirty)
    end if

  end subroutine field_destructor_scalar

  !> Destroy a 1d array of <code>field_type</code> instances.
  subroutine field_destructor_array1d(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none
    type(field_type), intent(inout)    :: self(:)
    integer(i_def) :: i
    integer(i_def) :: rc

    do i=lbound(self,1), ubound(self,1)
      nullify(self(i)%vspace)
      if(associated(self(i)%data)) then
        call ESMF_ArrayDestroy(self(i)%esmf_array, noGarbage=.TRUE., rc=rc)
        if (rc /= ESMF_SUCCESS ) &
          call log_event( "ESMF failed to destroy a 1d array of fields.", &
                          LOG_LEVEL_ERROR )
      end if
      if(allocated(self(i)%halo_dirty)) then
        deallocate(self(i)%halo_dirty)
      end if
    end do

  end subroutine field_destructor_array1d

  !> Destroy a 2d array of <code>field_type</code> instances.
  subroutine field_destructor_array2d(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none
    type(field_type), intent(inout)    :: self(:,:)
    integer(i_def) :: i,j
    integer(i_def) :: rc

    do i=lbound(self,1), ubound(self,1)
      do j=lbound(self,2), ubound(self,2)
        nullify(self(i,j)%vspace)
        if(associated(self(i,j)%data)) then
          call ESMF_ArrayDestroy(self(i,j)%esmf_array, noGarbage=.TRUE., rc=rc)
          if (rc /= ESMF_SUCCESS ) &
            call log_event( "ESMF failed to destroy a 2d array of fields.", &
                            LOG_LEVEL_ERROR )
        end if
        if(allocated(self(i,j)%halo_dirty)) then
          deallocate(self(i,j)%halo_dirty)
        end if
      end do
    end do

  end subroutine field_destructor_array2d

  !> Assignment operator between field_type pairs.
  !>
  !> @param[out] dest   field_type lhs
  !> @param[in]  source field_type rhs
  subroutine field_type_assign(dest, source)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none
    class(field_type), intent(out)     :: dest
    class(field_type), intent(in)      :: source

    integer(i_halo_index), allocatable :: global_dof_id(:)
    integer(i_def) :: rc
    integer(i_def) :: halo_start, halo_finish

    type (mesh_type), pointer   :: mesh => null()

    dest%vspace => source%vspace
    dest%ospace = source%ospace
    dest%write_field_method => source%write_field_method

    allocate(global_dof_id(source%vspace%get_last_dof_halo()))
    call source%vspace%get_global_dof_id(global_dof_id)

    halo_start  = source%vspace%get_last_dof_owned()+1
    halo_finish = source%vspace%get_last_dof_halo()
    !If this is a serial run (no halos), halo_start is out of bounds - so fix it
    if(halo_start > source%vspace%get_last_dof_halo())then
      halo_start  = source%vspace%get_last_dof_halo()
      halo_finish = source%vspace%get_last_dof_halo() - 1
    end if
    ! Create an ESMF array - this allows us to perform halo exchanges
    ! This call allocates the memory for the field data - we can extract a
    ! pointer to that allocated memory next
    dest%esmf_array = &
      ESMF_ArrayCreate( distgrid=source%vspace%get_distgrid(), &
                        typekind=ESMF_TYPEKIND_R8, &
                        haloSeqIndexList= &
                                      global_dof_id( halo_start:halo_finish ), &
                        rc=rc )

    ! Extract and store the pointer to the fortran array
    if (rc == ESMF_SUCCESS) &
      call ESMF_ArrayGet(array=dest%esmf_array, farrayPtr=dest%data, rc=rc)

    if (rc /= ESMF_SUCCESS) &
      call log_event( 'ESMF failed to allocate space for field data.', &
                      LOG_LEVEL_ERROR )

    deallocate(global_dof_id)

    ! Create a flag for holding whether a halo depth is dirty or not
    ! and initialise it with a copy of the source data
    mesh=>dest%vspace%get_mesh()
    allocate(dest%halo_dirty(mesh%get_halo_depth()))
    dest%halo_dirty(:) = 1

    if(source%extant) then
       dest%data(:) = source%data(:)
       dest%halo_dirty(:)=source%halo_dirty(:)
    end if
    dest%extant = .true. 

  end subroutine field_type_assign

  !> Setter for field write behaviour
  !>
  !> @param [in] pointer to procedure implementing write method 
  subroutine set_write_field_behaviour(self, write_field_behaviour)
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

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  !> Function to get mesh information from the field.
  !>
  !> @return Mesh object
  function get_mesh(self) result(mesh)

    implicit none

    class (field_type) :: self
    type (mesh_type), pointer :: mesh

    mesh => self%vspace%get_mesh()

    return
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
    integer(i_def), pointer :: map( : )

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
    type( function_space_type ), pointer :: vspace

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
  subroutine read_restart( self, file_name)
    use field_io_ncdf_mod,    only : field_io_ncdf_type

    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(IN)              :: file_name

    type(field_io_ncdf_type), allocatable :: ncdf_file

    allocate(ncdf_file)

    call ncdf_file%file_open( file_name )

    call ncdf_file%read_field_data ( self % data(:) )

    call ncdf_file%file_close()

    deallocate(ncdf_file)

  end subroutine read_restart

  !> Writes a checkpoint netCDF file
  !>
  subroutine write_checkpoint( self, file_name )
    use field_io_ncdf_mod,    only : field_io_ncdf_type

    implicit none

    class( field_type ),  target, intent( inout ) :: self
    character(len=*),     intent(IN)              :: file_name

    type(field_io_ncdf_type), allocatable :: ncdf_file

    allocate(ncdf_file)

    call ncdf_file%file_new( file_name )

    call ncdf_file%write_field_data ( self % data(:) )

    call ncdf_file%file_close()

    deallocate(ncdf_file)

  end subroutine write_checkpoint



  !! Perform a blocking halo exchange operation on the field
  !!
  subroutine halo_exchange( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(ESMF_RouteHandle) :: haloHandle
    integer(i_def) :: rc
    type (mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'attempt to exchange halos with depth out of range.', &
                      LOG_LEVEL_ERROR )

    haloHandle=self%vspace%get_haloHandle(depth)
    call ESMF_ArrayHalo( self%esmf_array, &
                         routehandle=haloHandle, &
                         routesyncflag=ESMF_ROUTESYNC_BLOCKING, &
                         rc=rc )

    if (rc /= ESMF_SUCCESS) call log_event( &
       'ESMF failed to perform the halo exchange.', &
       LOG_LEVEL_ERROR )

    ! Halo exchange is complete so set the halo dirty flag to say it
    ! is clean (or more accurately - not dirty)
    self%halo_dirty(1:depth) = 0

  end subroutine halo_exchange

  !! Start a halo exchange operation on the field
  !!
  subroutine halo_exchange_start( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(ESMF_RouteHandle) :: haloHandle
    integer(i_def) :: rc
    type (mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'attempt to exchange halos with depth out of range.', &
                      LOG_LEVEL_ERROR )

    haloHandle=self%vspace%get_haloHandle(depth)
    call ESMF_ArrayHalo( self%esmf_array, &
                         routehandle=haloHandle, &
                         routesyncflag=ESMF_ROUTESYNC_NBSTART, &
                         rc=rc )

    if (rc /= ESMF_SUCCESS) call log_event( &
       'ESMF failed to start the halo exchange.', &
       LOG_LEVEL_ERROR )

  end subroutine halo_exchange_start

  !! Wait for a halo exchange to complete
  !!
  subroutine halo_exchange_finish( self, depth )

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    class( field_proxy_type ), target, intent(inout) :: self
    integer(i_def), intent(in) :: depth
    type(ESMF_RouteHandle) :: haloHandle
    integer(i_def) :: rc
    type (mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'attempt to exchange halos with depth out of range.', &
                      LOG_LEVEL_ERROR )

    haloHandle=self%vspace%get_haloHandle(depth)
    call ESMF_ArrayHalo( self%esmf_array, &
                         routehandle=haloHandle, &
                         routesyncflag=ESMF_ROUTESYNC_NBWAITFINISH, &
                         rc=rc )

    if (rc /= ESMF_SUCCESS) call log_event( &
       'ESMF failed to finish the halo exchange.', &
       LOG_LEVEL_ERROR )

    ! Halo exchange is complete so set the halo dirty flag to say it
    ! is clean (or more accurately - not dirty)
    self%halo_dirty(1:depth) = 0

  end subroutine halo_exchange_finish

  !! Start performing a global sum operation on the field
  !!
  function get_sum(self) result (answer)

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: answer

    type(ESMF_VM) :: vm
    integer(i_def) :: rc

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
! Currently ESMF has only implemented blocking reductions. Using anything
! other than ESMF_SYNC_BLOCKING for syncflag results in an error
    call ESMF_VMAllFullReduce(vm, &
                              self%data, &
                              answer, &
                              self%vspace%get_last_dof_owned(), &
                              ESMF_REDUCE_SUM, &
                              syncflag = ESMF_SYNC_BLOCKING, &
                              rc=rc)
  end function get_sum

  !! Start the calculation of the global minimum of the field
  !!
  function get_min(self) result (answer)

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: answer

    type(ESMF_VM) :: vm
    integer(i_def) :: rc

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
! Currently ESMF has only implemented blocking reductions. Using anything
! other than ESMF_SYNC_BLOCKING for syncflag results in an error
    call ESMF_VMAllFullReduce(vm, &
                              self%data, &
                              answer, &
                              self%vspace%get_last_dof_owned(), &
                              ESMF_REDUCE_MIN, &
                              syncflag = ESMF_SYNC_BLOCKING, &
                              rc=rc)
  end function get_min

  !! Start the calculation of the global maximum of the field
  !!
  function get_max(self) result (answer)

    class(field_proxy_type), intent(in) :: self

    real(r_def) :: answer

    type(ESMF_VM) :: vm
    integer(i_def) :: rc

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
! Currently ESMF has only implemented blocking reductions. Using anything
! other than ESMF_SYNC_BLOCKING for syncflag results in an error
    call ESMF_VMAllFullReduce(vm, &
                              self%data, &
                              answer, &
                              self%vspace%get_last_dof_owned(), &
                              ESMF_REDUCE_MAX, &
                              syncflag = ESMF_SYNC_BLOCKING, &
                              rc=rc)
  end function get_max

  !! Wait for any current (non-blocking) reductions (sum, max, min) to complete
  !!
  !! Currently, ESMF has only implemented blocking reductions, so there is
  !! no need to ever call this subroutine. It is left in here to complete the
  !! API so when non-blocking reductions are implemented, we can support them
  subroutine reduction_finish(self)

    class(field_proxy_type), intent(in) :: self

    type(ESMF_VM) :: vm
    integer(i_def) :: rc
    logical(l_def) :: is_dirty_tmp

    is_dirty_tmp=self%is_dirty(1)    ! reduction_finish currently does nothing.
                                    ! The "self" that is passed in automatically
                                    ! to a type-bound subroutine is not used -
                                    ! so the compilers complain -  have to use
                                    ! it for something harmless.

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    call ESMF_VMCommWaitAll(vm=vm, rc=rc)

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
    type (mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'call to is_dirty() with depth out of range.', &
                      LOG_LEVEL_ERROR )    

    dirtiness = .false.
    if(self%halo_dirty(depth) == 1)dirtiness = .true.

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
    type (mesh_type), pointer   :: mesh => null()

    mesh=>self%vspace%get_mesh()
    if( depth > mesh%get_halo_depth() ) &
      call log_event( 'Error in field: '// &
                      'call to set_clean() with depth out of range.', &
                      LOG_LEVEL_ERROR )    

    self%halo_dirty(1:depth) = 0

  end subroutine set_clean

end module field_mod
