!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides access to Oasis coupling functionality
!>
module coupling_mod

#ifdef MCT
  use mod_oasis,                     only : oasis_init_comp,        &
                                            oasis_get_localcomm, oasis_abort,  &
                                            oasis_terminate, oasis_enddef,     &
                                            oasis_def_var, oasis_def_partition,&
                                            oasis_out, prism_ok, nnamcpl,      &
                                            namsrcfld, namdstfld, oasis_in,    &
                                            oasis_get_ncpl, oasis_get_freqs,   &
                                            prism_real
#endif

  use constants_mod,                 only : i_def, r_def, str_def, i_halo_index
  use field_collection_iterator_mod, only : field_collection_iterator_type
  use field_collection_mod,          only : field_collection_type
  use field_parent_mod,              only : field_parent_type
  use fs_continuity_mod,             only : W3
  use function_space_mod,            only : function_space_type
  use function_space_collection_mod, only : function_space_collection
  use key_value_collection_mod,      only : key_value_collection_type
  use key_value_mod,                 only : abstract_value_type
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_ERROR,   &
                                            LOG_LEVEL_INFO,   &
                                            LOG_LEVEL_DEBUG
  use mesh_mod,                      only : mesh_type
  use mpi_mod,                       only : mpi_type, &
                                            lfric_comm_type
  use sort_mod,                      only : bubble_sort

  implicit none

  private
  public get_coupling_fields, get_coupling_from_collection

  ! Maximum number of model components that can be coupled together
  integer(i_def), parameter    :: nmax = 8
  ! Format for writing the category number
  character(len=6), parameter :: cpl_fmt = "(i2.2)"
  ! Suffix for multi-category fields (multidata fields)
  character(len=4), parameter  :: cpl_cat = "_cat"

  public cpl_cat

  type, extends(abstract_value_type), public :: coupling_type

    private
    ! OASIS component id
    integer(i_def)                        :: comp_id
    ! OASIS partition id for 0d coupling
    integer(i_def)                        :: part_0d_id
    ! OASIS partition id for 2d coupling
    integer(i_def)                        :: part_2d_id
    ! Total length of the data used in 2d coupling
    integer(i_def)                        :: cpl_size
    ! Index to sort data for sending
    integer(i_def), allocatable           :: local_index(:)

  contains

    procedure, public :: initialise
    procedure, public :: define_partitions
    procedure, public :: define_variables
    procedure, public :: get_local_index
    procedure, public :: finalise

  end type coupling_type

contains


  !> @brief Initialises OASIS coupler
  !>
  !> @param [out] comm_out Communicator returned from OASIS to run the model in
  !> @param [in]  comm_in  Input communicator that OASIS can split
  !> @param [inout] comm_is_split Returns true if the MPI Comm has been split
  !
  subroutine initialise(self, cpl_name, comm_out, comm_in, comm_is_split)
    implicit none
    class(coupling_type),  intent(inout) :: self
    character(*),          intent(in)    :: cpl_name
    type(lfric_comm_type), intent(out)   :: comm_out
    type(lfric_comm_type), intent(in)    :: comm_in
    logical,               intent(inout) :: comm_is_split
#ifdef MCT
    integer(i_def)                :: comm  ! Communicator returned from Oasis
    integer(i_def)                :: kinfo ! error return by OASIS

    call oasis_init_comp (self%comp_id,  &
                         trim(cpl_name), &
                         kinfo,          &
                         commworld=comm_in%get_comm_mpi_val())

    if (kinfo .NE. prism_ok) then
      call oasis_abort(self%comp_id, trim(cpl_name), 'initialise')
    endif

    call oasis_get_localcomm(comm, kinfo)
    call comm_out%set_comm_mpi_val(comm)

    if (kinfo .NE. prism_ok) then
      call oasis_abort(self%comp_id, trim(cpl_name), 'initialise')
    endif

    comm_is_split = .true.

#else
    call comm_out%set_comm_mpi_val( -1 )
    comm_is_split = .false.

#endif
  end subroutine initialise


  !>@brief Defines the types of grid that can be used for coupling
  !> @param [in,out] mpi        The MPI object
  !> @param [in]     twod_mesh  2D mesh on which fields are defined (W3)
  !>
  subroutine define_partitions( self, mpi, twod_mesh )
    implicit none

    class(coupling_type), intent(inout)     :: self
    type(mpi_type), intent(inout)           :: mpi
    type( mesh_type ),  intent(in), pointer :: twod_mesh

#ifdef MCT
    ! Function space of fields used in coupling
    type(function_space_type), pointer          :: cpl_fs
    ! Pointer to the global indices from the mesh
    integer(i_halo_index), pointer              :: global_index_ptr(:)
    ! Global index for the first mesh level
    integer(i_def), allocatable                 :: global_index(:)
    ! Vector describing the 2d local grid partition in the global index space
    integer(i_def), allocatable                 :: ig_paral_2d(:)
    ! Vector describing the 0d local grid partition
    integer(i_def)                              :: ig_paral_0d(3)
    ! Returned OASIS error code
    integer(i_def)                              :: kinfo
    ! Loop index
    integer(i_def)                              :: i

    ! Coupling only works for 2d, zeroth-order, cell-centred fields
    ! - so create a function space for such a field
    if (twod_mesh%get_nlayers() > 1) then
      write(log_scratch_space,'(2A)') "define_coupling_partitions:", &
        "Currently, coupling only supports 2D meshes"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif
    cpl_fs => function_space_collection%get_fs( twod_mesh, 0, 0, W3 )
    self%cpl_size = cpl_fs%get_last_dof_owned()

    allocate(global_index(self%cpl_size))
    allocate(self%local_index(self%cpl_size))

    global_index_ptr => cpl_fs%get_global_dof_id()

    ! Convert global indices to integers, if possible
    if (maxval(global_index) > int(huge(i_def), i_halo_index)) then
      write(log_scratch_space,'(3A)') "define_coupling_partitions: ", &
         "Too many points for the coupler to deal with"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    else
      global_index(1:self%cpl_size) = &
                      int(global_index_ptr(1:self%cpl_size), i_def)
    endif

    do i = 1, self%cpl_size
      self%local_index(i) = i
    enddo

    ! Lookup used to sort 2d field indices to improve OASIS performance
    call bubble_sort(self%cpl_size, global_index, self%local_index)

    ! Set up a partition for 2d coupling
    allocate(ig_paral_2d(2+self%cpl_size))
    ig_paral_2d(1) = 4
    ig_paral_2d(2) = self%cpl_size

    do i = 1, self%cpl_size
      ig_paral_2d(i + 2) = global_index(i) + 1
    enddo

    deallocate(global_index)

    call oasis_def_partition (self%part_2d_id, ig_paral_2d, kinfo)

    !Set up a partition for 0d coupling (only couple 0d data (scalars) from PE0)
    ig_paral_0d(1)=0
    ig_paral_0d(2)=0
    if (mpi%get_comm_rank() == 0 ) then
      ig_paral_0d(3)=1
    else
      ig_paral_0d(3)=0
    endif
    call oasis_def_partition (self%part_0d_id, ig_paral_0d, kinfo)

    deallocate(ig_paral_2d)
#else
    write(log_scratch_space, * ) &
         "define_partitions: to use OASIS, cpp directive MCT must be set"
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif

  end subroutine define_partitions


  !>@brief Defines the variables that will be set to/received from the coupler
  !>
  !> @param [in,out] cpl_snd_2d field collection with fields to send
  !> @param [in,out] cpl_rcv_2d field collection with fields to receive
  !> @param [in,out] cpl_snd_0d field collection with fields to send as scalars
  !
  subroutine define_variables( self, cpl_snd_2d, cpl_rcv_2d, cpl_snd_0d )
    implicit none

    class(coupling_type), intent(inout)         :: self
    type( field_collection_type ), intent(inout):: cpl_snd_2d
    type( field_collection_type ), intent(inout):: cpl_rcv_2d
    type( field_collection_type ), intent(inout):: cpl_snd_0d

#ifdef MCT
    ! Iterator over field collection
    type( field_collection_iterator_type)       :: iter
    ! Pointer to a abstrct field parent type
    class( field_parent_type ), pointer         :: field_iter
    ! Name for transient fields (receive)
    character(str_def)                          :: var_name
    ! Name of field with level information for transient field (receive)
    character(str_def)                          :: var_name_lev
    ! Function space of fields used in coupling
    type(function_space_type), pointer          :: cpl_fs
    ! Number of multi-data fields
    integer(i_def)                              :: ndata
    ! Index for different do loops
    integer(i_def)                              :: i
    ! Name of the level of a multi-category field
    character(len=2)                            :: cpl_catno
    ! Rank/bundle information
    integer(i_def)                              :: var_nodims(2)
    ! Dimension of 2d coupled fields
    integer(i_def)                              :: var_shape_2d(2)
    ! Dimension of 2d coupled fields
    integer(i_def)                              :: var_shape_0d(1)
    ! Id for transient fields (receive)
    integer(i_def)                              :: var_id
    ! Number of coupling components the data will be sent to
    integer(i_def)                              :: ncpl
    ! Error code returned by oasis routine
    integer(i_def)                              :: kinfo
    ! Coupling frequency of each model
    integer(i_def)                              :: cpl_freqs(nmax)

    var_nodims(1) = 1 ! rank of coupling field
    var_nodims(2) = 1 ! number of bundles in coupling field (always 1)
    var_shape_2d(1) = 1
    var_shape_2d(2) = self%cpl_size
    var_shape_0d(1) = 1

    call iter%initialise(cpl_rcv_2d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      var_name = trim(adjustl(field_iter%get_name()))
      cpl_fs => field_iter%get_function_space()
      ndata = cpl_fs%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
          write(cpl_catno, cpl_fmt) i
          var_name_lev = trim(var_name)//cpl_cat//cpl_catno
          call oasis_def_var( var_id,             &
                              trim(var_name_lev), &
                              self%part_2d_id,    &
                              var_nodims,         &
                              oasis_in,           &
                              var_shape_2d,       &
                              prism_real,         &
                              kinfo)
          call field_iter%set_cpl_id(var_id, i)

          write(log_scratch_space, '(A)' ) &
                    "cpl_define: field "//trim(var_name_lev)//" receive"
          call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        enddo
      else
        call oasis_def_var( var_id,          &
                            trim(var_name),  &
                            self%part_2d_id, &
                            var_nodims,      &
                            oasis_in,        &
                            var_shape_2d,    &
                            prism_real,      &
                            kinfo)
        call field_iter%set_cpl_id(var_id, 1)

        write(log_scratch_space, '(A)' ) &
                     "cpl_define: field "//trim(var_name)//" receive"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      endif
      field_iter   => null()
    end do

    call iter%initialise(cpl_snd_2d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      var_name = trim(adjustl(field_iter%get_name()))
      cpl_fs => field_iter%get_function_space()
      ndata = cpl_fs%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
          write(cpl_catno, cpl_fmt) i
          var_name_lev = trim(var_name)//cpl_cat//cpl_catno
          call oasis_def_var( var_id,             &
                              trim(var_name_lev), &
                              self%part_2d_id,    &
                              var_nodims,         &
                              oasis_out,          &
                              var_shape_2d,       &
                              prism_real,         &
                              kinfo)
          call field_iter%set_cpl_id(var_id, i)

          write(log_scratch_space, '(A)' ) &
                       "cpl_define: field "//trim(var_name_lev)//" send"
          call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        enddo
      else
        call oasis_def_var( var_id,          &
                            trim(var_name),  &
                            self%part_2d_id, &
                            var_nodims,      &
                            oasis_out,       &
                            var_shape_2d,    &
                            prism_real,      &
                            kinfo)
        call field_iter%set_cpl_id(var_id, 1)

        write(log_scratch_space, '(A)' ) &
                          "cpl_define: field "//trim(var_name)//" send"
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      endif
      field_iter   => null()
    end do

    call iter%initialise(cpl_snd_0d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      var_name     = trim(adjustl(field_iter%get_name()))
      call oasis_def_var( var_id,          &
                          trim(var_name),  &
                          self%part_0d_id, &
                          var_nodims,      &
                          oasis_out,       &
                          var_shape_0d,    &
                          prism_real,      &
                          kinfo)
      call field_iter%set_cpl_id(var_id, 1)

      write(log_scratch_space, '(A)' ) &
                       "cpl_define: field "//trim(var_name)//" send"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
    end do

    call oasis_enddef (kinfo)

    ! Check that each field has the same the coupling frequency for
    ! all components the field is sent to
    call iter%initialise(cpl_snd_2d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      var_name = trim(adjustl(field_iter%get_name()))
      cpl_fs => field_iter%get_function_space()
      ndata = cpl_fs%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
          write(cpl_catno, cpl_fmt) i
          var_name_lev = trim(var_name)//cpl_cat//cpl_catno
          var_id = field_iter%get_cpl_id(i)
          call oasis_get_ncpl(var_id, ncpl, kinfo)
          call oasis_get_freqs(var_id, oasis_out, ncpl, &
                               cpl_freqs(1:ncpl), kinfo)
          if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
            write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
                   trim(var_name_lev),                                 &
                   " has different coupling frequencies for different components"
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          endif
        enddo
      else
        var_id = field_iter%get_cpl_id(1)
        call oasis_get_ncpl(var_id, ncpl, kinfo)
        call oasis_get_freqs(var_id, oasis_out, ncpl, &
                             cpl_freqs(1:ncpl), kinfo)
        if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
          write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
                 trim(var_name),                                     &
                 " has different coupling frequencies for different components"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
      endif
    end do
    call iter%initialise(cpl_snd_0d)
    do
      if (.not.iter%has_next())exit
      field_iter => iter%next()
      var_name     = trim(adjustl(field_iter%get_name()))
      var_id = field_iter%get_cpl_id(1)

      call oasis_get_ncpl(var_id, ncpl, kinfo)
      call oasis_get_freqs(var_id, oasis_out, ncpl, &
                           cpl_freqs(1:ncpl), kinfo)
      if (maxval(cpl_freqs(1:ncpl)) /= minval(cpl_freqs(1:ncpl))) then
        write(log_scratch_space, '(3A)' ) "ERROR: coupling field ", &
               trim(var_name),                                     &
               " has different coupling frequencies for different components"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      endif
    end do
#else
    write(log_scratch_space, * ) &
               "define_variables: to use OASIS, cpp directive MCT must be set"
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif

  end subroutine define_variables


  function get_local_index(self) result(local_index)
    implicit none
    class(coupling_type), intent(inout), target :: self

    integer(i_def), pointer             :: local_index(:)

    local_index => self%local_index

  end function get_local_index


  !> @brief Finalises coupler
  !
  subroutine finalise(self)
    implicit none
    class(coupling_type), intent(inout) :: self
#ifdef MCT
    integer(i_def) :: kinfo           ! error flag from OASIS

    if(allocated(self%local_index)) deallocate(self%local_index)

    kinfo = prism_ok
    call oasis_terminate(kinfo)
    if (kinfo .NE. prism_ok) then
      write(log_scratch_space,'(A, I4)') &
          "finalise: oasis_terminate error: ", kinfo
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      call oasis_abort(self%comp_id, 'finalise','abort1')
    else
      write(log_scratch_space,'(A)') "finalise : oasis_terminated OK"
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
    endif
#else
    write(log_scratch_space, * ) &
          "finalise: to use OASIS, cpp directive MCT must be set"
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif

  end subroutine finalise

  !-----------------------------------------------------------------------------
  ! Non-type-bound helper functions
  !-----------------------------------------------------------------------------

  !> @brief Helper function returns lists of fields to be sent to and
  !>        received from the coupler - from the coupler configuration
  subroutine get_coupling_fields(send_field_names, recv_field_names)
    implicit none
    character(str_def), allocatable, intent(out) :: send_field_names(:)
    character(str_def), allocatable, intent(out) :: recv_field_names(:)
#ifdef MCT
    integer(i_def)                               :: nfield

    allocate(send_field_names(nnamcpl))
    allocate(recv_field_names(nnamcpl))

    do nfield=1,nnamcpl
      send_field_names(nfield) = trim(adjustl(namsrcfld(nfield)))
      recv_field_names(nfield) = trim(adjustl(namdstfld(nfield)))
    end do
#else
    write(log_scratch_space, * ) &
          "get_coupling_fields: to use OASIS, cpp directive MCT must be set"
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif
  end subroutine get_coupling_fields

  !-----------------------------------------------------------------------------
  !> @brief Helper function to extract a concrete coupling object from a
  !>        key-value collection
  !> @param[in] collection The key-value collection to extract from
  !> @param[in] name       The name of the coupling object to extract
  !> @return    coupling   The requested coupling object
  function get_coupling_from_collection(collection, name) result(coupling)

  implicit none

    type(key_value_collection_type), intent(in) :: collection
    character(*),                    intent(in) :: name

    type(coupling_type), pointer        :: coupling

    class(abstract_value_type), pointer :: abstract_value

    call collection%get_value(trim(name), abstract_value)
    select type(abstract_value)
      type is (coupling_type)
        coupling => abstract_value
      class default
        write(log_scratch_space, * ) &
          "Error: the value called "//trim(name)//" is not a coupling object"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

  end function get_coupling_from_collection

end module coupling_mod
