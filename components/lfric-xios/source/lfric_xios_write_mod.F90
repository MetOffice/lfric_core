!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!>  @brief    Module for field writing routines
!>  @details  Holds all routines for writing LFRic fields
!>
module lfric_xios_write_mod

  use clock_mod,            only: clock_type
  use constants_mod,        only: i_def, l_def, r_def, str_def, r_second, &
                                  str_max_filename, EPS
  use lfric_xios_constants_mod, &
                            only: dp_xios, xios_max_int
  use linked_list_mod,      only: linked_list_item_type
  use field_real32_mod,     only: field_real32_type, field_real32_proxy_type
  use field_real64_mod,     only: field_real64_type, field_real64_proxy_type
  use io_value_mod,         only: io_value_type
  use key_value_mod,        only: key_value_type, abstract_key_value_type, &
                                  abstract_value_type
  use key_value_collection_mod, &
                            only: key_value_collection_type
  use key_value_collection_iterator_mod, &
                            only: key_value_collection_iterator_type
  use field_parent_mod,     only: field_parent_proxy_type
  use field_collection_iterator_mod, &
                            only: field_collection_iterator_type
  use field_collection_mod, only: field_collection_type
  use field_parent_mod,     only: field_parent_type
  use fs_continuity_mod,    only: W3
  use io_mod,               only: ts_fname
  use integer_field_mod,    only: integer_field_type, integer_field_proxy_type
  use lfric_xios_utils_mod, only: prime_io_mesh_is
  use lfric_xios_format_mod, &
                            only: format_field
  use mesh_mod,             only: mesh_type
  use model_clock_mod,      only: model_clock_type
  use log_mod,              only: log_event,         &
                                  log_scratch_space, &
                                  LOG_LEVEL_INFO,    &
                                  LOG_LEVEL_DEBUG,   &
                                  LOG_LEVEL_WARNING, &
                                  LOG_LEVEL_ERROR
  use lfric_string_mod,     only: split_string
  use timing_mod,           only: start_timing, stop_timing, &
                                  tik, LPROF
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,  only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr,   &
                                  xios_is_valid_field,  &
                                  lfric_xios_mock_pull_in
#else
  use lfric_xios_mock_mod,  only: lfric_xios_mock_pull_in
  use xios,                 only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr,   &
                                  xios_is_valid_field
#endif

  implicit none

  private
  public :: checkpoint_write_xios,    &
            write_field_generic,      &
            write_empty_field,        &
            checkpoint_write_value,   &
            write_value_generic,      &
            write_state,              &
            write_checkpoint,         &
            create_checkpoint_list

contains

!> @brief Write io_value data via XIOS
!> @details This routine assumes there is a XIOS field defined
!>          with a field id the same as the io_value id
!> @param[in,out] io_value The io_value to write data from
!>
subroutine write_value_generic(io_value, value_name)
  class(io_value_type), intent(inout) :: io_value
  character(*), optional, intent(in) :: value_name

  integer(i_def) :: array_dims
  character(:), allocatable :: value_id

  if (present(value_name)) then
    value_id = value_name
  else
    value_id = io_value%io_id
  end if

  array_dims = size(io_value%data)
  if ( xios_is_valid_field(trim(value_id)) ) then
    call xios_send_field( trim(value_id), &
                          reshape(io_value%data, (/ 1, array_dims /)) )
  else
    call log_event( 'No XIOS field with id="'//trim(io_value%io_id)//'" is defined', &
                    LOG_LEVEL_ERROR )
  end if

end subroutine write_value_generic

!>  @brief  Write field data to UGRIDs via XIOS
!>
!>  @param[in]     field_name       Field name (for error reporting only)
!>  @param[in]     field_proxy      A field proxy to be written
!>
subroutine write_field_generic(field_name, field_proxy)
  use lfric_xios_diag_mod,        only:  get_field_domain_ref, field_is_active
  implicit none

  character(len=*), intent(in) :: field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: undf
  integer(i_def) :: hdim          ! horizontal dimension, domain size
  integer(i_def) :: vdim          ! vertical dimension
  real(dp_xios), allocatable :: xios_data(:)
  logical(l_def) :: legacy
  integer(tik)   :: timing_id

  ! If the field is not active in xios at this timestep, exit this routine
  ! without doing anything
  if (.not. field_is_active(field_name, .true.)) return

  if ( LPROF ) call start_timing(timing_id, 'lfric_xios.write_fldg')

  undf = field_proxy%vspace%get_last_dof_owned() ! total dimension

  vdim = field_proxy%vspace%get_ndata() * size(field_proxy%vspace%get_levels())

  hdim = undf/vdim

  ! detect field with legacy checkpointing domain
  legacy = (index(get_field_domain_ref(field_name), 'checkpoint_') == 1)

  ! sanity check
  if (.not. legacy .and. .not. (hdim*vdim == undf)) then
    call log_event('assertion failed for field ' // field_name                &
      // ': hdim*vdim == undf', log_level_error)
  end if

  allocate(xios_data(undf))

  call format_field(xios_data, field_name, field_proxy, vdim, hdim, legacy)

  if (legacy) then
    call xios_send_field( field_name, reshape(xios_data, (/ 1, undf /)) )
  else
    call xios_send_field( field_name, reshape(xios_data, (/vdim, hdim/)) )
    ! The shape is only necessary for the mock implementation, and
    ! the only thing that matters is the product of the dimensions.
  end if

  deallocate(xios_data)

  if ( LPROF ) call stop_timing(timing_id, 'lfric_xios.write_fldg')

end subroutine write_field_generic

!>  @brief  Graceful failure if an empty field is attempted to be written
!>
!>  @param[in]     field_name       Field name
!>  @param[in]     field_proxy      A field proxy to be written
!>
subroutine write_empty_field(field_name, field_proxy)
  implicit none

  character(len=*), intent(in) :: field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  ! Note that this routine simply outputs an informative warning.
  ! Future versions may force an error by logging to LOG_LEVEL_ERROR.
  write(log_scratch_space,'(2A)') &
        "Attempt to write an empty field: ", field_name
  call log_event(log_scratch_space, LOG_LEVEL_WARNING)


end subroutine write_empty_field

!> @brief Checkpoint an io_value with XIOS
!> @details This routine assumes there is an XIOS field
!>          with the "checkpoint_" prefix
!> @param[in,out] io_value The io_value to write data from
!>
subroutine checkpoint_write_value(io_value, value_name)
  class(io_value_type), intent(inout) :: io_value
  character(*), optional, intent(in) :: value_name

  character(str_def) :: checkpoint_id
  integer(i_def)     :: array_dims

  if(present(value_name)) then
    checkpoint_id = trim(value_name)
  else
    checkpoint_id = trim(io_value%io_id)
  end if
  array_dims = size(io_value%data)
  if ( xios_is_valid_field(trim(checkpoint_id)) ) then
    call xios_send_field( trim(checkpoint_id), &
                          reshape(io_value%data, (/ 1, array_dims /)) )
  else
    call log_event( 'No XIOS field with id="'//trim(checkpoint_id)//'" is defined', &
                    LOG_LEVEL_ERROR )
  end if

end subroutine checkpoint_write_value

!>  @brief    I/O handler for writing an XIOS netcdf checkpoint
!>  @details  Note this routine accepts a filename but doesn't use it - this is
!>            to keep the interface the same for all methods
!>
!>  @param[in]      xios_field_name  XIOS identifier for the field
!>  @param[in]      file_name        Name of the file to write into
!>  @param[in,out]  field_proxy      A field proxy to be written
!>
subroutine checkpoint_write_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  character(len=*),               intent(in) :: file_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def)             :: undf
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  allocate(send_field(undf))

  ! Different field kinds are selected to access data
  select type(field_proxy)

    type is (field_real32_proxy_type)
    send_field = field_proxy%data(1:undf)

    type is (field_real64_proxy_type)
    send_field = field_proxy%data(1:undf)

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data(1:undf)) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    send_field = real( field_proxy%data(1:undf), dp_xios )

    class default
    call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  call xios_send_field("checkpoint_"//trim(xios_field_name), reshape (send_field, (/1, undf/)))

end subroutine checkpoint_write_xios

!>  @brief    Write a collection of fields
!>  @details  Iterate over a field collection and write each field if it is
!>            enabled for writing
!>
!>  @param[in]           state   A collection of fields
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine write_state(state, prefix, suffix)

  implicit none

  type(field_collection_type), intent(inout) :: state
  character(len=*), optional,  intent(in)    :: prefix
  character(len=*), optional,  intent(in)    :: suffix

  type(field_collection_iterator_type) :: iter
  character(str_def)                   :: xios_field_id

  class(field_parent_type), pointer :: fld => null()

  ! Create the iter iterator on the state collection
  call iter%initialise(state)
  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_real32_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if
      type is (field_real64_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if
      type is (integer_field_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if

    end select
  end do

  nullify(fld)

end subroutine write_state

!>  @brief    Write a checkpoint from a collection of fields
!>  @details  Iterate over a field collection and checkpoint each field
!>            if it is enabled for checkpointing
!>
!>  @param[in]  fields  Fields to checkpoint.
!>  @param[in]  values  Values to checkpoint.
!>  @param[in]  clock  Model time
!>  @param[in]  checkpoint_stem_name  The checkpoint file stem name
!>  @param[in]  checkpoint_times The checkpoint times
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine write_checkpoint( fields, values, clock, checkpoint_stem_name, &
                             checkpoint_times, prefix, suffix )

  implicit none

  type(field_collection_type), intent(inout) :: fields
  type(key_value_collection_type), intent(inout) :: values
  type(model_clock_type),      intent(in)    :: clock
  character(len=*),            intent(in)    :: checkpoint_stem_name
  real(r_second),              intent(in)    :: checkpoint_times(:)
  character(len=*), optional,  intent(in)    :: prefix
  character(len=*), optional,  intent(in)    :: suffix

  type(field_collection_iterator_type) :: iter
  type(key_value_collection_iterator_type) :: val_iter
  class(key_value_type), pointer :: kv
  class(abstract_value_type), pointer :: abstract_val

  class(field_parent_type), pointer    :: fld => null()

  character(str_def)                   :: xios_field_id
  character(str_def)                   :: field_prefix
  character(:), allocatable            :: split_stem_name(:)

  if(checkpoint_time(clock, checkpoint_times)) then
    ! Create the field prefix from the checkpoint stem name and current time step
    split_stem_name = split_string( trim(checkpoint_stem_name), '/' )
    write(field_prefix,'(A,A,I10.10,A)') &
          trim(split_stem_name(size(split_stem_name))),"_", &
          clock%get_step(), "_"
    call iter%initialise(fields)
    do
       if ( .not.iter%has_next() ) exit
       fld => iter%next()
       ! Construct the XIOS field ID from the LFRic field name and optional arguments
       xios_field_id = trim(adjustl(fld%get_name()))
       if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
       if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))
       select type(fld)
       type is (field_real32_type)
          if ( fld%can_checkpoint() ) then
             write(log_scratch_space,'(2A)') &
                  "Checkpointing ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_checkpoint( xios_field_id,      &
                                        trim(ts_fname(checkpoint_stem_name, &
                                        "",                                 &
                                        xios_field_id,      &
                                        clock%get_step(),                   &
                                        "")) )
          else if ( fld%can_write() ) then
             write(log_scratch_space,'(2A)') &
                  "Writing checkpoint for ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_field( trim(field_prefix) // trim(xios_field_id) )
          else
             call log_event( 'Writing not set up for '// xios_field_id, &
                            LOG_LEVEL_INFO )
          end if
       type is (field_real64_type)
          if ( fld%can_checkpoint() ) then
             write(log_scratch_space,'(2A)') &
                  "Checkpointing ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_checkpoint( xios_field_id,      &
                                        trim(ts_fname(checkpoint_stem_name, &
                                        "",                                 &
                                        xios_field_id,      &
                                        clock%get_step(),                   &
                                        "")) )
          else if ( fld%can_write() ) then
             write(log_scratch_space,'(2A)') &
                  "Writing checkpoint for ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_field( trim(field_prefix) // trim(xios_field_id) )
          else
             call log_event( 'Writing not set up for '// xios_field_id, &
                            LOG_LEVEL_INFO )
          end if
       type is (integer_field_type)
          if ( fld%can_checkpoint() ) then
             write(log_scratch_space,'(2A)') &
                  "Checkpointing ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_checkpoint( trim(adjustl(fld%get_name()) ),     &
                                        trim(ts_fname(checkpoint_stem_name, &
                                        "",                                 &
                                        xios_field_id,      &
                                        clock%get_step(),                   &
                                        "")) )
          else if ( fld%can_write() ) then
             write(log_scratch_space,'(2A)') &
                  "Writing checkpoint for ", xios_field_id
             call log_event(log_scratch_space, LOG_LEVEL_INFO)
             call fld%write_field( trim(field_prefix) // trim(xios_field_id) )
          else
             call log_event( 'Writing not set up for '// xios_field_id, &
                  LOG_LEVEL_INFO )
          end if
       class default
          call log_event('write_checkpoint:Invalid type of field, not supported supported',LOG_LEVEL_ERROR)
       end select
    end do

    ! Loop over values collection and checkpoint as appropriate
    call val_iter%initialise(values)
    do
      if (.not. val_iter%has_next()) exit
      kv => val_iter%next()
      select type(kv_typed => kv)
        type is (abstract_key_value_type)
          abstract_val => kv_typed%value
          select type (io_value_object => abstract_val)
            type is (io_value_type)
              if(io_value_object%can_write_checkpoint()) then
                call log_event( 'Writing checkpoint for ' // &
                                trim(io_value_object%io_id), &
                                LOG_LEVEL_INFO )
                call io_value_object%write_checkpoint( &
                        trim(field_prefix) // trim(io_value_object%io_id))
              end if
          end select
      end select
    end do
  end if

  nullify(fld)
  nullify(abstract_val)
  nullify(kv)

end subroutine write_checkpoint

!> @brief Check if the current time is a checkpoint time
!> @param[in] clock The model clock
!> @param[in] checkpoint_times The checkpoint times
!> @return is_checkpoint_time Logical flag indicating if current time is a checkpoint time
function checkpoint_time(clock, checkpoint_times) result(is_checkpoint_time)
  implicit none
  type(model_clock_type), intent(in) :: clock
  real(r_second), intent(in) :: checkpoint_times(:)
  logical(l_def) :: is_checkpoint_time
  integer(i_def) :: i

  is_checkpoint_time = .false.
  do i = 1, size(checkpoint_times)
    if ( clock%get_step() == clock%steps_from_seconds(checkpoint_times(i)) ) then
      is_checkpoint_time = .true.
      exit
    end if
  end do

end function checkpoint_time

!> @brief Create a list of checkpoint times, adding the end time if required
!> @param[in] clock The model clock
!> @param[in] checkpoint_times_input The input list of checkpoint times
!> @param[in] lcheckpoint_end Logical flag to indicate whether to add the end time
!> @param[out] checkpoint_times_output The output list of checkpoint times
subroutine create_checkpoint_list(clock, &
                                  checkpoint_times_input, &
                                  lcheckpoint_end, &
                                  checkpoint_times_output)
  implicit none
  type(model_clock_type), intent(in) :: clock
  real(r_second), intent(in) :: checkpoint_times_input(:)
  logical, intent(in) :: lcheckpoint_end
  real(r_second), allocatable, intent(out) :: checkpoint_times_output(:)

  integer(i_def) :: number_of_checkpoint_times

  number_of_checkpoint_times = size(checkpoint_times_input)
  if(lcheckpoint_end) then
    ! Check if the end time is already in the list (within a tolerance)
    if (any(abs(checkpoint_times_input - clock%seconds_from_steps(clock%get_last_step())) <= EPS)) then
      allocate(checkpoint_times_output, source=checkpoint_times_input)

    else
      allocate(checkpoint_times_output(number_of_checkpoint_times + 1))
      checkpoint_times_output(1:number_of_checkpoint_times) = checkpoint_times_input
      checkpoint_times_output(number_of_checkpoint_times + 1) = &
        clock%seconds_from_steps(clock%get_last_step())
    end if
  else
    allocate(checkpoint_times_output, source=checkpoint_times_input)
  endif

end subroutine create_checkpoint_list

end module lfric_xios_write_mod
