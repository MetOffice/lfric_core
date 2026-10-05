!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Provides functionality to output field values to the log. This will be
!>        used for debugging, so there may not be any usage from production
!>        code

module log_field_mod

  use abstract_external_field_mod, only: abstract_external_field_type
  use constants_mod,               only: i_def
  use field_mod,                   only: field_type, field_proxy_type
  use field_parent_mod,            only: field_parent_type
  use function_space_mod,          only: function_space_type
  use log_mod,                     only: log_event,       &
                                         log_level_trace, &
                                         log_level_error, &
                                         log_scratch_space
  implicit none

  private

  ! Private type can only be instantiated from within this module
  type, extends(abstract_external_field_type) :: log_field_external_type
    !> Log level at which to write the field data
    integer(i_def) :: log_level = log_level_trace
    !> First dof of the section to be logged
    integer(i_def) :: start_dof
    !> Last dof of the section to be logged
    integer(i_def) :: end_dof
  contains
    !> Initialises the object
    procedure, public :: initialise
    !> Set the logging level for the next message
    procedure, public :: set_log_level
    !> Copy data from the LFRic field and pass it to the logger
    procedure, public :: copy_from_lfric
    !> Dummy required by abstract - just produces an error message
    procedure, public :: copy_to_lfric
  end type  log_field_external_type

public log_field

contains

  !> @brief Initialises the external field for logging field data
  !> @param [in] lfric_field_ptr Pointer to an lfric field
  subroutine initialise( self, lfric_field_ptr, start_dof, end_dof )
  implicit none

  class(log_field_external_type),     intent(inout) :: self
  type(field_type), pointer, intent(in)    :: lfric_field_ptr
  integer(i_def), optional, intent(in) :: start_dof
  integer(i_def), optional, intent(in) :: end_dof

  type(function_space_type), pointer :: function_space

  function_space => lfric_field_ptr%get_function_space()
  self%end_dof = function_space%get_undf()
  self%start_dof = 1
  if(present(start_dof))self%start_dof = start_dof
  if(present(end_dof))self%end_dof = end_dof

  call self%abstract_external_field_initialiser(lfric_field_ptr)

  end subroutine initialise


  subroutine set_log_level(self, log_level)
  implicit none

  class(log_field_external_type), intent(inout) :: self
  integer(i_def),                 intent(in)    :: log_level

  self%log_level = log_level

  end subroutine set_log_level

  !>@brief Outputs field values to the log
  !>@param return_code Optional return code from the copy_from procedure
  subroutine copy_from_lfric(self, return_code)
  implicit none

  class(log_field_external_type), intent(inout) :: self
  integer(i_def), intent(out), optional :: return_code

  class(field_parent_type), pointer :: abstract_field
  type(field_type), pointer :: field
  type(field_proxy_type):: fieldp
  integer(i_def) :: df

  abstract_field => self%get_lfric_field_ptr()
  select type (abstract_field)
    type is (field_type)
    field => abstract_field
    class default
      call log_event( &
        'Failed to log field. Object being logged is not a real field', &
        log_level_error)
  end select
  fieldp = field%get_proxy()

  do df = self%start_dof, self%end_dof
    write( log_scratch_space, '( I6, E16.8 )' ) df,fieldp%data( df )
    call log_event( log_scratch_space, self%log_level )
  end do

  if(present(return_code))return_code = 0

  end subroutine copy_from_lfric


  !>@brief Dummy required by abstract - just produced an error message
  !>@param return_code The return code from the copy_to procedure
  subroutine copy_to_lfric( self, return_code )
  implicit none

  class(log_field_external_type), intent(inout) :: self
  integer(i_def), intent(out), optional :: return_code

  call log_event( "ERROR: log_field_external_type has no copy_to_lfric functionality", &
                  log_level_error )

  if(present(return_code))return_code = 1

  end subroutine copy_to_lfric


  !>@brief Public wrapper to hide all the mechanics of the external field.
  !>       The user can then simply call this subroutine to log field values
  !>@param field The field to output dofs from
  !>@param log_level Optional logging level to log the output to.
  !>                 If not present the output is logged at log_level_trace
  !>@param label Optional label to write to the log before the data
  !>@param start_dof If only a section of the field is to be logged, this
  !>                 is the first dof of the section to be logged
  !>@param end_dof   If only a section of the field is to be logged, this
  !>                 is the last dof of the section to be logged
  subroutine log_field(field, log_level, label, start_dof, end_dof)
  implicit none

  type(field_type), target, intent(in) :: field
  integer(i_def), optional, intent(in) :: log_level
  character(*),   optional, intent(in) :: label
  integer(i_def), optional, intent(in) :: start_dof
  integer(i_def), optional, intent(in) :: end_dof
  type(field_type), pointer            :: fieldp
  type(log_field_external_type)        :: field_log
  integer(i_def)                       :: level

  if(present(log_level))then
    level = log_level
  else
    level = log_level_trace
  end if
  if(present(label)) call log_event( label, level )

  fieldp => field
  call field_log%initialise(fieldp, start_dof, end_dof)
  call field_log%set_log_level(level)
  call field_log%copy_from_lfric()

  end subroutine log_field

end module log_field_mod
