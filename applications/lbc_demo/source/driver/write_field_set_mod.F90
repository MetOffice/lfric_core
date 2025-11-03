 !-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module write_field_set_mod

  use constants_mod,                 only: str_def
  use field_collection_mod,          only: field_collection_type
  use field_collection_iterator_mod, only: field_collection_iterator_type
  use field_mod,                     only: field_type
  use field_parent_mod,              only: field_parent_type
  use integer_field_mod,             only: integer_field_type

  use log_mod, only: log_event,         &
                     log_scratch_space, &
                     log_level_info

  implicit none

  private

  public :: write_field_set

contains
!===============================================================================

subroutine write_field_set( fld_collection, id_suffix )

  implicit none

  type(field_collection_type), pointer,  intent(inout) :: fld_collection
  character(str_def),          optional, intent(in)    :: id_suffix

  type( field_collection_iterator_type) :: iter

  class(field_parent_type), pointer :: field
  type(field_type),         pointer :: tmp_real_field
  type(integer_field_type), pointer :: tmp_int_field

  character(str_def) :: name
  character(str_def) :: suffix

  if (present(id_suffix)) then
    suffix = trim(id_suffix)
  else
    suffix = ''
  end if

  nullify(field)
  call iter%initialise( fld_collection )

  do

    if ( iter%has_next() ) then
      field => iter%next()
    else
      exit
    end if

    select type (field)

    type is (field_type)
      tmp_real_field => field
      name = trim(tmp_real_field%get_name())//trim(suffix)
      call tmp_real_field%write_field(trim(name))

    type is (integer_field_type)
      tmp_int_field => field
      name = trim(tmp_int_field%get_name())//trim(suffix)
      call tmp_int_field%write_field(trim(name))

    end select

    write(log_scratch_space,'(A)') &
        'Writing field to '//trim(name)
    call log_event( log_scratch_space, log_level_info )

  end do

end subroutine write_field_set

end module write_field_set_mod
