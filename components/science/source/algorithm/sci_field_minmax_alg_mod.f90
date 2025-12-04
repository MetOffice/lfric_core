!-----------------------------------------------------------------------------
! Copyright (c) 2023,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of field min-max for real32, real64
!         and int32 fields

module sci_field_minmax_alg_mod
  use, intrinsic :: iso_fortran_env, only : real32, real64, int32
  use sci_psykal_builtin_light_mod,  only : invoke_real32_field_min_max,       &
                                            invoke_real64_field_min_max,       &
                                            invoke_int32_field_min_max,        &
                                            invoke_real32_local_field_min_max, &
                                            invoke_real64_local_field_min_max, &
                                            invoke_int32_local_field_min_max

  use log_mod,                       only : log_event, log_scratch_space, &
                                            log_at_level
  use field_real32_mod,              only : field_real32_type
  use field_real64_mod,              only : field_real64_type
  use field_int32_mod,               only : field_int32_type

  implicit none

  interface get_field_minmax
     module procedure get_field_minmax_real32, &
                      get_field_minmax_real64, &
                      get_field_minmax_int32
  end interface get_field_minmax

  interface get_local_field_minmax
     module procedure get_local_field_minmax_real32, &
                      get_local_field_minmax_real64, &
                      get_local_field_minmax_int32
  end interface get_local_field_minmax

  interface log_field_minmax
     module procedure log_field_minmax_real32, &
                      log_field_minmax_real64, &
                      log_field_minmax_int32
  end interface log_field_minmax

contains
  !> Returns minimum and maximum of the data values of a field
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_field_minmax_real32( field, fmin, fmax )
    implicit none
    type(field_real32_type), intent(in) :: field
    real(kind=real32),      intent(out) :: fmin, fmax

    ! call the invoke in the PSy layer
    call invoke_real32_field_min_max( fmin, fmax, field )
  end subroutine get_field_minmax_real32

  !> Returns minimum and maximum of the data values of a field
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_field_minmax_real64( field, fmin, fmax )
    implicit none
    type(field_real64_type), intent(in) :: field
    real(kind=real64),      intent(out) :: fmin, fmax

    ! call the invoke in the psy layer
    call invoke_real64_field_min_max( fmin, fmax, field )
  end subroutine get_field_minmax_real64

  !> Returns minimum and maximum of the data values of a field
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_field_minmax_int32( field, fmin, fmax )
    implicit none
    type(field_int32_type), intent(in)  :: field
    integer(kind=int32),    intent(out) :: fmin, fmax

    ! call the invoke in the psy layer
    call invoke_int32_field_min_max( fmin, fmax, field )
  end subroutine get_field_minmax_int32

  !> Returns minimum and maximum of the field data values local to the rank
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_local_field_minmax_real32( field, fmin, fmax )
    implicit none
    type(field_real32_type), intent(in)  :: field
    real(kind=real32),       intent(out) :: fmin, fmax

    ! call the invoke in the PSy layer
    call invoke_real32_local_field_min_max( fmin, fmax, field )
  end subroutine get_local_field_minmax_real32

  !> Returns minimum and maximum of the field data values local to the rank
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_local_field_minmax_real64( field, fmin, fmax )
    implicit none
    type(field_real64_type), intent(in)  :: field
    real(kind=real64),       intent(out) :: fmin, fmax

    ! call the invoke in the psy layer
    call invoke_real64_local_field_min_max( fmin, fmax, field )
  end subroutine get_local_field_minmax_real64

  !> Returns minimum and maximum of the field data values local to the rank
  !> @param[in] field The field for which the min and max are required
  !> @param[out] fmin The minimum of the field
  !> @param[out] fmax The minimum of the field
  subroutine get_local_field_minmax_int32( field, fmin, fmax )
    implicit none
    type(field_int32_type), intent(in)  :: field
    integer(kind=int32),    intent(out) :: fmin, fmax

    ! call the invoke in the psy layer
    call invoke_int32_local_field_min_max( fmin, fmax, field )
  end subroutine get_local_field_minmax_int32

  !> Logs the minimum and maximum values of a field to the log
  !> @param[in] log_level The logging level at which to write the message
  !> @param[in] label Text that will be written  along with the min/max
  !> @param[in]  field The field for which the min and max are required
  subroutine log_field_minmax_real32( log_level, label, field )
    integer,                     intent(in) :: log_level
    character(len = *),          intent(in) :: label
    type(field_real32_type),     intent(in) :: field
    real(kind=real32)                       :: fmin, fmax

    ! If we aren't going to log the min and max then we don't need to
    ! do any further work here.
    if ( .not. log_at_level(log_level) ) return

    ! Calculate min and max field values for real32 field
    call invoke_real32_field_min_max( fmin, fmax, field )
    write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
           "Min/max ", trim(label), " = ", fmin, fmax
    call log_event( log_scratch_space, log_level )
  end subroutine log_field_minmax_real32

  !> Logs the minimum and maximum values of a field to the log
  !> @param[in] log_level The logging level at which to write the message
  !> @param[in] label Text that will be written  along with the min/max
  !> @param[in]  field The field for which the min and max are required
  subroutine log_field_minmax_real64( log_level, label, field )
    integer,                     intent(in) :: log_level
    character(len = *),          intent(in) :: label
    type(field_real64_type),     intent(in) :: field
    real(kind=real64)                       :: fmin, fmax

    ! If we aren't going to log the min and max then we don't need to
    ! do any further work here.
    if ( .not. log_at_level(log_level) ) return

    ! Calculate min and max field values for real64 field
    call invoke_real64_field_min_max( fmin, fmax, field )
    write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
          "Min/max ", trim(label), " = ", fmin, fmax
    call log_event( log_scratch_space, log_level )
  end subroutine log_field_minmax_real64

  !> Logs the minimum and maximum values of a field to the log
  !> @param[in] log_level The logging level at which to write the message
  !> @param[in] label Text that will be written  along with the min/max
  !> @param[in]  field The field for which the min and max are required
  subroutine log_field_minmax_int32( log_level, label, field )
    integer,                     intent(in) :: log_level
    character(len = *),          intent(in) :: label
    type(field_int32_type),      intent(in) :: field
    integer(kind=int32)                     :: fmin, fmax

    ! If we aren't going to log the min and max then we don't need to
    ! do any further work here.
    if ( .not. log_at_level(log_level) ) return

    ! Calculate min and max field values for real32 field
    call invoke_int32_field_min_max( fmin, fmax, field )
    write( log_scratch_space, '( A, A, A, 2I16 )' ) &
           "Min/max ", trim(label), " = ", fmin, fmax
    call log_event( log_scratch_space, log_level )
  end subroutine log_field_minmax_int32

end module sci_field_minmax_alg_mod

