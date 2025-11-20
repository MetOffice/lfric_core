!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Collection of small utility procedures for lfric-xios
!>
module lfric_xios_utils_mod

  use constants_mod,        only: i_def, r_def, str_def, str_long
  use file_mod,             only: FILE_OP_OPEN, FILE_MODE_READ
  use lfric_ncdf_dims_mod,  only: lfric_ncdf_dims_type
  use lfric_ncdf_field_mod, only: lfric_ncdf_field_type
  use lfric_ncdf_file_mod,  only: lfric_ncdf_file_type
  use lfric_mpi_mod,        only: global_mpi
  use lfric_xios_field_mod, only: lfric_xios_field_type
  use log_mod,              only: log_event, log_scratch_space, &
                                  LOG_LEVEL_ERROR, LOG_LEVEL_INFO, &
                                  LOG_LEVEL_TRACE
  use mesh_mod,             only: mesh_type
  use xios,                 only: xios_date, xios_duration,        &
                                  xios_get_time_origin,            &
                                  xios_get_year_length_in_seconds, &
                                  xios_date_convert_to_seconds,    &
                                  operator(<), operator(+)


  implicit none
  private
  public :: parse_date_as_xios, seconds_from_date, &
            set_prime_io_mesh, prime_io_mesh_is,   &
            read_time_data

  integer(i_def), private, allocatable :: prime_io_mesh_ids(:)

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Interpret a string as an XIOS date object.
  !> Expected format is yyyy-mm-dd hh:mm:ss
  !>
  !> @param [in] date_str The string representation of the date
  !> @result     date_obj The xios_date object represented by the input string
  !>
  function parse_date_as_xios( date_str ) result( date_obj )
    implicit none
    character(*), intent(in) :: date_str
    type(xios_date)          :: date_obj
    integer(i_def)           :: y, mo, d, h, mi, s, size

    size = len(date_str)

    ! Indexing from end to support arbitrarily long year
    read( date_str(1      :size-15), * ) y
    read( date_str(size-13:size-12), * ) mo
    read( date_str(size-10:size-9 ), * ) d
    read( date_str(size-7 :size-6 ), * ) h
    read( date_str(size-4 :size-3 ), * ) mi
    read( date_str(size-1 :size   ), * ) s

    date_obj = xios_date( y, mo, d, h, mi, s )

  end function parse_date_as_xios

  !> @brief  Wrapper around xios_date_convert_to_seconds to nullify XIOS bug.
  !>
  !> @param[in] date            The input xios_date object
  !> @result    date_in_seconds The resulting seconds converted from "date"
  function seconds_from_date(date) result(date_in_seconds)

    implicit none

    type(xios_date), intent(in) :: date

    type(xios_date)        :: time_origin
    integer(i_def)         :: year_diff
    real(r_def)            :: date_in_seconds

    integer(i_def), parameter :: length360d = 31104000

    call xios_get_time_origin(time_origin)

    ! Get time in seconds from XIOS dates - due to a bug in XIOS, non-360day
    ! calendars do not return the correct values around the time origin so a
    ! workaround is implemented below. Also due to a bug in XIOS calendar types
    ! cannot be identified except by the number of seconds per year
    if ( date < time_origin .and. &
          xios_get_year_length_in_seconds(date%year) /= length360d ) &
      then
      year_diff = date%year - time_origin%year
      date_in_seconds = real(xios_date_convert_to_seconds(date), r_def) + &
        ( real(xios_get_year_length_in_seconds(date%year), r_def) * &
          real(year_diff, r_def) )
    else
      date_in_seconds = real(xios_date_convert_to_seconds(date), r_def)
    end if

  end function seconds_from_date

  !> @brief Registers a mesh to be used as the primary I/O mesh
  !> @param[in] mesh  The mesh object to be registered
  subroutine set_prime_io_mesh( mesh )

    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(i_def), allocatable :: mesh_id_list(:)

    ! Set up array of ints to hold mesh ids, bring in previous mesh IDs if
    ! already present and add the new mesh ID to the new array
    if (allocated(prime_io_mesh_ids)) then
      allocate(mesh_id_list(size(prime_io_mesh_ids) + 1))
      mesh_id_list(1:size(prime_io_mesh_ids)) = prime_io_mesh_ids
      mesh_id_list(size(prime_io_mesh_ids)+1) = mesh%get_id()
    else
      allocate(mesh_id_list(1))
      mesh_id_list(1) = mesh%get_id()
    end if

    ! Make the new array the main array
    call move_alloc(mesh_id_list, prime_io_mesh_ids)

  end subroutine set_prime_io_mesh

  function prime_io_mesh_is( mesh ) result( mesh_is_prime_io )

    implicit none

    type(mesh_type), intent(in) :: mesh

    logical :: mesh_is_prime_io

    if (.not. allocated(prime_io_mesh_ids) ) then
      mesh_is_prime_io = .false.
    else
      mesh_is_prime_io = (any(prime_io_mesh_ids == mesh%get_id()))
    end if

  end function prime_io_mesh_is


  !>  @brief  Read time data into time axis using NetCDF tools
  !!
  !> @param[in] file_path  The path to the file containing the data controlled
  !!                       by this object
  function read_time_data(file_path) result(time_data)

    implicit none

    character(len=*),      intent(in)   :: file_path

    type(xios_date), allocatable :: time_data(:)

    ! Local variables for XIOS interface
    integer(i_def)                :: n_t, t, n_t_buf(1)
    real(r_def), allocatable      :: input_data(:)
    character(str_def)            :: time_units, ref_date_str, var_id, dim_id
    character(str_def)            :: time_meta_buf(2)
    character(str_long)           :: unit_attr
    character(str_def), parameter :: valid_units(4) = &
                                  (/'seconds', 'days   ', 'hours  ', 'months '/)
    type(xios_duration)           :: ref_time, mean_time, month_duration
    type(xios_date)               :: ref_date
    integer(i_def),  parameter :: len_date = 18 ! Length of CF date is 18 characters
    integer(i_def),  parameter :: len_delim = 7 ! Delimiter between unit and
                                                    ! date is 7 characters

    ! NetCDF reading variables
    type(lfric_ncdf_file_type)  :: file_ncdf
    type(lfric_ncdf_dims_type)  :: time_dim
    type(lfric_ncdf_field_type) :: time_var

    if (global_mpi%get_comm_rank() == 0) then
      file_ncdf = lfric_ncdf_file_type( trim(file_path)//".nc", &
                                        open_mode=FILE_OP_OPEN, &
                                        io_mode=FILE_MODE_READ )

      ! Some JULES surface ancils have non-CF-compliant time representation
      ! so we need to account for that for the time being
      if (file_ncdf%contains_var("time")) then
        var_id = "time"
        dim_id = "time"
      else if (file_ncdf%contains_var("month_number")) then
        var_id = "month_number"
        dim_id = "month_number"
      else
        call log_event( "Invalid representation of time in file ["// &
                        trim(file_path)//"]", log_level_error)
      end if

      ! Get size of time axis from file
      time_dim = lfric_ncdf_dims_type(trim(dim_id), file_ncdf)
      n_t = time_dim%get_size()
      allocate( input_data( n_t ) )
      allocate( time_data( n_t ) )

      ! Read the time data from the ancil file
      time_var = lfric_ncdf_field_type(trim(var_id), file_ncdf)
      call time_var%read_data(input_data)

      ! Read time units and reference date from file
      if (trim(var_id) == "time") then
        unit_attr = time_var%get_char_attribute("units")
        time_units = unit_attr( 1 : len(trim(unit_attr))-len_date-len_delim )
        if ( .not. any( valid_units == trim(adjustl(time_units)) ) ) then
          write( log_scratch_space,'(A,A)' ) "Invalid units of ["//trim(time_units)// &
                                            "] for time axis in file: "// trim(file_path)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        ref_date_str = unit_attr( len(trim(unit_attr))-len_date : len(trim(unit_attr)) )
      else if (trim(var_id) == "month_number") then
        ! Non CF files are given the bare minimum treatment
        time_units = "months"
        ref_date_str = "1970-01-01 00:00:00"
      end if

      call file_ncdf%close_file()

    end if ! comm_rank == 0

    ! Broadcast data across MPI ranks
    n_t_buf(1) = n_t
    time_meta_buf = [ref_date_str, time_units]
    call global_mpi%broadcast(n_t_buf, 1, 0)
    if (global_mpi%get_comm_rank() /= 0) allocate(input_data(n_t_buf(1)))
    call global_mpi%broadcast(input_data, n_t_buf(1), 0)
    call global_mpi%broadcast(time_meta_buf, size(time_meta_buf, 1)*str_def, 0)

    ref_date_str = time_meta_buf(1)
    time_units = time_meta_buf(2)
    ref_date = parse_date_as_xios(trim(adjustl(ref_date_str)))

    ! Convert input time data to xios_date type
    do t = 1, size(input_data)
      ref_time = xios_duration(0, 0, 0, 0, 0, 0)
      if ( trim(adjustl(time_units)) == "seconds" ) then
        ref_time%second = input_data(t)
      else if ( trim(adjustl(time_units)) == "hours" ) then
        ref_time%hour = input_data(t)
      else if ( trim(adjustl(time_units)) == "days" ) then
        ref_time%day = input_data(t)
      else if ( trim(adjustl(time_units)) == "months" ) then
        ! Offset months backwards to account for monthly mean
        ref_time%month = input_data(t) - 1
      end if
      time_data(t) = ref_date + ref_time
    end do

    ! Correct "months" data to be monthly mean - centred on middle of month.
    ! This can't be done above as the conversion to seconds can only be done
    ! with xios_date objects, not xios_durations
    if ( time_units == "months" ) then
      mean_time = xios_duration(0, 0, 0, 0, 0, 0)
      month_duration = xios_duration(0, 1, 0, 0, 0, 0)
      do t = 1, n_t_buf(1)
        mean_time%second = ( &
          seconds_from_date(time_data(t)+month_duration) - &
          seconds_from_date(time_data(t)) ) / 2
        time_data(t) = time_data(t) + mean_time
      end do
    end if

  end function read_time_data

end module lfric_xios_utils_mod