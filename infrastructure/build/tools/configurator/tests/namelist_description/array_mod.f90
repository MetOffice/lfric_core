!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the aerial namelist.
!>
module aerial_config_mod

  use constants_mod, only: i_def, &
                           r_def, &
                           str_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key
  use wibble_mod, only: esize

  implicit none

  private
  public :: read_aerial_namelist, postprocess_aerial_namelist, &
            aerial_is_loadable, aerial_is_loaded, &
            aerial_reset_load_status, &
            aerial_multiples_allowed, aerial_final, &
            get_aerial_nml, get_new_aerial_nml

  integer(i_def), parameter, public :: max_array_size = 500

  character(str_def), public, protected :: absolute(5) = cmdi
  integer(i_def), public, protected, allocatable :: inlist(:)
  integer(i_def), public, protected :: lsize = imdi
  real(r_def), public, protected, allocatable :: outlist(:)
  integer(i_def), public, protected, allocatable :: unknown(:)

  character(*), parameter :: listname = 'aerial'
  character(str_def) :: profile_name = cmdi

  logical, parameter :: multiples_allowed = .false.

  logical :: nml_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !> @param [in] scan .true. if reading namelist to acquire scalar
  !>                  values which may possbly be required for
  !>                  array sizing during postprocessing.
  !>
  subroutine read_aerial_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan )

  end subroutine read_aerial_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    integer(i_def) :: buffer_integer_i_def(1)

    namelist /aerial/ absolute, &
                      inlist, &
                      lsize, &
                      outlist, &
                      unknown

    integer(i_def) :: condition

    if (allocated(inlist)) deallocate(inlist)
    allocate( inlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    if (allocated(outlist)) deallocate(outlist)
    allocate( outlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    if (allocated(unknown)) deallocate(unknown)
    allocate( unknown(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    absolute = cmdi
    inlist = imdi
    lsize = imdi
    outlist = rmdi
    unknown = imdi

    if (local_rank == 0) then

      read( file_unit, nml=aerial, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_def(1) = lsize

    call global_mpi%broadcast( buffer_integer_i_def, 1, 0 )

    lsize = buffer_integer_i_def(1)

    call global_mpi%broadcast( absolute, size(absolute, 1)*str_def, 0 )
    call global_mpi%broadcast( inlist, size(inlist, 1), 0 )
    call global_mpi%broadcast( outlist, size(outlist, 1), 0 )
    call global_mpi%broadcast( unknown, size(unknown, 1), 0 )

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_aerial_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(5)

    call members(1)%initialise( &
                'absolute', absolute )

    call members(2)%initialise( &
                'inlist', inlist )

    call members(3)%initialise( &
                'lsize', lsize )

    call members(4)%initialise( &
                'outlist', outlist )

    call members(5)%initialise( &
                'unknown', unknown )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_aerial_nml

  !> @brief Returns a <<aerial_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<aerial_nml_type>> with current namelist contents.
  function get_new_aerial_nml() result(namelist_obj)

    use aerial_nml_mod, only: aerial_nml_type

    implicit none

    type(aerial_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(5)

    call members(1)%initialise( &
                'absolute', absolute )

    call members(2)%initialise( &
                'inlist', inlist )

    call members(3)%initialise( &
                'lsize', lsize )

    call members(4)%initialise( &
                'outlist', outlist )

    call members(5)%initialise( &
                'unknown', unknown )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_aerial_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_aerial_namelist()

    use constants_mod, only: i_def

    implicit none

    integer(i_def) :: condition
    integer(i_def) :: array_size


    integer(i_def), allocatable :: new_inlist(:)
    real(r_def), allocatable :: new_outlist(:)
    integer(i_def), allocatable :: new_unknown(:)
    integer(i_def) :: index_unknown

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !
    ! Arrays are re-sized to fit data.
    !
    condition  = 0
    array_size = 0


    array_size = lsize
    if (array_size == imdi) then
      write(log_scratch_space, '(A)') &
          '"aerial:inlist" not allocated, '// &
          'deferred size "lsize" '//   &
          'has not been specified.'
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      array_size = 0
    end if
    allocate( new_inlist(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_inlist(:array_size) = inlist(:array_size)
    call move_alloc( new_inlist, inlist )
    if (allocated(new_inlist)) deallocate( new_inlist)

    array_size = esize
    if (array_size == imdi) then
      write(log_scratch_space, '(A)') &
          '"aerial:outlist" not allocated, '// &
          'deferred size "esize" '//   &
          'has not been specified.'
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      array_size = 0
    end if
    allocate( new_outlist(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_outlist(:array_size) = outlist(:array_size)
    call move_alloc( new_outlist, outlist )
    if (allocated(new_outlist)) deallocate( new_outlist)

    do index_unknown=ubound(unknown, 1), 1, -1
      if (unknown(index_unknown) /= imdi) exit
    end do
    array_size = index_unknown
    allocate( new_unknown(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_unknown(:array_size) = unknown(:array_size)
    call move_alloc( new_unknown, unknown )
    if (allocated(new_unknown)) deallocate( new_unknown)


  end subroutine postprocess_aerial_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function aerial_is_loadable()

    implicit none

    logical :: aerial_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      aerial_is_loadable = .true.
    else
      aerial_is_loadable = .false.
    end if

  end function aerial_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function aerial_is_loaded()

    implicit none

    logical :: aerial_is_loaded

    aerial_is_loaded = nml_loaded

  end function aerial_is_loaded

  !> Are multiple aerial namelists allowed to be read?
  !>
  !> @return True If multiple aerial namelists are
  !>              permitted.
  !>
  function aerial_multiples_allowed()

    implicit none

    logical :: aerial_multiples_allowed

    aerial_multiples_allowed = multiples_allowed

  end function aerial_multiples_allowed

  !> Resets the load status to allow
  !> aerial namelist to be read.
  !>
  subroutine aerial_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine aerial_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine aerial_final()

    implicit none

    absolute = cmdi
    lsize = imdi

    if ( allocated(inlist) ) deallocate(inlist)
    if ( allocated(outlist) ) deallocate(outlist)
    if ( allocated(unknown) ) deallocate(unknown)

    return
  end subroutine aerial_final


end module aerial_config_mod
