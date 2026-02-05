!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the mirth namelist.
!>
module mirth_config_mod

  use constants_mod, only: i_def, &
                           str_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key
  use random_config_mod, only: biggles

  implicit none

  private
  public :: read_mirth_namelist, postprocess_mirth_namelist, &
            mirth_is_loadable, mirth_is_loaded, &
            mirth_reset_load_status, &
            mirth_multiples_allowed, mirth_final, &
            get_mirth_nml, get_new_mirth_nml

  integer(i_def), parameter, public :: max_array_size = 500

  character(str_def), public, protected, allocatable :: chortle(:)
  character(str_def), public, protected :: chuckle = cmdi
  character(str_def), public, protected :: guffaw(3) = cmdi
  character(str_def), public, protected, allocatable :: hysterics(:)

  character(*), parameter :: listname = 'mirth'
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
  subroutine read_mirth_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan )

  end subroutine read_mirth_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    character(str_def) :: buffer_character_str_def(1)

    namelist /mirth/ chortle, &
                     chuckle, &
                     guffaw, &
                     hysterics

    integer(i_def) :: condition

    if (allocated(hysterics)) deallocate(hysterics)
    allocate( hysterics(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "hysterics"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    if (allocated(chortle)) deallocate(chortle)
    allocate( chortle(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "chortle"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    chortle = cmdi
    chuckle = cmdi
    guffaw = cmdi
    hysterics = cmdi

    if (local_rank == 0) then

      read( file_unit, nml=mirth, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_character_str_def(1) = chuckle

    call global_mpi%broadcast( buffer_character_str_def, 1*str_def, 0 )

    chuckle = buffer_character_str_def(1)

    call global_mpi%broadcast( chortle, size(chortle, 1)*str_def, 0 )
    call global_mpi%broadcast( guffaw, size(guffaw, 1)*str_def, 0 )
    call global_mpi%broadcast( hysterics, size(hysterics, 1)*str_def, 0 )

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_mirth_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(4)

    call members(1)%initialise( &
                'chortle', chortle )

    call members(2)%initialise( &
                'chuckle', chuckle )

    call members(3)%initialise( &
                'guffaw', guffaw )

    call members(4)%initialise( &
                'hysterics', hysterics )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_mirth_nml

  !> @brief Returns a <<mirth_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<mirth_nml_type>> with current namelist contents.
  function get_new_mirth_nml() result(namelist_obj)

    use mirth_nml_mod, only: mirth_nml_type

    implicit none

    type(mirth_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(4)

    call members(1)%initialise( &
                'chortle', chortle )

    call members(2)%initialise( &
                'chuckle', chuckle )

    call members(3)%initialise( &
                'guffaw', guffaw )

    call members(4)%initialise( &
                'hysterics', hysterics )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_mirth_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_mirth_namelist()

    use constants_mod, only: i_def

    implicit none

    integer(i_def) :: condition
    integer(i_def) :: array_size


    character(str_def), allocatable :: new_hysterics(:)
    integer(i_def) :: index_hysterics
    character(str_def), allocatable :: new_chortle(:)

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !
    ! Arrays are re-sized to fit data.
    !
    condition  = 0
    array_size = 0


    do index_hysterics=ubound(hysterics, 1), 1, -1
      if (hysterics(index_hysterics) /= cmdi) exit
    end do
    array_size = index_hysterics
    allocate( new_hysterics(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "hysterics"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_hysterics(:array_size) = hysterics(:array_size)
    call move_alloc( new_hysterics, hysterics )
    if (allocated(new_hysterics)) deallocate( new_hysterics)

    array_size = biggles
    if (array_size == imdi) then
      write(log_scratch_space, '(A)') &
          '"mirth:chortle" not allocated, '// &
          'deferred size "biggles" '//   &
          'has not been specified.'
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      array_size = 0
    end if
    allocate( new_chortle(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "chortle"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_chortle(:array_size) = chortle(:array_size)
    call move_alloc( new_chortle, chortle )
    if (allocated(new_chortle)) deallocate( new_chortle)


  end subroutine postprocess_mirth_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function mirth_is_loadable()

    implicit none

    logical :: mirth_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      mirth_is_loadable = .true.
    else
      mirth_is_loadable = .false.
    end if

  end function mirth_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function mirth_is_loaded()

    implicit none

    logical :: mirth_is_loaded

    mirth_is_loaded = nml_loaded

  end function mirth_is_loaded

  !> Are multiple mirth namelists allowed to be read?
  !>
  !> @return True If multiple mirth namelists are
  !>              permitted.
  !>
  function mirth_multiples_allowed()

    implicit none

    logical :: mirth_multiples_allowed

    mirth_multiples_allowed = multiples_allowed

  end function mirth_multiples_allowed

  !> Resets the load status to allow
  !> mirth namelist to be read.
  !>
  subroutine mirth_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine mirth_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine mirth_final()

    implicit none

    chuckle = cmdi
    guffaw = cmdi

    if ( allocated(chortle) ) deallocate(chortle)
    if ( allocated(hysterics) ) deallocate(hysterics)

    return
  end subroutine mirth_final


end module mirth_config_mod
