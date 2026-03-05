!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the cheese namelist.
!>
module cheese_config_mod

  use constants_mod, only: i_def, &
                           r_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, FUDGE, imdi, rmdi, str_def, unset_key

  implicit none

  private
  public :: read_cheese_namelist, postprocess_cheese_namelist, &
            cheese_is_loadable, cheese_is_loaded, &
            cheese_reset_load_status, &
            cheese_multiples_allowed, cheese_final, &
            get_cheese_nml, get_new_cheese_nml

  real(r_def), public, protected :: fred = rmdi
  real(r_def), public, protected :: wilma = rmdi

  character(*), parameter :: listname = 'cheese'
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
  subroutine read_cheese_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan )

  end subroutine read_cheese_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    real(r_def) :: buffer_real_r_def(1)

    namelist /cheese/ fred

    integer(i_def) :: condition

    fred = rmdi
    wilma = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=cheese, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = fred

    call global_mpi%broadcast( buffer_real_r_def, 1, 0 )

    fred = buffer_real_r_def(1)

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_cheese_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(2)

    call members(1)%initialise( &
                'fred', fred )

    call members(2)%initialise( &
                'wilma', wilma )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_cheese_nml

  !> @brief Returns a <<cheese_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<cheese_nml_type>> with current namelist contents.
  function get_new_cheese_nml() result(namelist_obj)

    use cheese_nml_mod, only: cheese_nml_type

    implicit none

    type(cheese_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(2)

    call members(1)%initialise( &
                'fred', fred )

    call members(2)%initialise( &
                'wilma', wilma )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_cheese_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_cheese_namelist()

    use constants_mod, only: i_def

    implicit none

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !! Parameter name wilma: derived by computation
    wilma = fred * FUDGE

  end subroutine postprocess_cheese_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function cheese_is_loadable()

    implicit none

    logical :: cheese_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      cheese_is_loadable = .true.
    else
      cheese_is_loadable = .false.
    end if

  end function cheese_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function cheese_is_loaded()

    implicit none

    logical :: cheese_is_loaded

    cheese_is_loaded = nml_loaded

  end function cheese_is_loaded

  !> Are multiple cheese namelists allowed to be read?
  !>
  !> @return True If multiple cheese namelists are
  !>              permitted.
  !>
  function cheese_multiples_allowed()

    implicit none

    logical :: cheese_multiples_allowed

    cheese_multiples_allowed = multiples_allowed

  end function cheese_multiples_allowed

  !> Resets the load status to allow
  !> cheese namelist to be read.
  !>
  subroutine cheese_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine cheese_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine cheese_final()

    implicit none

    fred = real(rmdi,r_def)
    wilma = real(rmdi,r_def)

    return
  end subroutine cheese_final


end module cheese_config_mod
