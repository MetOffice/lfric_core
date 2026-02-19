!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the teapot namelist.
!>
module teapot_config_mod

  use constants_mod, only: i_def, &
                           r_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key

  implicit none

  private
  public :: read_teapot_namelist, postprocess_teapot_namelist, &
            teapot_is_loadable, teapot_is_loaded, &
            teapot_reset_load_status, &
            teapot_multiples_allowed, teapot_final, &
            get_teapot_nml, get_new_teapot_nml

  real(r_def), public, protected :: bar = rmdi
  real(r_def), public, protected :: foo = rmdi
  real(r_def), public, protected :: fum = rmdi

  character(*), parameter :: listname = 'teapot'
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
  subroutine read_teapot_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan )

  end subroutine read_teapot_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    real(r_def) :: buffer_real_r_def(2)

    namelist /teapot/ foo, &
                      fum

    integer(i_def) :: condition

    bar = rmdi
    foo = rmdi
    fum = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=teapot, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = foo
    buffer_real_r_def(2) = fum

    call global_mpi%broadcast( buffer_real_r_def, 2, 0 )

    foo = buffer_real_r_def(1)
    fum = buffer_real_r_def(2)

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_teapot_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(3)

    call members(1)%initialise( &
                'bar', bar )

    call members(2)%initialise( &
                'foo', foo )

    call members(3)%initialise( &
                'fum', fum )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_teapot_nml

  !> @brief Returns a <<teapot_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<teapot_nml_type>> with current namelist contents.
  function get_new_teapot_nml() result(namelist_obj)

    use teapot_nml_mod, only: teapot_nml_type

    implicit none

    type(teapot_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(3)

    call members(1)%initialise( &
                'bar', bar )

    call members(2)%initialise( &
                'foo', foo )

    call members(3)%initialise( &
                'fum', fum )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_teapot_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_teapot_namelist()

    use constants_mod, only: i_def

    implicit none

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !! Parameter name bar: derived by computation
    bar = foo ** 2

  end subroutine postprocess_teapot_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function teapot_is_loadable()

    implicit none

    logical :: teapot_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      teapot_is_loadable = .true.
    else
      teapot_is_loadable = .false.
    end if

  end function teapot_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function teapot_is_loaded()

    implicit none

    logical :: teapot_is_loaded

    teapot_is_loaded = nml_loaded

  end function teapot_is_loaded

  !> Are multiple teapot namelists allowed to be read?
  !>
  !> @return True If multiple teapot namelists are
  !>              permitted.
  !>
  function teapot_multiples_allowed()

    implicit none

    logical :: teapot_multiples_allowed

    teapot_multiples_allowed = multiples_allowed

  end function teapot_multiples_allowed

  !> Resets the load status to allow
  !> teapot namelist to be read.
  !>
  subroutine teapot_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine teapot_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine teapot_final()

    implicit none

    bar = real(rmdi,r_def)
    foo = real(rmdi,r_def)
    fum = real(rmdi,r_def)

    return
  end subroutine teapot_final


end module teapot_config_mod
