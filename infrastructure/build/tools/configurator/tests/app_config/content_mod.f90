!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   Defines \<config_type\> object.
!> @details A container object that holds namelist
!>          objects of various types).
!>
!>          Access pattern will differ for namelist types that are permitted
!>          to have multiple instances within the configuration.
!>
module config_mod

  use constants_mod,   only: i_def, l_def, str_def, cmdi
  use log_mod,         only: log_event, log_scratch_space, &
                             log_level_error, log_level_warning
  use linked_list_mod, only: linked_list_type, linked_list_item_type

  use namelist_mod,            only: namelist_type
  use namelist_collection_mod, only: namelist_collection_type

  use foo_nml_mod, only: foo_nml_type
  use bar_nml_mod, only: bar_nml_type
  use moo_nml_mod, only: moo_nml_type
  use pot_nml_mod, only: pot_nml_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that stores namelists of an application configuration
  !-----------------------------------------------------------------------------
  type, public :: config_type

    private

    !> The name of the namelist collection if provided.
    character(:), allocatable :: config_name

    !> Whether object has been initialised or not
    logical :: isinitialised = .false.

    !> The name of the namelist collection if provided.
    character(str_def), allocatable :: nml_fullnames(:)

    ! Single instance namelists
    type(foo_nml_type), public, allocatable :: foo
    type(moo_nml_type), public, allocatable :: moo

    ! Namelists which may have multiple instances.
    ! These are accesed via the associated
    ! <namelist name>_list methods.
    type(linked_list_type), public, allocatable :: bar
    type(linked_list_type), public, allocatable :: pot

  contains

    procedure, public :: initialise
    procedure, public :: name
    procedure, public :: add_namelist
    procedure, public :: contents
    procedure, public :: n_namelists
    procedure, public :: namelist_exists

    procedure, public :: bar_list
    procedure, public :: pot_list

    procedure, public :: clear

    final :: config_destructor

    procedure, private :: update_contents

  end type config_type

contains


!> @brief Initialises application configuration.
!> @param [in] name  Optional: The name given to the configuration.
!=====================================================================
subroutine initialise(self, name)

  implicit none

  class(config_type), intent(inout) :: self

  character(*), optional, intent(in) :: name

  if (self%isinitialised) then
    write(log_scratch_space, '(A)')       &
        'Application configuration: [' // &
         trim(self%config_name)        // &
        '] has already been initiaised.'
    call log_event(log_scratch_space, log_level_error)
  end if

  if (present(name)) then
    self%config_name = trim(name)
  else
    self%config_name = cmdi
  end if

  self%isinitialised = .true.

end subroutine initialise


!> @brief Installs a new namelist object into the configuration.
!> @param [in] namelist_obj The extended namelist type object. Only
!>                          extended namelist types defined by the
!>                          application metadata file will be accepted.
!===================================================================
subroutine add_namelist(self, namelist_obj)

  implicit none

  class(config_type), intent(inout) :: self

  class(namelist_type), intent(in) :: namelist_obj

  character(:), allocatable :: name
  character(:), allocatable :: profile_name
  character(:), allocatable :: full_name

  ! Check namelist name is valid, if not then exit with error
  full_name    = namelist_obj%get_full_name()
  profile_name = namelist_obj%get_profile_name()
  name         = namelist_obj%get_listname()

  select type(namelist_obj)

  type is( foo_nml_type )
    ! Multiple instances: NOT ALLOWED
    if (self%namelist_exists(trim(name))) then
      write(log_scratch_space, '(A)') &
          trim(name) // ' namelist already allocated.'
      call log_event(log_scratch_space, log_level_error)
    else
      allocate(self%foo, source=namelist_obj)
      call self%update_contents(trim(name))
    end if

  type is( moo_nml_type )
    ! Multiple instances: NOT ALLOWED
    if (self%namelist_exists(trim(name))) then
      write(log_scratch_space, '(A)') &
          trim(name) // ' namelist already allocated.'
      call log_event(log_scratch_space, log_level_error)
    else
      allocate(self%moo, source=namelist_obj)
      call self%update_contents(trim(name))
    end if

  type is ( bar_nml_type )
    ! Multiple instances: ALLOWED
    if (trim(profile_name) == cmdi) then
      write(log_scratch_space, '(A)') 'Ignoring ' // trim(name) // &
          ' namelist: missing profile name.'
      call log_event(log_scratch_space, log_level_warning)
    else if (self%namelist_exists(trim(full_name))) then
      write(log_scratch_space, '(A)') trim(name) // &
          ' namelist (' // trim(profile_name) // '), already allocated.'
      call log_event(log_scratch_space, log_level_error)
    else
      if (.not. allocated(self%bar)) then
        allocate(self%bar)
      end if
      call self%bar%insert_item( namelist_obj )
      call self%update_contents(namelist_obj%get_full_name())
    end if

  type is ( pot_nml_type )
    ! Multiple instances: ALLOWED
    if (trim(profile_name) == cmdi) then
      write(log_scratch_space, '(A)') 'Ignoring ' // trim(name) // &
          ' namelist: missing profile name.'
      call log_event(log_scratch_space, log_level_warning)
    else if (self%namelist_exists(trim(full_name))) then
      write(log_scratch_space, '(A)') trim(name) // &
          ' namelist (' // trim(profile_name) // '), already allocated.'
      call log_event(log_scratch_space, log_level_error)
    else
      if (.not. allocated(self%pot)) then
        allocate(self%pot)
      end if
      call self%pot%insert_item( namelist_obj )
      call self%update_contents(namelist_obj%get_full_name())
    end if

  class default
     write(log_scratch_space, '(A)')                  &
         ' Undefined namelist type(' // trim(name) // &
         '), for this configuration.'
     call log_event(log_scratch_space, log_level_error)

  end select

end subroutine add_namelist


!> @brief Check if a namelist is present the collection.
!> @param [in] name         The name of the namelist to be checked.
!> @param [in] profile_name Optional: In the case of namelists which
!>                          are permitted to have multiple instances,
!>                          the profile name distiguishes the instances
!>                          of <name> namelists.
!> @return exists           Flag stating if namelist is present or not
!=====================================================================
function namelist_exists(self, name, profile_name) result(exists)

  implicit none

  class(config_type), intent(in) :: self

  character(*),           intent(in) :: name
  character(*), optional, intent(in) :: profile_name

  logical(l_def) :: exists

  integer(i_def)     :: i
  character(str_def) :: full_name

  exists = .false.

  if (allocated(self%nml_fullnames)) then

    if (present(profile_name)) then
      full_name = trim(name)//':'//trim(profile_name)
    else
      full_name = trim(name)
    end if

    do i=1, size(self%nml_fullnames)
      if (trim(self%nml_fullnames(i)) == trim(full_name)) then
        exists = .true.
        exit
      end if
    end do
  end if

end function namelist_exists


!> @brief Returns a pointer to an instance of <bar_nml_type>.
!> @param [in] profile_name  Profile name used to identify the
!>                           instance of <bar_nml_type>.
!> @return bar_nml_obj  Pointer to the requested namelist object.
!=====================================================================
function bar_list(self, profile_name) result(bar_nml_obj)

  implicit none

  class(config_type), intent(in) :: self
  character(*),       intent(in) :: profile_name

  type(bar_nml_type), pointer :: bar_nml_obj

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop
  character(str_def) :: payload_name

  nullify(bar_nml_obj)
  nullify(loop)

  loop => self%bar%get_head()
  do
    ! If the list is empty or the end of the list was
    ! reached without finding the namelist, fail with
    ! an error.
    if (.not. associated(loop)) then
      write(log_scratch_space, '(A)') &
          'Instance ' // trim(profile_name) // ' of ' // &
          'bar_nml_type ' // &
          'not found in configuration.'
      call log_event(log_scratch_space, log_level_error)
    end if

    ! Otherwise 'cast' to a bar_namelist_type
    select type(payload => loop%payload)
    type is (bar_nml_type)
      payload_name = payload%get_profile_name()
      if (trim(profile_name) == trim(payload_name)) then
        bar_nml_obj => payload
        exit
      end if
    end select

    loop => loop%next
  end do

end function bar_list

!> @brief Returns a pointer to an instance of <pot_nml_type>.
!> @param [in] profile_name  Profile name used to identify the
!>                           instance of <pot_nml_type>.
!> @return pot_nml_obj  Pointer to the requested namelist object.
!=====================================================================
function pot_list(self, profile_name) result(pot_nml_obj)

  implicit none

  class(config_type), intent(in) :: self
  character(*),       intent(in) :: profile_name

  type(pot_nml_type), pointer :: pot_nml_obj

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop
  character(str_def) :: payload_name

  nullify(pot_nml_obj)
  nullify(loop)

  loop => self%pot%get_head()
  do
    ! If the list is empty or the end of the list was
    ! reached without finding the namelist, fail with
    ! an error.
    if (.not. associated(loop)) then
      write(log_scratch_space, '(A)') &
          'Instance ' // trim(profile_name) // ' of ' // &
          'pot_nml_type ' // &
          'not found in configuration.'
      call log_event(log_scratch_space, log_level_error)
    end if

    ! Otherwise 'cast' to a pot_namelist_type
    select type(payload => loop%payload)
    type is (pot_nml_type)
      payload_name = payload%get_profile_name()
      if (trim(profile_name) == trim(payload_name)) then
        pot_nml_obj => payload
        exit
      end if
    end select

    loop => loop%next
  end do

end function pot_list

!> @brief Queries config_type for the total number of namelists stored.
!> @return answer  The number of namelists stored
!=====================================================================
function n_namelists(self) result(answer)

  implicit none

  class(config_type), intent(in) :: self

  integer(i_def) :: answer

  answer = 0
  if (allocated(self%nml_fullnames)) then
    answer = size(self%nml_fullnames)
  end if

end function n_namelists

!> @brief Queries the name of config_type.
!> @return name  The name identifying this namelist collection
!>               on initialisation.
!=====================================================================
function name(self) result(answer)

  implicit none

  class(config_type), intent(in) :: self
  character(:), allocatable :: answer

  answer = self%config_name

end function name

!> @brief Extracts namelist names in config_type.
!> @param  listname        Optional: if specified, returns entries
!>                         begining with this string.
!> @return namelist_names  Array of unique names of namelists in the
!>                         collection.
!=====================================================================
function contents(self, listname) result(namelist_names)

  implicit none

  class(config_type), intent(in) :: self

  character(*), optional, intent(in) :: listname

  character(str_def), allocatable :: namelist_names(:)

  character(str_def), allocatable :: tmp(:)
  character(str_def) :: tmp_str
  integer(i_def) :: n_found, i, start_index

  if (allocated(namelist_names)) deallocate(namelist_names)

  n_found = 0
  if (present(listname)) then

    allocate(tmp(size(self%nml_fullnames)))

    do i=1, size(self%nml_fullnames)
      if (index(trim(self%nml_fullnames(i)), trim(listname)) > 0) then
        tmp_str      = trim(self%nml_fullnames(i))
        start_index  = index(tmp_str, ':')
        n_found      = n_found + 1_i_def
        tmp(n_found) = trim(tmp_str(start_index+1:))
      end if
    end do

    allocate(namelist_names(n_found))
    namelist_names = tmp(1:n_found)
    deallocate(tmp)

  else

    allocate(namelist_names, source=self%nml_fullnames)

  end if

end function contents


!> @brief Clears all items from the namelist collection.
!=====================================================================
subroutine clear(self)

  implicit none

  class(config_type), intent(inout) :: self

  ! Namlists which may have multiple instances per configuration
  if (allocated(self%bar)) call self%bar%clear()
  if (allocated(self%pot)) call self%pot%clear()

  if (allocated(self%foo)) deallocate(self%foo)
  if (allocated(self%bar)) deallocate(self%bar)
  if (allocated(self%moo)) deallocate(self%moo)
  if (allocated(self%pot)) deallocate(self%pot)

  if (allocated(self%nml_fullnames)) deallocate(self%nml_fullnames)

  self%config_name = cmdi
  self%isinitialised = .false.

end subroutine clear


!> @brief Destructor for the namelist collection
!=====================================================================
subroutine config_destructor(self)

  implicit none

  type(config_type), intent(inout) :: self

  call self%clear()

end subroutine config_destructor


!> @brief Adds namelist identifier to the to list on namelists stored.
!> @param [in] nml_full_name  Namelists identifier to be added.
!=====================================================================
subroutine update_contents(self, nml_full_name)

  implicit none

  class(config_type), intent(inout) :: self

  character(*), intent(in) :: nml_full_name

  character(str_def), allocatable :: tmp_str(:)
  integer(i_def) :: n_entries

  if (allocated(self%nml_fullnames)) then

    n_entries = size(self%nml_fullnames)
    allocate(tmp_str, source=self%nml_fullnames)
    deallocate(self%nml_fullnames)
    allocate(self%nml_fullnames(n_entries+1))
    self%nml_fullnames(1:n_entries) = tmp_str(:)
    self%nml_fullnames(n_entries+1) = nml_full_name

  else

    allocate(self%nml_fullnames(1))
    self%nml_fullnames(1) = trim(nml_full_name)

  end if

end subroutine update_contents

end module config_mod
