!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> A module containing utility functions to contain events to be attached to
!> the ticking of a model clock
!>
module event_mod

  use clock_mod,     only : clock_type
  use constants_mod, only : i_def

  implicit none

  private

  ! Abstract class used for timestep event interface
  type, public, abstract :: event_actor_type
  contains
  end type event_actor_type

  ! Abstract interface for timestep event actions
  abstract interface
    subroutine event_action( actor, clock )
      import event_actor_type, clock_type
      class(event_actor_type), intent(inout) :: actor
      class(clock_type),       intent(in)    :: clock
    end subroutine event_action
  end interface

  ! Type containing an event that will be attached to the model clock's tick
  type, public :: event_type
    private
    procedure(event_action), nopass, pointer :: action => null()
    class(event_actor_type),         pointer :: actor  => null()
    logical                                  :: active = .false.
  contains
    procedure, public :: is_active
    procedure, public :: happens
  end type event_type

  interface event_type
    procedure event_constructor
  end interface event_type

  public :: event_action

contains

  !> Constructs an event object
  !>
  !> @param[in] input_action The procedure that will be called during the event
  !> @param[in] input_actor  The object that will perform the action
  function event_constructor(input_action, input_actor) result(new_event)

    implicit none

    procedure(event_action), pointer, intent(in) :: input_action
    class(event_actor_type), target,  intent(in) :: input_actor
    type(event_type) :: new_event

    new_event%action => input_action
    new_event%actor => input_actor
    new_event%active = .true.

  end function event_constructor

  !> Returns true if the event has been constructed and is ready to run
  function is_active(this) result(l_active)

    implicit none

    class(event_type), intent(inout) :: this
    logical                          :: l_active

    l_active = this%active

  end function is_active

  !> Runs an instance of the event
  !>
  !> @param[in] clock The model clock
  subroutine happens(this, clock)

    implicit none

    class(event_type), intent(inout) :: this
    class(clock_type), intent(in)    :: clock

    call this%action( this%actor, clock )

  end subroutine happens

end module event_mod
