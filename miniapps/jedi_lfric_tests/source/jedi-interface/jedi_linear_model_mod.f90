!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the linear model (tlm)
!>
!> @details This module includes a class that handles the linear models time
!>          stepping. The class includes init, step and final methods for the
!>          Tangent Linear (TL) and Adjoint (AD). These methods are used by the
!>          forecastTL and forecastAD methods also included within the class.
!>          The forecast methods require a linear-state produced by pre-running
!>          the non-linear model and storing  the result in a trajectory
!>          object. The set_trajectory method is included to provide the means
!>          to create and populate the linear state fields. In JEDI, the
!>          forecast* methods are defined in the OOPS base class and the init,
!>          step and final are defined in the model interface (LFRIC-JEDI). An
!>          included forecast application (jedi_tlm_forecast_tl) uses the model
!>          forecastTL method to propagate the increment.
!>
module jedi_linear_model_mod

  use constants_mod,                 only : i_def, str_def
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_state_mod,                only : jedi_state_type
  use jedi_increment_mod,            only : jedi_increment_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use linear_state_trajectory_mod, &
                                     only : linear_state_trajectory_type
  use field_collection_mod,          only : field_collection_type

  implicit none

  private

type, public :: jedi_linear_model_type
  private

  !> The model time step duration
  type ( jedi_duration_type )          :: time_step

  !> Trajectory of linear states obtained by running the non-linear model
  type( linear_state_trajectory_type ) :: linear_state_trajectory

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, public  :: set_trajectory

  procedure, private :: model_initTL
  procedure, private :: model_stepTL
  procedure, private :: model_finalTL

  procedure, private :: model_initAD
  procedure, private :: model_stepAD
  procedure, private :: model_finalAD

  !> Run a TL and AD forecasts
  procedure, public  :: forecastTL
  procedure, public  :: forecastAD

  !> Finalizer
  final              :: jedi_linear_model_destructor

end type jedi_linear_model_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_linear_model_type
!>
!> @param [in] config  The linear model configuration
subroutine initialise( self, config )

  use jedi_lfric_fake_nl_driver_mod, only : mesh
  use jedi_linear_model_config_mod,  only : jedi_linear_model_config_type

  implicit none

  class( jedi_linear_model_type ),    intent(inout) :: self
  type( jedi_linear_model_config_type ), intent(in) :: config

  self%time_step = config%time_step
  call self%linear_state_trajectory%initialise( config%forecast_length, &
                                                config%time_step,       &
                                                config%field_meta_data, &
                                                mesh )
end subroutine initialise

!> @brief    Set an instance of the trajectory
!>
!> @param [in] jedi_state The state to add to the trajectory
subroutine set_trajectory( self, jedi_state )

  implicit none

  class( jedi_linear_model_type ),   intent(inout) :: self
  type( jedi_state_type ),           intent(inout) :: jedi_state

  ! Local
  type( field_collection_type ), pointer :: next_linear_state
  character( len=str_def ), allocatable  :: variable_names(:)

  ! Create a new linear state and return a pointer to it
  call self%linear_state_trajectory%create_next_linear_state( &
                                              jedi_state%valid_time(), &
                                              next_linear_state )

  ! Get the names of the variables to copy
  call self%linear_state_trajectory%get_variable_names(variable_names)

  ! Copy data from the input state into the model_data
  call jedi_state%to_lfric_field_collection( variable_names, &
                                             next_linear_state )

end subroutine set_trajectory

!> @brief    Initialise the TL model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initTL(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(in) :: self
  type( jedi_increment_type ),          intent(inout) :: increment

  ! Create a model_data to propagate
  call increment%create_model_data()

end subroutine model_initTL

!> @brief    Step the TL model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepTL(self, increment)

  use trajectory_copy_fields_mod, only : trajectory_copy_fields

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: linear_state
  type( field_collection_type ), pointer :: depository

  nullify(depository)
  depository => increment%model_data%get_field_collection("depository")

  ! Get the linear state to be used with the linear model time-step
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      linear_state )

  !> @todo The linear model step would be called here but it is not yet
  !> available. The fields are stored in the model_data and includes the
  !> linear-state fields inside model_data. For each field, there is an
  !> increment to be propagated and a linear state that has the same name with
  !> “_ls” postpended. In the LFRic linear model, the linear-state is updated
  !>  via a file read but here the data will be copied from the linear-state
  !> trajectory.

  !> When the LFRic linear model is available, we will need to:
  !> 1. Copy the linear-state values from the linear-state into the model_data
  !> “_ls” fields. A copy method such as the one used here could do that, e.g.:
  !>  call trajectory_copy_fields(increment%model_data%depository, linear_state)
  !> 2. Call the linear step method with the model_data:
  !>  call linear_step( mesh,       &
  !>                    twod_mesh,  &
  !>                    model_data, &
  !>                    model_clock )
  !> 3. copy model data into Atlas fields

  ! For now, we can partially test this by simply copying the linear state into
  ! the model_data
  call trajectory_copy_fields( depository, linear_state )

  ! Copy model_data back to the Atlas fields
  call increment%from_model_data()

  ! Update the increment time
  call increment%update_time( self%time_step )

end subroutine model_stepTL

!> @brief    Finalise the TL model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalTL(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(in) :: self
  type( jedi_increment_type ),          intent(inout) :: increment

  !> @todo we will have a destroy call here to compliment the create in init

end subroutine model_finalTL

!> @brief    Initialise the AD model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initAD(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(in) :: self
  type( jedi_increment_type ),          intent(inout) :: increment

  !> @todo add the adjoint init method

end subroutine model_initAD

!> @brief    Step the AD model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepAD(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),             intent(inout) :: increment

  !> @todo add the adjoint step method

end subroutine model_stepAD

!> @brief    Finalise the TL model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalAD(self, increment)

  implicit none

  class( jedi_linear_model_type ), target, intent(in) :: self
  type( jedi_increment_type ),          intent(inout) :: increment

  !> @todo add the adjoint final method

end subroutine model_finalAD

!> @brief    Finalize the jedi_pseudo_model_type
!>
subroutine jedi_linear_model_destructor(self)

  implicit none

  type(jedi_linear_model_type), intent(inout) :: self

end subroutine jedi_linear_model_destructor

!------------------------------------------------------------------------------
! OOPS defined forecastTL and forecastAD methods
!------------------------------------------------------------------------------

!> @brief    Run a forecastTL using the model init, step and final
!>
!> @param [inout] increment       The Increment object to propagate
!> @param [in]    forecast_length The duration of the forecastTL
subroutine forecastTL( self, increment, forecast_length )

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment
  type( jedi_duration_type ),         intent(in) :: forecast_length

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! End time
  end_time = increment%valid_time() + forecast_length

  ! Initialize the model
  call self%model_initTL( increment )
  ! Initialize the post processor and call first process
  ! call post_processor%pp_init( increment )
  ! call post_processor%process( increment )

  ! Loop until date_time_end
  do while ( end_time%is_ahead( increment%valid_time() ) )
    call self%model_stepTL( increment )
    ! call post_processor%process( increment )
  end do

  ! Finalize model and post processor
  ! call post_processor%pp_init( increment )
  call self%model_finalTL( increment )

end subroutine forecastTL

!> @brief    Run a forecastAD using the model init, step and final
!>
!> @param [inout] increment       The Increment object to propagate
!> @param [in]    forecast_length The duration of the forecastAD
subroutine forecastAD(self, increment, forecast_length)

  implicit none

  class( jedi_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),     intent(inout) :: increment
  type( jedi_duration_type ),         intent(in) :: forecast_length

  ! Local
  type( jedi_datetime_type ) :: end_time

  ! End time
  end_time = increment%valid_time() + forecast_length

  ! Initialize the model
  call self%model_initAD( increment )
  ! Initialize the post processor and call first process
  ! call post_processor%pp_init( increment )
  ! call post_processor%process( increment )

  ! Loop until date_time_end
  do while ( end_time%is_ahead( increment%valid_time() ) )
    call self%model_stepAD( increment )
    ! call post_processor%process( increment )
  end do

  ! Finalize model and post processor
  ! call post_processor%pp_init( increment )
  call self%model_finalAD( increment )

end subroutine forecastAD

end module jedi_linear_model_mod
