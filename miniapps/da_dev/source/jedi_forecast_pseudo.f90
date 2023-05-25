!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_forecast_pseudo program

!> @brief Main program for running pseudo model forecast with jedi emulator
!>        objects.

!> @details Setup and run a pseudo model forecast using the jedi emulator
!>          objects. The jedi objects are constructed via an initialiser call
!>          and the forecast is handled by the model object.
!>
program jedi_forecast_pseudo

  use constants_mod,     only : i_def, i_native
  use da_dev_driver_mod, only : finalise_model

  ! Data types and methods to get/store configurations
  use jedi_state_config_mod,        only : jedi_state_config_type
  use jedi_pseudo_model_config_mod, only : jedi_pseudo_model_config_type
  use cli_mod,                      only : get_initial_filename

  ! Jedi emulator objects
  use jedi_run_mod,          only : jedi_run_type
  use jedi_geometry_mod,     only : jedi_geometry_type
  use jedi_state_mod,        only : jedi_state_type
  use jedi_pseudo_model_mod, only : jedi_pseudo_model_type


  implicit none

  ! Jedi objects
  type(jedi_geometry_type)     :: jedi_geometry
  type(jedi_state_type)        :: jedi_state
  type(jedi_pseudo_model_type) :: jedi_model
  type(jedi_run_type)          :: jedi_run

  ! Emulator configs
  type(jedi_state_config_type)        :: jedi_state_config
  type(jedi_pseudo_model_config_type) :: jedi_pseudo_model_config
  integer( kind=i_def )               :: date_time_duration

  ! Local
  character(:), allocatable :: filename
  integer( kind=i_native )  :: model_communicator
  character(*), parameter   :: program_name = "jedi_forecast_pseudo"

  ! Infrastructure config
  call get_initial_filename( filename )

  ! Run object - handles initialization and finalization of required infrastructure
  ! Initialize external libraries such as XIOS
  call jedi_run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call jedi_run%initialise_infrastructure( filename, model_communicator )

  ! Config for the jedi emulator objects
  ! State config
  call jedi_state_config%initialise( use_nl_model = .false. )

  ! Model config
  call jedi_pseudo_model_config%initialise()

  ! Forecast config - duration of forecast / seconds
  date_time_duration = 5_i_def

  ! Geometry
  call jedi_geometry%initialise()

  ! State
  call jedi_state%initialise( jedi_geometry, jedi_state_config )

  ! Model
  call jedi_model%initialise( jedi_pseudo_model_config )

  ! Run app via model class
  call jedi_model%forecast( jedi_state, date_time_duration )

  ! To provide KGO
  call finalise_model( program_name, jedi_state%io_collection )

end program jedi_forecast_pseudo
