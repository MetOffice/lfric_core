!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_tlm_forecast_tl program

!> @brief Main program for running linear forecast with jedi emulator
!>        objects.
!>
!> @details Setup and run a linear model forecastTL using the JEDI
!>          emulator objects. The linear state trajectory is provided via the
!>          pseudo model forecast. The jedi objects are constructed via an
!>          initialiser call and the forecasts are handled by the model
!>          objects.
!>
program jedi_tlm_forecast_tl

  use constants_mod,           only : PRECISION_REAL, i_def, i_native
  use log_mod,                 only : log_event, log_scratch_space, &
                                      LOG_LEVEL_ALWAYS
  use field_collection_mod,    only : field_collection_type

  ! Data types and methods to get/store configurations
  use jedi_state_config_mod,        only : jedi_state_config_type
  use jedi_pseudo_model_config_mod, only : jedi_pseudo_model_config_type
  use jedi_increment_config_mod,    only : jedi_increment_config_type
  use jedi_linear_model_config_mod, only : jedi_linear_model_config_type
  use cli_mod,                      only : get_initial_filename

  ! Jedi emulator objects
  use jedi_checksum_mod,            only : output_checksum
  use jedi_lfric_duration_mod,      only : jedi_duration_type
  use jedi_run_mod,                 only : jedi_run_type
  use jedi_geometry_mod,            only : jedi_geometry_type
  use jedi_state_mod,               only : jedi_state_type
  use jedi_increment_mod,           only : jedi_increment_type
  use jedi_pseudo_model_mod,        only : jedi_pseudo_model_type
  use jedi_linear_model_mod,        only : jedi_linear_model_type
  use jedi_post_processor_traj_mod, only : jedi_post_processor_traj_type

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )            :: jedi_geometry
  type( jedi_state_type )               :: jedi_state
  type( jedi_increment_type )           :: jedi_increment
  type( jedi_pseudo_model_type )        :: jedi_psuedo_model
  type( jedi_linear_model_type )        :: jedi_linear_model
  type( jedi_run_type )                 :: jedi_run
  type( jedi_post_processor_traj_type ) :: pp_traj

  ! Emulator object configs
  type( jedi_state_config_type )        :: jedi_state_config
  type( jedi_increment_config_type )    :: jedi_increment_config
  type( jedi_pseudo_model_config_type ) :: jedi_pseudo_model_config
  type( jedi_linear_model_config_type ) :: jedi_linear_model_config
  type(jedi_duration_type)              :: forecast_length

  ! Local
  character(:),              allocatable :: filename
  integer( kind=i_native )               :: model_communicator
  character(*),                parameter :: program_name = "jedi_tlm_forecast_tl"
  type( field_collection_type ), pointer :: depository => null()

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

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
  call jedi_state_config%initialise( use_pseudo_model = .true. )

  ! Increment config
  call jedi_increment_config%initialise()

  ! Linear Model config
  call jedi_linear_model_config%initialise()

  ! Model config
  call jedi_pseudo_model_config%initialise()

  ! Forecast length
  call forecast_length%init( 'P0DT6H0M0S' )

  ! Create geometry
  call jedi_geometry%initialise()

  ! Create state
  call jedi_state%initialise( program_name, jedi_geometry, jedi_state_config )

  ! Create increment
  call jedi_increment%initialise( jedi_geometry, jedi_increment_config )

  ! Create linear model
  call jedi_linear_model%initialise( jedi_linear_model_config )

  ! Initialise trajectory post processor with instance of jedi_linear_model
  call pp_traj%initialise( jedi_linear_model )

  ! Create non-linear model
  call jedi_psuedo_model%initialise( jedi_pseudo_model_config )

  ! Run non-linear model forecast to populate the trajectory object
  call jedi_psuedo_model%forecast( jedi_state, forecast_length, pp_traj )

  ! Run the linear model TL forecast
  call jedi_linear_model%forecastTL( jedi_increment, forecast_length )

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  ! To provide KGO
  depository => jedi_increment%model_data%get_field_collection("depository")
  call output_checksum( program_name, depository )

end program jedi_tlm_forecast_tl
