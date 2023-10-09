!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a configuration for the JEDI linear model
!>        emulator.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct a JEDI linear model emulator. An initialiser is included
!>          and this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_linear_model_config_mod

  use constants_mod,             only : i_def, str_def, l_def
  use jedi_lfric_field_meta_mod, only : jedi_lfric_field_meta_type
  use fs_continuity_mod,         only : W3, Wtheta
  use jedi_lfric_datetime_mod,   only : jedi_datetime_type
  use jedi_lfric_duration_mod,   only : jedi_duration_type

  implicit none

  private

type, public :: jedi_linear_model_config_type

  !> The field meta data
  type( jedi_lfric_field_meta_type ) :: field_meta_data
  !> The TLM forecast_length
  type( jedi_duration_type )         :: forecast_length
  !> The TLM time-step
  type( jedi_duration_type )         :: time_step

contains

  !> Initialiser.
  procedure :: initialise

  !> jedi_state_config finalizer
  final     :: jedi_state_config_destructor

end type jedi_linear_model_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_linear_model_config_type
!>
subroutine initialise( self )

  implicit none

  class( jedi_linear_model_config_type ), intent(inout) :: self

  ! Local
  integer( kind=i_def ), parameter :: nvars = 2
  character( len=str_def )         :: variable_names(nvars)
  integer( kind=i_def )            :: variable_function_spaces(nvars)
  logical( kind=l_def )            :: variable_is_2d(nvars)

  ! Configuration inputs

  ! Init model time-step and forecast length
  call self%time_step%init( 'P0DT1H0M0S' )
  call self%forecast_length%init( 'P0DT6H0M0S' )

  ! Setup arrays required for field_meta_data
  ! Variable names
  variable_names(1) = "theta"
  variable_names(2) = "rho"
  ! Variable function spaces
  variable_function_spaces(1) = Wtheta
  variable_function_spaces(2) = W3
  ! Variable is_2d
  variable_is_2d(1) = .false.
  variable_is_2d(2) = .false.

  call self%field_meta_data%initialise( variable_names,           &
                                        variable_function_spaces, &
                                        variable_is_2d )

end subroutine initialise

!> @brief    Finalizer for jedi_linear_model_config_type
!>
subroutine jedi_state_config_destructor( self )

  implicit none

  type( jedi_linear_model_config_type ), intent(inout)    :: self

end subroutine jedi_state_config_destructor

end module jedi_linear_model_config_mod
