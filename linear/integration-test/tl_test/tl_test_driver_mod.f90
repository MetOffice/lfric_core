!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   Drives the execution of the tangent linear model tests.
!>@details The tests are initialised and finalised using a similar
!!         method to gungho, but with the addition of the linearisation state.
module tl_test_driver_mod

  use constants_mod,              only : i_def, i_native, imdi
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model,          &
                                         finalise_infrastructure,   &
                                         finalise_model
  use gungho_model_data_mod,      only : model_data_type,       &
                                         create_model_data,     &
                                         initialise_model_data, &
                                         finalise_model_data
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         LOG_LEVEL_ALWAYS
  use mesh_mod,                   only : mesh_type
  use model_clock_mod,            only : model_clock_type
  use linear_model_data_mod,      only : linear_create_ls,  &
                                         linear_init_ls
  use tl_test_kinetic_energy_gradient_mod, only : test_kinetic_energy_gradient
  use tl_test_advect_density_field_mod,    only : test_advect_density_field
  use tl_test_advect_theta_field_mod,      only : test_advect_theta_field
  use tl_test_vorticity_mod,               only : test_vorticity_advection
  use tl_test_project_eos_pressure_mod,    only : test_project_eos_pressure
  use tl_test_sample_eos_pressure_mod,     only : test_sample_eos_pressure
  use tl_test_hydrostatic_mod,             only : test_hydrostatic
  use tl_test_pressure_grad_bd_mod,        only : test_pressure_gradient_bd
  use tl_test_rk_alg_mod,                  only : test_rk_alg
  use tl_test_transport_control_mod,       only : test_transport_control
  use tl_test_rhs_sample_eos_mod,          only : test_rhs_sample_eos
  use tl_test_rhs_project_eos_mod,         only : test_rhs_project_eos
  use tl_test_rhs_alg_mod,                 only : test_rhs_alg
  use tl_test_semi_imp_alg_mod,            only : test_semi_imp_alg
  use tl_test_timesteps_alg_mod,           only : test_timesteps

  implicit none

  private
  public initialise,                  &
         finalise,                    &
         run_timesteps,               &
         run_kinetic_energy_gradient, &
         run_advect_density_field,    &
         run_advect_theta_field,      &
         run_vorticity_advection,     &
         run_project_eos_pressure,    &
         run_sample_eos_pressure,     &
         run_hydrostatic,             &
         run_pressure_gradient_bd,    &
         run_rk_alg,                  &
         run_rhs_alg,                 &
         run_rhs_project_eos,         &
         run_rhs_sample_eos,          &
         run_semi_imp_alg,            &
         run_transport_control

  type(model_data_type)               :: model_data
  type(model_clock_type), allocatable :: model_clock

  type(mesh_type), pointer :: mesh              => null()
  type(mesh_type), pointer :: twod_mesh         => null()
  type(mesh_type), pointer :: shifted_mesh      => null()
  type(mesh_type), pointer :: double_level_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief     Sets up the required state in preparation for run.
  !>@param[in] filename            Name of the file containing the desired
  !!                               configuration
  subroutine initialise( filename, program_name )

    implicit none

    character(*), intent(in) :: filename
    character(*), intent(in) :: program_name

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( filename,          &
                                    program_name,      &
                                    mesh,              &
                                    twod_mesh,         &
                                    shifted_mesh,      &
                                    double_level_mesh, &
                                    model_data,        &
                                    model_clock )

    ! Instantiate the fields stored in model_data
    call create_model_data( model_data, &
                            mesh,       &
                            twod_mesh,  &
                            model_clock )

    ! Instantiate the linearisation state
    call linear_create_ls( model_data, &
                           mesh,       &
                           twod_mesh )

    ! Initialise the fields stored in the model_data prognostics. This needs
    ! to be done before initialise_model.
    call initialise_model_data( model_data, model_clock, mesh, twod_mesh )

    ! Model configuration initialisation
    call initialise_model( mesh,  &
                           model_data )

    ! Initialise the linearisation state
    call linear_init_ls( mesh,       &
                         twod_mesh,  &
                         model_data, &
                         model_clock )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear kinetic energy gradient kernel.
  subroutine run_timesteps()

    implicit none

    call test_timesteps( model_data, &
                         mesh,       &
                         twod_mesh,  &
                         model_clock )

  end subroutine run_timesteps

  subroutine run_kinetic_energy_gradient()

    implicit none

    call test_kinetic_energy_gradient( model_data, &
                                       mesh,       &
                                       twod_mesh )

  end subroutine run_kinetic_energy_gradient

  subroutine run_advect_density_field()

    implicit none

    call test_advect_density_field( model_data,  &
                                    model_clock, &
                                    mesh,        &
                                    twod_mesh )

  end subroutine run_advect_density_field

  subroutine run_advect_theta_field()

    implicit none

    call test_advect_theta_field( model_data,  &
                                  model_clock, &
                                  mesh,        &
                                  twod_mesh )

  end subroutine run_advect_theta_field

  subroutine run_vorticity_advection()

    implicit none

    call test_vorticity_advection( model_data, &
                                   mesh,       &
                                   twod_mesh )

  end subroutine run_vorticity_advection

  subroutine run_project_eos_pressure()

    implicit none

    call test_project_eos_pressure( model_data, &
                                    mesh,       &
                                    twod_mesh )

  end subroutine run_project_eos_pressure

  subroutine run_sample_eos_pressure()

    implicit none

    call test_sample_eos_pressure( model_data,  &
                                   mesh,        &
                                   twod_mesh )

  end subroutine run_sample_eos_pressure

  subroutine run_hydrostatic()

    implicit none

    call test_hydrostatic( model_data, &
                           mesh,       &
                           twod_mesh )

  end subroutine run_hydrostatic

  subroutine run_pressure_gradient_bd()

    implicit none

    call test_pressure_gradient_bd( model_data, &
                                    mesh,       &
                                    twod_mesh )

  end subroutine run_pressure_gradient_bd

  subroutine run_rk_alg()

    implicit none

    call test_rk_alg( model_data, &
                      mesh,       &
                      twod_mesh,  &
                      model_clock )

  end subroutine run_rk_alg

  subroutine run_transport_control()

    implicit none

    call test_transport_control( model_data,  &
                                 model_clock, &
                                 mesh,        &
                                 twod_mesh )

  end subroutine run_transport_control

  subroutine run_semi_imp_alg()

    implicit none

    call test_semi_imp_alg( model_data, &
                            mesh,       &
                            twod_mesh,  &
                            model_clock )

  end subroutine run_semi_imp_alg

  subroutine run_rhs_alg()

    implicit none

    call test_rhs_alg( model_data, &
                       mesh,       &
                       twod_mesh )

  end subroutine run_rhs_alg

  subroutine run_rhs_project_eos()

    implicit none

    call test_rhs_project_eos( model_data,  &
                               mesh,        &
                               twod_mesh )

  end subroutine run_rhs_project_eos

  subroutine run_rhs_sample_eos()

    implicit none

    call test_rhs_sample_eos( model_data,  &
                              mesh,        &
                              twod_mesh )

  end subroutine run_rhs_sample_eos

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise( program_name )

    implicit none

    character(*), intent(in) :: program_name

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Model configuration finalisation
    call finalise_model( model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module tl_test_driver_mod
