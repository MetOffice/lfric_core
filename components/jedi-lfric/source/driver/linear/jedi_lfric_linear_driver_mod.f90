!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to drive the execution of
!>           the (tangent) linear model from LFRic-JEDI
!>
module jedi_lfric_linear_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use constants_mod,              only : i_def, imdi, r_def
  use driver_io_mod,              only : get_io_context
  use extrusion_mod,              only : TWOD
  use field_mod,                  only : field_type
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use jedi_lfric_tlm_model_mod,   only : initialise_infrastructure, &
                                         initialise_model,          &
                                         finalise_infrastructure,   &
                                         finalise_model
  use gungho_init_fields_mod,     only : create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_modeldb_mod,         only : modeldb_type
  use gungho_step_mod,            only : gungho_step
  use init_fd_prognostics_mod,    only : init_fd_prognostics_dump
  use initial_output_mod,         only : write_initial_output
  use initialization_config_mod,  only : ls_option,                 &
                                         ls_option_file,            &
                                         init_option,               &
                                         init_option_fd_start_dump, &
                                         coarse_aerosol_ancil
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS,  &
                                         LOG_LEVEL_INFO
  use linear_model_data_mod,      only : linear_create_ls, &
                                         linear_init_ls,   &
                                         linear_init_pert
  use linear_model_mod,           only : initialise_linear_model, &
                                         finalise_linear_model
  use linear_diagnostics_driver_mod, &
                                  only : linear_diagnostics_driver
  use linear_step_mod,            only : linear_step
  use linear_data_algorithm_mod,  only : update_ls_file_alg
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use model_clock_mod,            only : model_clock_type
  use multires_coupling_config_mod, &
                                  only : aerosol_mesh_name
  use create_tl_prognostics_mod,  only : create_tl_prognostics

  implicit none

  private
  public initialise, step, finalise

  type( mesh_type ), pointer :: mesh                 => null()
  type( mesh_type ), pointer :: twod_mesh            => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine initialise( program_name, modeldb )

    implicit none

    character(*),           intent(in)    :: program_name
    type(modeldb_type),     intent(inout) :: modeldb

    class(io_context_type), pointer :: io_context        => null()
    type( mesh_type ),      pointer :: aerosol_mesh      => null()
    type( mesh_type ),      pointer :: aerosol_twod_mesh => null()

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( program_name,       &
                                    modeldb%model_data, &
                                    modeldb%clock,      &
                                    modeldb%mpi )

    ! Get primary and 2D meshes for initialising model data
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! If aerosol data is on a different mesh, get this
    if (coarse_aerosol_ancil) then
      ! For now use the coarsest mesh
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      aerosol_twod_mesh => mesh_collection%get_mesh(aerosol_mesh, TWOD)
      write( log_scratch_space,'(A,A)' ) "aerosol mesh name:", aerosol_mesh%get_mesh_name()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    else
      aerosol_mesh => mesh
      aerosol_twod_mesh => twod_mesh
    end if

    call modeldb%values%add_key_value( 'temperature_correction_rate', &
                                       0.0_r_def )
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,         &
                            mesh, twod_mesh, &
                            aerosol_mesh, aerosol_twod_mesh )

    ! Instantiate the fields required to read the initial
    ! conditions from a file.
    if ( init_option == init_option_fd_start_dump ) then
      call create_tl_prognostics( mesh, twod_mesh,               &
                                  modeldb%model_data%fd_fields,  &
                                  modeldb%model_data%depository)
    end if

    ! Instantiate the linearisation state
    call linear_create_ls( modeldb, mesh )

    ! Initialise the fields stored in the model_data
    if ( init_option == init_option_fd_start_dump ) then
      call init_fd_prognostics_dump( modeldb%model_data%fd_fields )
    else
      call initialise_model_data( modeldb, mesh, twod_mesh )
    end if

    ! Model configuration initialisation
    call initialise_model( mesh,  &
                           modeldb%model_data )

    ! Initialise the linearisation state
    call linear_init_ls( mesh, twod_mesh, modeldb )

    ! Initialise the linear model perturbation state
    call linear_init_pert( mesh,      &
                           twod_mesh, &
                           modeldb%model_data )

    ! Initial output
    io_context => get_io_context()
    call write_initial_output( modeldb, mesh, twod_mesh, &
                               io_context, nodal_output_on_w3 )

    ! Linear model configuration initialisation
    call initialise_linear_model( mesh,        &
                                  modeldb%model_data )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Timesteps the model, calling the desired timestepping algorithm
  !         based upon the configuration
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine step( modeldb, model_clock )

    implicit none

    type(modeldb_type),     intent(inout) :: modeldb
    type(model_clock_type), intent(inout) :: model_clock

    if ( ls_option == ls_option_file ) then
      call update_ls_file_alg( modeldb%model_axes%ls_times_list, &
                               model_clock,                      &
                               modeldb%model_data%ls_fields,     &
                               modeldb%model_data%ls_mr,         &
                               modeldb%model_data%ls_moist_dyn )
    end if

    call linear_step( mesh,       &
                      twod_mesh,  &
                      modeldb   , &
                      model_clock )

    if ( ( mod(model_clock%get_step(), diagnostic_frequency) == 0 ) &
          .and. ( write_diag ) ) then

      ! Calculation and output diagnostics
      call gungho_diagnostics_driver( modeldb, mesh, twod_mesh, &
                                      nodal_output_on_w3 )

      call linear_diagnostics_driver( mesh,                &
                                      modeldb%model_data,  &
                                      model_clock,         &
                                      nodal_output_on_w3 )
    end if

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine finalise( program_name, modeldb, model_clock )

    implicit none

    character(*),           intent(in)    :: program_name
    type(modeldb_type),     intent(inout) :: modeldb
    type(model_clock_type), intent(inout) :: model_clock

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! ! Output the fields stored in the model_data (checkpoint and dump)
    ! call output_model_data( model_data, model_clock )

    ! Model configuration finalisation
    call finalise_model( modeldb%model_data, &
                         program_name )

    call finalise_linear_model()

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb%model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module jedi_lfric_linear_driver_mod
