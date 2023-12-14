!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use constants_mod,              only : r_def
  use derived_config_mod,         only : l_esm_couple
  use driver_io_mod,              only : get_io_context
  use extrusion_mod,              only : TWOD
  use field_collection_mod,       only : field_collection_type
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_init_fields_mod,     only : create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_modeldb_mod,         only : modeldb_type
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_step_mod,            only : gungho_step
  use initial_output_mod,         only : write_initial_output
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use initialization_config_mod,  only : lbc_option,               &
                                         lbc_option_gungho_file,   &
                                         lbc_option_um2lfric_file, &
                                         ancil_option,             &
                                         ancil_option_updating,    &
                                         coarse_aerosol_ancil
  use init_gungho_lbcs_alg_mod,   only : update_lbcs_file_alg
  use log_mod,                    only : log_event,           &
                                         log_level_always,    &
                                         log_level_error,     &
                                         log_level_info,      &
                                         log_scratch_space
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use mpi_mod,                    only : mpi_type
  use multires_coupling_config_mod, &
                                  only : aerosol_mesh_name
#ifdef UM_PHYSICS
  use variable_fields_mod,        only : update_variable_fields
  use lfric_xios_time_axis_mod,   only : regridder
  use intermesh_mappings_alg_mod, only : map_scalar_field
  use update_ancils_alg_mod,      only : update_ancils_alg
  use gas_calc_all_mod,           only : gas_calc_all
#endif
#ifdef COUPLED
  use esm_couple_config_mod,      only : l_esm_couple_test
  use coupler_mod,                only : cpl_snd, cpl_rcv, cpl_fld_update
  use coupler_diagnostics_mod,    only : save_sea_ice_frac_previous
#endif

  implicit none

  private
  public initialise, step, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine initialise( program_name, modeldb, calendar )

    use io_context_mod,         only : io_context_type
    use lfric_xios_context_mod, only : lfric_xios_context_type

    implicit none

    character(*),         intent(in)    :: program_name
    type(modeldb_type),   intent(inout) :: modeldb
    class(calendar_type), intent(in)    :: calendar

    class(io_context_type), pointer :: io_context => null()
    type(mesh_type),        pointer :: mesh              => null()
    type(mesh_type),        pointer :: twod_mesh         => null()
    type(mesh_type),        pointer :: aerosol_mesh      => null()
    type(mesh_type),        pointer :: aerosol_twod_mesh => null()

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( modeldb, &
                                    calendar )

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

    ! Rate of temperature adjustment for energy correction
    call modeldb%values%add_key_value( 'temperature_correction_rate', &
                                       0.0_r_def )
    ! Total mass of dry atmosphere used for energy correction
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    ! Total energy of moist atmosphere for calculating energy correction
    call modeldb%values%add_key_value( 'total_energy', 0.0_r_def )
    ! Total energy of moist atmosphere at previous energy correction step
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,         &
                            mesh, twod_mesh, &
                            aerosol_mesh, aerosol_twod_mesh )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( modeldb, mesh, twod_mesh )

    ! Initial output
    io_context => get_io_context()
    call write_initial_output( modeldb, mesh, twod_mesh, &
                               io_context, nodal_output_on_w3 )

    ! Model configuration initialisation
    call initialise_model( mesh, modeldb )

#ifdef COUPLED
    ! Placeholder for ESM coupling initialisation code.
    ! Check we have a value for related namelist control variable
    write(log_scratch_space,'(A,L1)') program_name//': Couple flag l_esm_couple_test: ', &
                                     l_esm_couple_test
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
#endif

    nullify(mesh, twod_mesh, aerosol_mesh, aerosol_twod_mesh)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Timestep the model, calling the desired timestepping algorithm
  !>        based upon the configuration.
  !> @param [in,out] modeldb   The structure that holds model state
  !>
  subroutine step( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(mesh_type), pointer :: mesh      => null()
    type(mesh_type), pointer :: twod_mesh => null()

#ifdef COUPLED
    type( field_collection_type ), pointer :: depository => null()
#endif

#ifdef UM_PHYSICS
    procedure(regridder), pointer :: regrid_operation => null()

    regrid_operation => map_scalar_field
#endif
    ! Get primary and 2D meshes
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

#ifdef COUPLED
    if(l_esm_couple) then

       write(log_scratch_space, &
             '(A, I0)') 'Coupling timestep: ', modeldb%clock%get_step() - 1
       call log_event( log_scratch_space, LOG_LEVEL_INFO )

       depository => modeldb%fields%get_field_collection("depository")

       call save_sea_ice_frac_previous(depository)

       ! Receive all incoming (ocean/seaice fields) from the coupler
       call cpl_rcv( modeldb%model_data%cpl_rcv,    &
                     depository, &
                     modeldb%clock )

       ! Send all outgoing (ocean/seaice driving fields) to the coupler
       call cpl_snd( modeldb%model_data%cpl_snd,    &
                     depository, &
                     modeldb%clock )

    endif
#endif

    if ( lbc_option == lbc_option_gungho_file .or. &
         lbc_option == lbc_option_um2lfric_file) then

      call update_lbcs_file_alg( modeldb%model_axes%lbc_times_list, &
                                 modeldb%clock, modeldb%model_data%lbc_fields )
    endif

    ! Perform a timestep
    call gungho_step( mesh, twod_mesh, modeldb, modeldb%clock )

    ! Use diagnostic output frequency to determine whether to write
    ! diagnostics on this timestep

    if ( ( mod(modeldb%clock%get_step(), diagnostic_frequency) == 0 ) &
         .and. ( write_diag ) ) then

      ! Calculation and output diagnostics
      call gungho_diagnostics_driver( modeldb, mesh, twod_mesh, &
                                      nodal_output_on_w3 )
    end if

#ifdef UM_PHYSICS
    ! Update time-varying ancillaries
    ! This is done last in the timestep, because the time data of the
    ! ancillaries needs to be that of the start of timestep, but the
    ! first thing that happens in a timestep is that the clock ticks to the
    ! end of timestep date.
    if ( ancil_option == ancil_option_updating ) then
      call update_variable_fields( modeldb%model_axes%ancil_times_list, &
                                   modeldb%clock,                       &
                                   modeldb%model_data%ancil_fields,     &
                                   regrid_operation )
      call update_ancils_alg( modeldb%model_axes%ancil_times_list, &
                              modeldb%clock,                       &
                              modeldb%model_data%ancil_fields,     &
                              modeldb%model_data%surface_fields)
    end if

    ! Update the time varying trace gases
    call gas_calc_all()
#endif

    nullify(mesh, twod_mesh)

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !> @param [in]     program_name An identifier given to the model begin run
  !> @param [in,out] modeldb      The structure that holds model state
  !>
  subroutine finalise( program_name, modeldb )

    implicit none

    character(*), intent(in)          :: program_name
    type(modeldb_type), intent(inout) :: modeldb

#ifdef COUPLED
    type( field_collection_type ), pointer :: depository => null()

    if (l_esm_couple) then
       depository => modeldb%fields%get_field_collection("depository")
       ! Ensure coupling fields are updated at the end of a cycle to ensure the values
       ! stored in and recovered from checkpoint dumps are correct and reproducible
       ! when (re)starting subsequent runs!
       call cpl_fld_update(modeldb%model_data%cpl_snd,    &
                           depository,                    &
                           modeldb%clock)
    endif
#endif

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( modeldb%model_data, modeldb%clock )

    ! Model configuration finalisation
    call finalise_model( modeldb,               &
                         modeldb%configuration, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure()

  end subroutine finalise

end module gungho_driver_mod
