!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the diagnostics miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module diagnostics_driver_mod

  use clock_mod,                     only : clock_type
  use constants_mod,                 only : i_def, i_native, str_def, r_def
  use diagnostics_configuration_mod, only : load_configuration, program_name
  use field_mod,                     only : field_type
  use field_parent_mod,              only : field_parent_type
  use field_collection_mod,          only : field_collection_type
  use fieldspec_collection_mod,      only : fieldspec_collection
  use diagnostics_model_data_mod,    only : model_data_type
  use io_config_mod,                 only : write_diag, &
                                            use_xios_io
  use io_context_mod,                only : io_context_type
  use local_mesh_collection_mod,     only : local_mesh_collection, &
                                            local_mesh_collection_type
  use log_mod,                       only : log_event, &
                                            log_set_level, &
                                            log_scratch_space, &
                                            initialise_logging, &
                                            finalise_logging, &
                                            LOG_LEVEL_ALWAYS, &
                                            LOG_LEVEL_ERROR, &
                                            LOG_LEVEL_WARNING, &
                                            LOG_LEVEL_INFO, &
                                            LOG_LEVEL_DEBUG, &
                                            LOG_LEVEL_TRACE
  use mesh_collection_mod,           only : mesh_collection, &
                                            mesh_collection_type
  use mesh_mod,                      only : mesh_type
  use mpi_mod,                       only : store_comm,    &
                                            get_comm_size, &
                                            get_comm_rank
  use simple_io_mod,                 only : initialise_simple_io
  use yaxt,                          only : xt_initialize, xt_finalize

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type), target :: model_data

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id

  type(mesh_type), pointer :: mesh      => null()
  type(mesh_type), pointer :: twod_mesh => null()


  class(io_context_type), allocatable :: io_context

  character(len = *), public, parameter :: xios_ctx = program_name
  character(len = *), public, parameter :: xios_id = "lfric_client"

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !>
  !> mostly boiler plate - note the init and seeding of the fields at the end
  !> of the function.
  !>
  subroutine initialise( filename, model_communicator )

    use convert_to_upper_mod,       only : convert_to_upper
    use driver_fem_mod,             only : init_fem
    use driver_mesh_mod,            only : init_mesh
    use fieldspec_xml_parser_mod,   only : populate_fieldspec_collection
    use init_diagnostics_mod,       only : init_diagnostics
    use lfric_xios_io_mod,          only : initialise_xios
    use logging_config_mod,         only : run_log_level, &
                                           key_from_run_log_level, &
                                           RUN_LOG_LEVEL_ERROR, &
                                           RUN_LOG_LEVEL_INFO, &
                                           RUN_LOG_LEVEL_DEBUG, &
                                           RUN_LOG_LEVEL_TRACE, &
                                           RUN_LOG_LEVEL_WARNING
    use mod_wait,                   only : init_wait
    use seed_diagnostics_mod,       only : seed_diagnostics
    use time_config_mod,            only : timestep_end,   &
                                           timestep_start, &
                                           calendar_start, &
                                           calendar_type,  &
                                           key_from_calendar_type
    use timestepping_config_mod,    only : dt, &
                                           spinup_period
    use diagnostics_miniapp_config_mod, only : iodef_path

    implicit none

    character(:),      intent(in), allocatable :: filename
    integer(i_native), intent(in) :: model_communicator

    character(len = *), parameter :: xios_ctx = "diagnostics"

    class(clock_type), pointer :: clock
    real(r_def)                :: dt_model

    integer(i_def)     :: total_ranks, local_rank, stencil_depth

    integer(i_native) :: log_level

    ! Store the MPI communicator for later use
    call store_comm( model_communicator )

    ! Initialise YAXT
    call xt_initialize( model_communicator )

    ! and get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank = get_comm_rank()

    call initialise_logging(local_rank, total_ranks, program_name)

    call load_configuration(filename)

    select case (run_log_level)
    case(RUN_LOG_LEVEL_ERROR)
        log_level = LOG_LEVEL_ERROR
    case(RUN_LOG_LEVEL_WARNING)
        log_level = LOG_LEVEL_WARNING
    case(RUN_LOG_LEVEL_INFO)
        log_level = LOG_LEVEL_INFO
    case(RUN_LOG_LEVEL_DEBUG)
        log_level = LOG_LEVEL_DEBUG
    case(RUN_LOG_LEVEL_TRACE)
        log_level = LOG_LEVEL_TRACE
    end select

    call log_set_level(log_level)

    write(log_scratch_space, '(A)')                            &
            'Runtime message logging severity set to log level: ' // &
                    convert_to_upper(key_from_run_log_level(run_log_level))
    call log_event(log_scratch_space, LOG_LEVEL_ALWAYS)

    !----------------------------------------------------------------------
    ! Model init
    !----------------------------------------------------------------------
    call log_event( 'Initialising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    allocate( local_mesh_collection, &
            source = local_mesh_collection_type() )
    allocate( mesh_collection, &
              source=mesh_collection_type() )

    stencil_depth = 1_i_def

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, stencil_depth, &
                    mesh, twod_mesh=twod_mesh )


    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    !----------------------------------------------------------------------
    ! IO init
    !----------------------------------------------------------------------

    call log_event("Populating fieldspec collection", LOG_LEVEL_INFO)
    call populate_fieldspec_collection(iodef_path)

    if (use_xios_io) then
      call log_event( "init XIOS", log_level_info )
      call initialise_xios( io_context,         &
                            xios_ctx,           &
                            model_communicator, &
                            mesh,               &
                            twod_mesh,          &
                            chi,                &
                            panel_id,           &
                            timestep_start,     &
                            timestep_end,       &
                            spinup_period,      &
                            dt,                 &
                            calendar_start,     &
                            key_from_calendar_type(calendar_type) )
    else
      call initialise_simple_io( io_context,     &
                                 timestep_start, &
                                 timestep_end,   &
                                 spinup_period,  &
                                 dt )
    end if

    clock => io_context%get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create and initialise prognostic fields
    call init_diagnostics( mesh, twod_mesh,                &
                           chi, panel_id, dt_model,        &
                           model_data, fieldspec_collection)

    call log_event("seed starting values", LOG_LEVEL_INFO)
    ! Seed values as this is a test!
    call seed_diagnostics(model_data)

    call log_event("finish init", LOG_LEVEL_INFO)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    use diagnostics_step_mod,       only : diagnostics_step

    implicit none

    class(clock_type), pointer :: clock

    clock => io_context%get_clock()

    ! standard timestepping from gungho
    do while (clock%tick())

        write(log_scratch_space, '("/", A, "\ ")') repeat("*", 76)
        call log_event(log_scratch_space, LOG_LEVEL_TRACE)
        write( log_scratch_space, &
               '(A,I0)' ) 'Start of timestep ', clock%get_step()
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        call log_event( 'Running ' // program_name // ' ...', &
                        LOG_LEVEL_ALWAYS )
        call diagnostics_step( mesh,       &
                               twod_mesh,  &
                               model_data, &
                               clock )

    end do

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    use checksum_alg_mod,  only : checksum_alg
    use configuration_mod, only : final_configuration
    use fieldspec_mod,     only : fieldspec_type

    implicit none

    type(field_collection_type), pointer :: depository
    class(field_type), pointer :: hex
    class(field_type), pointer :: mutable_numbers
    class(field_type), pointer :: mutable_categories
    class(field_type), pointer :: immutable_both

    !----------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------
    call log_event( 'Finalising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    depository => model_data%depository
    ! as with the run step this could use a specific checksum collection to
    ! control if it outputs a checksum for a given field
    hex => depository%get_field("colours__hex")
    mutable_numbers => depository%get_field("colours__mutable_numbers")
    mutable_categories => depository%get_field("colours__mutable_categories")
    immutable_both => depository%get_field("colours__immutable_both")
    call checksum_alg('diagnostics', &
            hex, trim(hex%get_name()), &
            mutable_numbers, trim(mutable_numbers%get_name()), &
            mutable_categories, trim(mutable_categories%get_name()), &
            immutable_both, trim(immutable_both%get_name()))

    !----------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------

    deallocate( io_context )

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    call log_event(program_name // ' completed.', LOG_LEVEL_ALWAYS)

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise

end module diagnostics_driver_mod
