!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> Drives the execution of the gravity_wave miniapp.
!>
module gravity_wave_driver_mod

  use clock_mod,                      only: clock_type
  use constants_mod,                  only: i_def, i_native
  use gravity_wave_infrastructure_mod,only: initialise_infrastructure, &
                                            finalise_infrastructure
  use gravity_wave_grid_mod,          only: initialise_grid
  use gravity_wave_io_mod,            only: initialise_io, &
                                            finalise_io
  use create_gravity_wave_prognostics_mod, &
                                      only: create_gravity_wave_prognostics
  use gravity_wave_constants_config_mod, &
                                      only: b_space,       &
                                            b_space_w0,    &
                                            b_space_w3,    &
                                            b_space_wtheta
  use field_mod,                      only: field_type
  use field_collection_mod,           only: field_collection_type
  use function_space_chain_mod,       only: function_space_chain_type
  use gravity_wave_mod,               only: program_name
  use gravity_wave_diagnostics_driver_mod, &
                                      only: gravity_wave_diagnostics_driver
  use gravity_wave_grid_mod,          only: initialise_grid
  use gravity_wave_infrastructure_mod, &
                                      only: initialise_infrastructure, &
                                            finalise_infrastructure
  use gravity_wave_io_mod,            only: initialise_io, &
                                            finalise_io
  use gw_init_fields_alg_mod,         only: gw_init_fields_alg
  use init_clock_mod,                 only: initialise_clock
  use init_gravity_wave_mod,          only: init_gravity_wave
  use runtime_constants_mod,          only: create_runtime_constants
  use step_gravity_wave_mod,          only: step_gravity_wave
  use final_gravity_wave_mod,         only: final_gravity_wave
  use boundaries_config_mod,          only: limited_area
  use log_mod,                        only: log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ALWAYS,   &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_TRACE,    &
                                            LOG_LEVEL_ERROR
  use read_methods_mod,               only: read_checkpoint
  use write_methods_mod,              only: write_checkpoint
  use io_config_mod,                  only: write_diag,           &
                                            checkpoint_read,      &
                                            checkpoint_write,     &
                                            diagnostic_frequency, &
                                            use_xios_io,          &
                                            nodal_output_on_w3,   &
                                            subroutine_timers
  use files_config_mod,               only: checkpoint_stem_name
  use timer_mod,                      only: init_timer, timer, output_timer
  use xios,                           only: xios_update_calendar

  implicit none

  private

  public initialise, run, finalise

  character(*), public, parameter   :: xios_ctx = 'gravity_wave'

  ! Depository of shared fields
  type( field_collection_type ), target :: depository

  ! Gravity wave model state
  type( field_collection_type ) :: prognostic_fields

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  ! Function space chains
  type(function_space_chain_type) :: multigrid_function_space_chain

  class(clock_type), allocatable :: clock

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

  implicit none

  character(*),      intent(in) :: filename
  integer(i_native), intent(in) :: model_communicator

  ! Initialise aspects of the infrastructure
  call initialise_infrastructure( model_communicator, &
                                  filename,           &
                                  program_name )

  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! The limited area version is unable to work with buoyancy in W0 when using
  ! a biperiodic mesh. However, this could be solved by breaking the continuity
  ! at the edges - so that the W0 dofs on the boundary are not shared, e.g. by
  ! using a non-periodic mesh.
  if ( limited_area ) then
    select case( b_space )
      case( b_space_w0 )
        call log_event( 'Limited area version currently unable to run with buoyancy in W0', LOG_LEVEL_ERROR )
      case( b_space_w3 )
        call log_event( 'Limited area version has not been tested with buoyancy in W3', LOG_LEVEL_ERROR )
      case( b_space_wtheta )
        call log_event( 'Running limited area version with buoyancy in Wtheta', LOG_LEVEL_INFO )
      case default
        call log_event( 'Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select
  endif

  if ( subroutine_timers ) then
    call init_timer()
    call timer(program_name)
  end if

  call initialise_clock( clock )

  multigrid_function_space_chain = function_space_chain_type()

  ! Initialise aspects of the grid
  call initialise_grid(mesh_id, twod_mesh_id, chi, &
                       multigrid_function_space_chain)

  !Initialise aspects of output
  call initialise_io( model_communicator, &
                      clock,              &
                      mesh_id,            &
                      twod_mesh_id,       &
                      chi,                &
                      xios_ctx )

  ! Create runtime_constants object. This in turn creates various things
  ! needed by the timestepping algorithms such as mass matrix operators, mass
  ! matrix diagonal fields and the geopotential field and limited area masks.
  call create_runtime_constants(mesh_id, twod_mesh_id, chi)

  ! Create the depository and prognostics field collections.
  ! Actually, here they will have the same contents (prognostics points to all
  ! the fields in the depository), but both are maintained to be consistent
  ! with more complex models
  depository=field_collection_type(name='depository')
  prognostic_fields=field_collection_type(name='prognostics')

  ! Create the prognostic fields
  call create_gravity_wave_prognostics(mesh_id, depository, prognostic_fields)

  ! Initialise prognostic fields
  if (checkpoint_read) then                 ! Recorded check point to start from
    call read_checkpoint(prognostic_fields, clock%get_step()-1)
  else                                      ! No check point to start from
    call gw_init_fields_alg(prognostic_fields)
  end if

  ! Initialise the gravity-wave model
  call init_gravity_wave( mesh_id, prognostic_fields )

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  if (clock%is_initialisation()) then
    if ( use_xios_io ) then
      ! Need to ensure calendar is initialised here as XIOS has no concept of
      ! timestep 0
      call xios_update_calendar(1)
    end if

    if ( write_diag ) then
      call gravity_wave_diagnostics_driver( mesh_id,           &
                                            prognostic_fields, &
                                            clock,             &
                                            nodal_output_on_w3)
    end if

  end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  write(log_scratch_space,'(A,I0,A)') 'Running '//program_name//' ...'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  !--------------------------------------------------------------------------
  ! Model step
  !--------------------------------------------------------------------------
  do while (clock%tick())

    ! Update XIOS calendar if we are using it for diagnostic output or
    ! checkpoint
    !
    if ( use_xios_io ) then
      call log_event( program_name//': Updating XIOS timestep', &
                      LOG_LEVEL_INFO )
      call xios_update_calendar( clock%get_step() )
    end if

    write( log_scratch_space, '("/",A,"\ ")' ) repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call step_gravity_wave(prognostic_fields)

    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\",A,"/ ")' ) repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if ( (mod(clock%get_step(), diagnostic_frequency) == 0) &
         .and. (write_diag) ) then

      call log_event("Gravity Wave: writing diagnostic output", LOG_LEVEL_INFO)

      call gravity_wave_diagnostics_driver( mesh_id,           &
                                            prognostic_fields, &
                                            clock,             &
                                            nodal_output_on_w3)
    end if

  end do

  end subroutine run


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

  implicit none

  !--------------------------------------------------------------------------
  ! Model finalise
  !--------------------------------------------------------------------------
  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  call final_gravity_wave(prognostic_fields, program_name)

  ! Write checkpoint/restart files if required
  if( checkpoint_write ) then
    call write_checkpoint(prognostic_fields, clock)
  end if

  if ( subroutine_timers ) then
    call timer(program_name)
    call output_timer()
  end if

  !--------------------------------------------------------------------------
  ! Driver layer finalise
  !--------------------------------------------------------------------------

  call finalise_io()

  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  call finalise_infrastructure()

  end subroutine finalise

end module gravity_wave_driver_mod
