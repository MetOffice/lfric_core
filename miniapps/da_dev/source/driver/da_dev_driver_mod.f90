!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the da_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module da_dev_driver_mod

  use checksum_alg_mod,         only: checksum_alg
  use constants_mod,            only: i_def, i_native, &
                                      PRECISION_REAL, r_def
  use clock_mod,                only: clock_type
  use driver_model_data_mod,    only: model_data_type
  use driver_log_mod,           only: init_logger, final_logger
  use driver_time_mod,          only: init_time, get_calendar
  use driver_mesh_mod,          only: init_mesh, final_mesh
  use driver_fem_mod,           only: init_fem, final_fem
  use driver_io_mod,            only: init_io, final_io, filelist_populator, &
                                      get_io_context
  use io_config_mod,            only: use_xios_io
  use io_context_mod,           only: io_context_type
  use field_mod,                only: field_type
  use da_dev_model_init_mod,    only: create_da_model_data, initialise_da_model_data
  use log_mod,                  only: log_event,          &
                                      log_scratch_space,  &
                                      LOG_LEVEL_ALWAYS,   &
                                      LOG_LEVEL_INFO
  use mesh_mod,                 only: mesh_type
  use extrusion_mod,            only: extrusion_type
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: mpi_type
  use da_dev_increment_alg_mod, only: da_dev_increment_alg
  use da_dev_init_files_mod,    only: init_da_dev_files
  use da_dev_config_mod,        only: write_data, test_field
  use da_dev_extrusion_mod,     only: create_extrusion
  use lfric_da_driver_mod,      only: init_da, final_da
  use lfric_xios_read_mod,      only: read_state
  use lfric_xios_write_mod,     only: write_state
  use lfric_xios_context_mod,   only: lfric_xios_context_type, advance
  use field_collection_mod,     only: field_collection_type

  implicit none

  private
  ! LFRic-API
  public initialise_lfric, finalise_lfric, step_lfric
  ! To run local LFRic mini-app
  public run, initialise_model, finalise_model

  type(model_clock_type), allocatable, public :: model_clock

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(mesh_type),  pointer, public      :: mesh      => null()
  type(mesh_type),  pointer, public      :: twod_mesh => null()

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise_lfric( program_name, mpi )

    implicit none

    character(len=*),  intent(in)    :: program_name
    class(mpi_type),   intent(inout) :: mpi

    procedure(filelist_populator), pointer :: fl_populator => null()
    class(io_context_type),        pointer :: model_io_context => null()
    class(extrusion_type),     allocatable :: extrusion

    call init_logger( mpi%get_comm(), program_name )

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create model clock and calendar
    call init_time( model_clock )

    ! Create the mesh
    allocate( extrusion, source=create_extrusion() )
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(), &
                    mesh, twod_mesh = twod_mesh, input_extrusion = extrusion )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Initialise I/O context
    fl_populator => init_da_dev_files
    call init_io( program_name, mpi%get_comm(), chi, panel_id, &
                  model_clock, get_calendar(), populate_filelist=fl_populator )

    ! Do initial step
    model_io_context => get_io_context()
    if (model_clock%is_initialisation()) then
      select type (model_io_context)
      type is (lfric_xios_context_type)
          call advance(model_io_context, model_clock)
      end select
    end if

  end subroutine initialise_lfric
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> initialise the model data (create and read/initialise).
  !>
  subroutine initialise_model( mpi, model_data )

    implicit none

    class(mpi_type),       intent(inout) :: mpi
    type(model_data_type), intent(inout) :: model_data

    ! Create prognostic fields in model_data
    call create_da_model_data( mesh, twod_mesh, model_data )

    ! Initialise DA interface
    call init_da( mpi%get_comm(), chi, panel_id, model_clock, get_calendar() )

    ! Initialise prognostic fields in model_data
    call initialise_da_model_data( model_data )

  end subroutine initialise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run( program_name, model_data )

    implicit none

    character(len=*), intent(in)         :: program_name
    type(model_data_type), intent(inout) :: model_data

    call log_event(program_name//": Run begins", LOG_LEVEL_INFO)

    do while( model_clock%tick() )
      call step_lfric( model_data )
    end do

    call log_event(program_name//": Run ends", LOG_LEVEL_INFO)

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs a single time step.
  !>
  subroutine step_lfric( model_data )

    implicit none
    type(model_data_type), intent(inout) :: model_data

    type(field_type),       pointer :: working_field    => null()
    class(io_context_type), pointer :: model_io_context => null()

    call model_data%depository%get_field( test_field, working_field )

    ! Switch to the main model I/O context
    model_io_context => get_io_context()
    call model_io_context%set_current()

    call da_dev_increment_alg( working_field )

    if ( write_data ) then
      call write_state( model_data%depository, prefix="write_" )
    end if

  end subroutine step_lfric

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finalise the model data after a run.
  subroutine finalise_model( program_name, field_collection  )

    implicit none

    character(len=*), intent(in)               :: program_name
    type(field_collection_type), intent(inout) :: field_collection

    type(field_type), pointer :: working_field => null()

    call field_collection%get_field( test_field, working_field )

    !---------------------------------------------------------------------------
    ! Model finalise
    !---------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg( program_name, working_field, test_field )

  end subroutine finalise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise_lfric( program_name )

    implicit none
    character(len=*), intent(in) :: program_name

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    call final_da()

    call final_io()

    call final_fem()

    call final_mesh()

    call final_logger( program_name )

  end subroutine finalise_lfric

end module da_dev_driver_mod
