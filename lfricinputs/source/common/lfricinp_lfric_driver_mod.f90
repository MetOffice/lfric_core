! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_lfric_driver_mod

USE constants_mod,              ONLY: i_def, imdi, r_second
USE log_mod,                    ONLY: log_event, log_scratch_space,            &
                                      LOG_LEVEL_INFO, LOG_LEVEL_ERROR,         &
                                      LOG_LEVEL_ALWAYS

! LFRic Modules
USE driver_mesh_mod,            ONLY: init_mesh
USE driver_fem_mod,             ONLY: init_fem
USE driver_log_mod,             ONLY: init_logger, final_logger
USE derived_config_mod,         ONLY: set_derived_config
USE extrusion_mod,              ONLY: extrusion_type
USE field_collection_mod,       ONLY: field_collection_type
USE field_mod,                  ONLY: field_type
USE gungho_extrusion_mod,       ONLY: create_extrusion
USE halo_comms_mod,             ONLY: initialise_halo_comms
USE model_clock_mod,            ONLY: model_clock_type
USE lfric_xios_context_mod,     ONLY: lfric_xios_context_type, advance
USE lfric_xios_driver_mod,      ONLY: lfric_xios_initialise, &
                                      lfric_xios_finalise
USE lfricinp_setup_io_mod,      ONLY: init_lfricinp_files
USE linked_list_mod,            ONLY: linked_list_type
USE mesh_mod,                   ONLY: mesh_type
USE lfricinp_runtime_constants_mod, ONLY: lfricinp_create_runtime_constants
USE step_calendar_mod,          ONLY: step_calendar_type

! Interface to mpi
USE mpi_mod,                    ONLY: initialise_comm, store_comm,             &
                                      get_comm_size, get_comm_rank,            &
                                      finalise_comm

! lfricinp modules
USE lfricinp_um_parameters_mod, ONLY: fnamelen

IMPLICIT NONE

PRIVATE
PUBLIC :: lfricinp_initialise_lfric, lfricinp_finalise_lfric, lfric_fields,    &
          io_context

CHARACTER(len=fnamelen) :: xios_id
! xios_ctx names needs to match iodef.xml file
CHARACTER(len=*), PARAMETER :: xios_ctx  = "gungho_atm"
CHARACTER(len=fnamelen) :: program_name

! MPI ranks
INTEGER(KIND=i_def), PUBLIC :: total_ranks
INTEGER(KIND=i_def), PUBLIC :: local_rank

INTEGER(KIND=i_def), PUBLIC :: comm = -999

! Coordinate field
TYPE(field_type), TARGET :: chi(3)
TYPE(field_type), TARGET :: panel_id

TYPE(mesh_type), PUBLIC, pointer :: mesh      => null()
TYPE(mesh_type), PUBLIC, pointer :: twod_mesh => null()

! Container for all input fields
TYPE(field_collection_type) :: lfric_fields

TYPE(model_clock_type), PUBLIC, ALLOCATABLE :: model_clock
type(lfric_xios_context_type),  ALLOCATABLE :: io_context

CONTAINS

SUBROUTINE lfricinp_initialise_lfric(program_name_arg,                         &
                                     lfric_nl_fname,                           &
                                     required_lfric_namelists,                 &
                                     time_origin,                              &
                                     first_step, last_step,                    &
                                     spinup_period, seconds_per_step)

! Description:
!  Initialises LFRic infrastructure, MPI, XIOS and halos.

IMPLICIT NONE

CHARACTER(LEN=*),    INTENT(IN) :: program_name_arg
CHARACTER(LEN=*),    INTENT(IN) :: lfric_nl_fname
CHARACTER(LEN=*),    INTENT(IN) :: required_lfric_namelists(:)
CHARACTER(LEN=*),    INTENT(IN) :: time_origin
INTEGER(KIND=i_def), INTENT(IN) :: first_step, last_step
REAL(r_second),      INTENT(IN) :: spinup_period
REAL(r_second),      INTENT(IN) :: seconds_per_step

CLASS(extrusion_type),    ALLOCATABLE :: extrusion
TYPE(step_calendar_type), ALLOCATABLE :: model_calendar
TYPE(linked_list_type),   POINTER     :: file_list => null()

! Set module variables
program_name = program_name_arg
xios_id = TRIM(program_name) // "_client"

! Initialise MPI and create the default communicator: mpi_comm_world
CALL initialise_comm(comm)

! Initialise xios
call lfric_xios_initialise( program_name, comm, .false. )

! Save LFRic's part of the split communicator for later use, and
! set the total number of ranks and the local rank of the split
! communicator
CALL store_comm(comm)
total_ranks = get_comm_size()
local_rank = get_comm_rank()

!Initialise halo functionality
CALL initialise_halo_comms( comm )

CALL load_configuration(lfric_nl_fname, required_lfric_namelists)

! Initialise logging system
CALL init_logger( comm, program_name )

WRITE(log_scratch_space, '(2(A,I0))') 'total ranks = ', total_ranks,           &
                            ', local_rank = ', local_rank
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Sets variables used interally by the LFRic infrastructure.
CALL set_derived_config( .TRUE. )

CALL log_event('Initialising mesh', LOG_LEVEL_INFO)

! Generate prime mesh extrusion
ALLOCATE(extrusion, source=create_extrusion())

CALL init_mesh(local_rank, total_ranks, mesh, twod_mesh, &
               input_extrusion=extrusion)

! Create FEM specifics (function spaces and chi field)
CALL log_event('Creating function spaces and chi', LOG_LEVEL_INFO)
CALL init_fem(mesh, chi, panel_id)

! XIOS domain initialisation
model_calendar = step_calendar_type(time_origin)
model_clock = model_clock_type( first_step, last_step, seconds_per_step, &
                                spinup_period )

allocate( io_context )
file_list => io_context%get_filelist()
CALL init_lfricinp_files(file_list)
CALL io_context%initialise( xios_ctx, comm, chi, panel_id, model_clock, &
                            model_calendar )
call advance(io_context, model_clock)

! Initialise runtime constants
CALL log_event('Initialising runtime constants', LOG_LEVEL_INFO)
CALL lfricinp_create_runtime_constants(mesh,                                   &
                              twod_mesh,                                       &
                              chi,                                             &
                              panel_id,                                        &
                              seconds_per_step)

END SUBROUTINE lfricinp_initialise_lfric

!------------------------------------------------------------------

SUBROUTINE load_configuration(lfric_nl, required_lfric_namelists)

! Description:
!  Reads lfric namelists and checks that all required namelists are present

USE configuration_mod, ONLY: read_configuration, ensure_configuration

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: lfric_nl

CHARACTER(*), INTENT(IN)  :: required_lfric_namelists(:)

LOGICAL              :: okay
LOGICAL, ALLOCATABLE :: success_map(:)
INTEGER              :: i

ALLOCATE(success_map(SIZE(required_lfric_namelists)))

CALL log_event('Loading '//TRIM(program_name)//' configuration ...',           &
               LOG_LEVEL_ALWAYS)

CALL read_configuration( lfric_nl )

okay = ensure_configuration(required_lfric_namelists, success_map)
IF (.NOT. okay) THEN
  WRITE(log_scratch_space, '(A)')                                              &
                         'The following required namelists were not loaded:'
  DO i = 1, SIZE(required_lfric_namelists)
    IF (.NOT. success_map(i))                                                  &
      log_scratch_space = TRIM(log_scratch_space) // ' '                       &
                          // required_lfric_namelists(i)
  END DO
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

DEALLOCATE(success_map)

END SUBROUTINE load_configuration

!-------------------------------------------------------------------------------

SUBROUTINE lfricinp_finalise_lfric()

! Description:
!  Call finalise routines for associated APIs and logging system

USE halo_comms_mod,            ONLY: finalise_halo_comms
USE log_mod,                   ONLY: log_event, LOG_LEVEL_INFO

! External libraries
USE xios,                      ONLY: xios_finalize
USE mpi_mod,                   ONLY: finalise_comm


IMPLICIT NONE

CALL log_event( 'Calling lfric finalise routines', LOG_LEVEL_INFO )

! Finalise halos, XIOS, etc.
CALL finalise_halo_comms()
deallocate( io_context )
CALL lfric_xios_finalise()

! Finalise the logging system. This has to be done before finallising MPI
! as logging is an MPI process.
!
CALL final_logger(program_name)

CALL finalise_comm()

END SUBROUTINE lfricinp_finalise_lfric

END MODULE lfricinp_lfric_driver_mod
