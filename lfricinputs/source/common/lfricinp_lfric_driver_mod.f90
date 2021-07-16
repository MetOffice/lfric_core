! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_lfric_driver_mod

USE constants_mod,              ONLY: i_def, imdi, r_second
USE log_mod,                    ONLY: log_event, log_scratch_space,            &
                                      LOG_LEVEL_INFO, LOG_LEVEL_ERROR,         &
                                      LOG_LEVEL_ALWAYS, initialise_logging,    &
                                      finalise_logging

! LFRic Modules
USE lfric_xios_io_mod,          ONLY: initialise_xios
USE clock_mod,                  ONLY: clock_type
USE create_mesh_mod,            ONLY: init_mesh
USE create_fem_mod,             ONLY: init_fem
USE derived_config_mod,         ONLY: set_derived_config
USE field_collection_mod,       ONLY: field_collection_type
USE field_mod,                  ONLY: field_type
USE io_context_mod,             ONLY: io_context_type
USE mod_wait,                   ONLY: init_wait
USE linked_list_mod,            ONLY: linked_list_type
USE lfric_xios_io_mod,          ONLY: populate_filelist_if
USE lfricinp_setup_io_mod,      ONLY: init_lfricinp_files
USE local_mesh_collection_mod,  ONLY: local_mesh_collection,                   &
                                      local_mesh_collection_type
USE mesh_collection_mod,        ONLY: mesh_collection,                         &
                                      mesh_collection_type

! Interface to mpi
USE mpi_mod,                    ONLY: initialise_comm, store_comm,             &
                                      get_comm_size, get_comm_rank,            &
                                      finalise_comm
! External libs
USE xios,                       ONLY: xios_finalize, xios_initialize
USE yaxt,                       ONLY: xt_initialize

! lfricinp modules
USE lfricinp_um_parameters_mod, ONLY: fnamelen

IMPLICIT NONE

PRIVATE
PUBLIC :: lfricinp_initialise_lfric, lfricinp_finalise_lfric, lfric_fields

CHARACTER(len=fnamelen) :: xios_id
! xios_ctx names needs to match iodef.xml file
CHARACTER(len=*), PARAMETER :: xios_ctx  = "gungho_atm"
CHARACTER(len=fnamelen) :: program_name

! MPI ranks
INTEGER(KIND=i_def), PUBLIC :: total_ranks
INTEGER(KIND=i_def), PUBLIC :: local_rank

INTEGER(KIND=i_def), PUBLIC :: comm = -999

! Coordinate field
TYPE(field_type), TARGET :: chi_xyz(3)
TYPE(field_type), TARGET :: chi_sph(3)
TYPE(field_type), TARGET :: panel_id

INTEGER(KIND=i_def), PUBLIC :: mesh_id      = imdi
INTEGER(KIND=i_def), PUBLIC :: twod_mesh_id = imdi

! Container for all input fields
TYPE(field_collection_type) :: lfric_fields

CLASS(io_context_type), allocatable :: io_context

CONTAINS

SUBROUTINE populate_file_list( file_list, clock )
  ! Description:
  !   Populates I/O context's list of interesting files.

  IMPLICIT NONE

  CLASS(linked_list_type), INTENT(INOUT) :: file_list
  CLASS(clock_type),       INTENT(IN)    :: clock

  CALL init_lfricinp_files( file_list, clock )

END SUBROUTINE populate_file_list


SUBROUTINE lfricinp_initialise_lfric(program_name_arg, &
                                     lfric_nl_fname,   &
                                     required_lfric_namelists)

! Description:
!  Initialises LFRic infrastructure, MPI, XIOS and YAXT.

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: program_name_arg
CHARACTER(LEN=*), INTENT(IN) :: lfric_nl_fname
CHARACTER(LEN=*), INTENT(IN) :: required_lfric_namelists(:)

CHARACTER(*),   PARAMETER :: first_step       = '1'
CHARACTER(*),   PARAMETER :: last_step        = '1'
REAL(r_second), PARAMETER :: spinup_period    = 0.0_r_second
REAL(r_second), PARAMETER :: seconds_per_step = 1.0_r_second

CLASS(clock_type), POINTER :: clock
LOGICAL                    :: advance

PROCEDURE(populate_filelist_if), POINTER :: populate_pointer

! Set module variables
program_name = program_name_arg
xios_id = TRIM(program_name) // "_client"

! Initialise MPI, yaxt and XIOS
! Initialise MPI and create the default communicator: mpi_comm_world
CALL initialise_comm(comm)
CALL init_wait()
CALL xios_initialize(xios_id, return_comm = comm)
! Save LFRic's part of the split communicator for later use
CALL store_comm(comm)
CALL xt_initialize(comm)

! Get rank information
total_ranks = get_comm_size()
local_rank = get_comm_rank()

! Initialise logging system
CALL initialise_logging(local_rank, total_ranks, program_name)

WRITE(log_scratch_space, '(2(A,I0))') 'total ranks = ', total_ranks, &
                            ', local_rank = ', local_rank
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

CALL log_event('Loading LFRic Infrastructure namelists', LOG_LEVEL_INFO)
CALL load_configuration(lfric_nl_fname, required_lfric_namelists)

! Sets variables used interally by the LFRic infrastructure.
CALL set_derived_config( .TRUE. )

CALL log_event('Initialising mesh', LOG_LEVEL_INFO)
ALLOCATE(local_mesh_collection, source = local_mesh_collection_type())
ALLOCATE(mesh_collection, source=mesh_collection_type() )

CALL init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

! Create FEM specifics (function spaces and chi field)
CALL log_event('Creating function spaces and chi', LOG_LEVEL_INFO)
CALL init_fem(mesh_id, chi_xyz, chi_sph, panel_id)

! XIOS domain initialisation
populate_pointer => populate_file_list
call initialise_xios( io_context,       &
                      xios_ctx,         &
                      comm,             &
                      mesh_id,          &
                      twod_mesh_id,     &
                      chi_sph,          &
                      panel_id,         &
                      first_step,       &
                      last_step,        &
                      spinup_period,    &
                      seconds_per_step, &
                      populate_filelist=populate_pointer )

! Advance XIOS calendar to enable dump writing
clock => io_context%get_clock()
advance = clock%tick()

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

CALL log_event('Loading '//TRIM(program_name)//' configuration ...',                 &
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

USE log_mod,                   ONLY: finalise_logging, LOG_LEVEL_INFO, &
                                     log_event
! External libraries
USE xios,                      ONLY: xios_finalize
USE yaxt,                      ONLY: xt_finalize
USE mpi_mod,                   ONLY: finalise_comm


IMPLICIT NONE

CALL log_event( 'Calling lfric finalise routines', LOG_LEVEL_INFO )

! Finalise YAXT, XIOS, MPI, etc.
CALL xt_finalize()
CALL xios_finalize()
CALL finalise_comm()

! Finalise the logging system
CALL finalise_logging()

END SUBROUTINE lfricinp_finalise_lfric


END MODULE lfricinp_lfric_driver_mod
