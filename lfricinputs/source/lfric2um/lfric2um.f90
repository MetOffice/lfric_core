! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM lfric2um

! lfricinputs modules
USE lfricinp_read_command_line_args_mod, ONLY: lfricinp_read_command_line_args
USE lfricinp_create_lfric_fields_mod,  ONLY: lfricinp_create_lfric_fields
USE lfricinp_setup_io_mod,  ONLY: init_io_setup

! lfric2um modules
USE lfric2um_namelists_mod, ONLY: lfric2um_nl_fname, lfric_nl_fname,           &
     io_nl_fname, lfric2um_config, required_lfric_namelists
USE lfricinp_lfric_driver_mod, ONLY: lfricinp_initialise_lfric,     &
                                     lfricinp_finalise_lfric, mesh, &
                                     twod_mesh,                     &
                                     lfric_fields
USE lfric2um_initialise_um_mod, ONLY: lfric2um_initialise_um, um_output_file
USE lfric2um_initialise_lfric2um_mod, ONLY: lfric2um_initialise_lfric2um
USE lfric2um_main_loop_mod, ONLY: lfric2um_main_loop
USE lfricinp_um_grid_mod, ONLY: um_grid
USE lfricinp_regrid_options_mod, ONLY: lfricinp_init_regrid_options
USE lfricinp_datetime_mod, ONLY : datetime_type

IMPLICIT NONE

TYPE(datetime_type) :: datetime

! Read command line args
CALL lfricinp_read_command_line_args(lfric2um_nl_fname, lfric_nl_fname,        &
                                     io_nl_fname)

! Set up IO file configuration
CALL init_io_setup(io_nl_fname)

! Read in global regrid options
CALL lfricinp_init_regrid_options(lfric2um_nl_fname)

! Load date and time information
CALL datetime % initialise()

! Initialise LFRic Infrastructure
CALL lfricinp_initialise_lfric(program_name_arg="lfric2um",                    &
     lfric_nl_fname=lfric_nl_fname,                                            &
     required_lfric_namelists = required_lfric_namelists,                      &
     time_origin = datetime % first_validity_time,                             &
     first_step = datetime % first_step,                                       &
     last_step = datetime % last_step,                                         &
     spinup_period = datetime % spinup_period,                                 &
     seconds_per_step = datetime % seconds_per_step)

! Initialise lfric2um
CALL lfric2um_initialise_lfric2um()

! Initialise UM Infrastructure
CALL lfric2um_initialise_um()

! Create LFRic field collection based on list of stashcodes
CALL lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields,       &
                                   lfric2um_config%stash_list, um_grid, &
                                   um_output_file )

! Main loop over fields to be read, regridded and written to output dump
CALL lfric2um_main_loop()

! Finalise LFRic infrastructure
CALL lfricinp_finalise_lfric()

END PROGRAM lfric2um
