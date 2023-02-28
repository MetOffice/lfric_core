! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM um2lfric

! lfricinputs modules
USE lfricinp_read_command_line_args_mod, ONLY: lfricinp_read_command_line_args
USE lfricinp_um_parameters_mod, ONLY: fnamelen
USE lfricinp_lfric_driver_mod, ONLY: lfricinp_initialise_lfric,                &
                lfricinp_finalise_lfric, mesh, twod_mesh, lfric_fields
USE lfricinp_ancils_mod, ONLY: lfricinp_create_ancil_fields, ancil_fields
USE lfricinp_create_lfric_fields_mod, ONLY: lfricinp_create_lfric_fields
USE lfricinp_um_grid_mod, ONLY: um_grid
USE lfricinp_initialise_um_mod, ONLY: lfricinp_initialise_um,                  &
    lfricinp_finalise_um, um_input_file
USE lfricinp_regrid_options_mod, ONLY: lfricinp_init_regrid_options
USE lfricinp_datetime_mod, ONLY : datetime_type
USE lfricinp_read_um_time_data_mod, ONLY: lfricinp_read_um_time_data
USE lfricinp_setup_io_mod,          ONLY: init_io_setup

! um2lfric modules
USE um2lfric_namelist_mod, ONLY: um2lfric_config, required_lfric_namelists
USE um2lfric_initialise_um2lfric_mod, ONLY: um2lfric_initialise_um2lfric
USE um2lfric_regrid_weights_mod, ONLY: um2lfric_regrid_weightsfile_ctl
USE um2lfric_init_masked_field_adjustments_mod,  ONLY:                         &
                                       um2lfric_init_masked_field_adjustments
USE um2lfric_regrid_and_output_data_mod, ONLY: um2lfric_regrid_and_output_data

! LFRic modules
USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO

IMPLICIT NONE

TYPE(datetime_type)        :: datetime
CHARACTER(LEN=fnamelen)    :: lfric_fname, um2lfric_fname, io_fname

CALL lfricinp_read_command_line_args(um2lfric_fname, lfric_fname, io_fname)

! Set up IO file configuration
CALL init_io_setup(io_fname)

! Read um2lfric configuration namelist
CALL um2lfric_config%load_namelist(um2lfric_fname)

! Read in global regrid options
CALL lfricinp_init_regrid_options(um2lfric_fname)

! Open the UM file
CALL log_event('Initialising UM input file', LOG_LEVEL_INFO)
CALL lfricinp_initialise_um(um2lfric_config%um_file)

! Load date and time information for requested stash items from um input file
CALL datetime % initialise()
CALL lfricinp_read_um_time_data(datetime, um_input_file,                       &
                                um2lfric_config%stash_list)

! Initialise LFRic Infrastructure
CALL lfricinp_initialise_lfric(program_name_arg="um2lfric",                    &
     lfric_nl_fname=lfric_fname,                                               &
     required_lfric_namelists = required_lfric_namelists,                      &
     time_origin = datetime % first_validity_time,                             &
     first_step = datetime % first_step,                                       &
     last_step = datetime % last_step,                                         &
     spinup_period = datetime % spinup_period,                                 &
     seconds_per_step = datetime % seconds_per_step)

! Initialise um2lfric
CALL um2lfric_initialise_um2lfric()

! Initialise LFRic field collection
CALL lfricinp_create_lfric_fields( mesh, twod_mesh, lfric_fields, &
                                   um2lfric_config%stash_list,    &
                                   um_grid, um_input_file )

! Initialise LFRic ancils field collection
CALL lfricinp_create_ancil_fields( ancil_fields, mesh, twod_mesh )

! Read in, process and partition regridding weights
CALL um2lfric_regrid_weightsfile_ctl()

! Now initialise masked points that requires post regridding adjustments
CALL um2lfric_init_masked_field_adjustments()

! Perform regridding and output data
CALL um2lfric_regrid_and_output_data(datetime)

! Unloads data from memory and closes UM input file
CALL lfricinp_finalise_um()

! Finalise YAXT, XIOS, MPI, logging
CALL lfricinp_finalise_lfric()

END PROGRAM um2lfric
