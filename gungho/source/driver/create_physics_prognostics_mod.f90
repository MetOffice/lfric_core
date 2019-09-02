!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type, &
                                             write_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,         &
                                             LOG_LEVEL_ERROR
  use radiation_config_mod,           only : n_radstep
  use section_choice_config_mod,      only : cloud, cloud_none, &
                                             radiation, radiation_none
  use time_config_mod,                only : timestep_start, timestep_end
  implicit none

  private
  public :: create_physics_prognostics

contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[inout] prognostic_fields A collection of the fields that make up the
  !>                                 prognostic variables in the model
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields
  !> @param[out]   cloud_fields Collection of FD cloud fields
  !> @param[out]   twod_fields Collection of two dimensional fields
  !> @param[out]   radstep_fields Collection of radiation timestep fields
  !> @param[out]   physics_incs Collection of physics increments
  !> @param[out]   jules_ancils Ancillary fields for Jules
  !> @param[out]   jules_prognostics Prognostic fields for Jules
  subroutine create_physics_prognostics( mesh_id, twod_mesh_id, &
                                         depository, &
                                         prognostic_fields, &
                                         derived_fields, cloud_fields, &
                                         twod_fields, radstep_fields, &
                                         physics_incs, &
                                         jules_ancils, jules_prognostics )

    implicit none

    integer(i_def), intent(in) :: mesh_id
    integer(i_def), intent(in) :: twod_mesh_id

    ! Collections of fields
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(out) :: twod_fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: derived_fields
    type(field_collection_type), intent(out) :: radstep_fields
    type(field_collection_type), intent(out) :: physics_incs
    type(field_collection_type), intent(out) :: jules_ancils
    type(field_collection_type), intent(out) :: jules_prognostics

    ! pointers to vector spaces
    type(function_space_type), pointer :: vector_space => null()

    type( field_type ), pointer :: theta => null()
 
    ! Each column of a higher-order discontinuous field will be used to
    ! represent multi-dimensional quantities like tiles, plant functional
    ! types and sea ice categories. Set parameters for the orders required:
    integer(i_def) :: tile_order = 2 ! Enough space for 27 tiles
    integer(i_def) :: pft_order  = 1 ! Enough space for 8 plant functional types
    integer(i_def) :: sice_order = 1 ! Enough space for 8 sea ice categories
    integer(i_def) :: soil_order = 1 ! Enough space for 8 soil levels

    integer(i_def) :: theta_space
    logical(l_def) :: checkpoint_restart_flag

    call log_event( 'Create physics prognostics', LOG_LEVEL_INFO )


    theta => prognostic_fields%get_field('theta')
    theta_space=theta%which_function_space()

    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta', LOG_LEVEL_ERROR )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements', LOG_LEVEL_ERROR )
    end if

    !========================================================================
    ! Here we create some field collections
    !========================================================================
    derived_fields  =  field_collection_type(name='derived_fields')
    checkpoint_restart_flag = .false.

    ! Wtheta fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)

    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'w_physics',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'rho_in_wth',     vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'wetrho_in_wth',  vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'exner_in_wth',   vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'w_physics_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'theta_star',     vector_space, checkpoint_restart_flag)

    ! W3 fields
    vector_space=> function_space_collection%get_fs(mesh_id, 0, W3)

    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u1_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u2_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u3_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'theta_in_w3',   vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'wetrho_in_w3',  vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u1_in_w3_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u2_in_w3_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u3_in_w3_star', vector_space, checkpoint_restart_flag)

    ! W2 fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, W2)

    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u_physics',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u_physics_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, depository, prognostic_fields, &
      'u_star',         vector_space, checkpoint_restart_flag)


    !========================================================================
    ! Here we create some 2d fields for the UM physics
    !========================================================================
    twod_fields = field_collection_type(name='twod_fields')
    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3)
    checkpoint_restart_flag = .true.

    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'tstar',   vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'zh',      vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'z0msea',  vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'conv_rain',  vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'conv_snow',  vector_space, checkpoint_restart_flag, twod=.true.)

    checkpoint_restart_flag = .false.
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'ntml',    vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'cumulus', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'cos_zenith_angle',   vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'lit_fraction',       vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'ls_rain',  vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, depository, prognostic_fields, &
      'ls_snow',  vector_space, checkpoint_restart_flag, twod=.true.)

    !========================================================================
    ! Here we create some cloud fields
    !========================================================================

    cloud_fields = field_collection_type(name='cloud_fields')
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)
    if (cloud /= cloud_none)then
      checkpoint_restart_flag = .true.
    end if

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'area_fraction',   vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'ice_fraction',    vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'liquid_fraction', vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'bulk_fraction',   vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'rh_crit_wth',     vector_space, checkpoint_restart_flag)
    ! convective cloud field prognostics
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'cca', vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'ccw', vector_space, checkpoint_restart_flag)

    !========================================================================
    ! Radiation timestep fields
    !========================================================================
    radstep_fields = field_collection_type(name='radstep_fields')
    if (radiation /= radiation_none) then
      ! Checkpoint unless both the first timestep of this run and the
      ! first timestep of the next run are radiation timesteps
      checkpoint_restart_flag = &
        mod(timestep_start-1, n_radstep) /= 0 .or. &
        mod(timestep_end,     n_radstep) /= 0
    else
      checkpoint_restart_flag = .false.
    end if

    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_heating_rate_rts', vector_space, checkpoint_restart_flag)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'lw_heating_rate_rts', vector_space, checkpoint_restart_flag)

    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3) 
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'lw_down_surf_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_down_surf_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_direct_surf_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_down_blue_surf_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_direct_blue_surf_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'cos_zenith_angle_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'lit_fraction_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'stellar_irradiance_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sin_stellar_declination_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'stellar_eqn_of_time_rts', vector_space, checkpoint_restart_flag, &
      twod=.true.)

    ! Multidimensional fields currently not able to be checkpointed.
    checkpoint_restart_flag = .false.
    vector_space=>function_space_collection%get_fs(twod_mesh_id, tile_order, W3)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'lw_up_tile_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_up_tile_rts', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(radstep_fields, depository, prognostic_fields, &
      'sw_up_blue_tile_rts', vector_space, checkpoint_restart_flag, twod=.true.)

    !========================================================================
    ! Increment values from individual physics parametrizations
    ! either for use by subsequent parametrizations or as diagnostics
    !========================================================================
    physics_incs = field_collection_type(name='physics_incs')
    checkpoint_restart_flag = .false. ! no need to dump any of these

    vector_space => function_space_collection%get_fs(mesh_id, 0, Wtheta)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dt_bl', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dmv_bl', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dt_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dmv_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dmcl_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dmcf_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dtl_mphys', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dmt_mphys', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'sw_heating_rate', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'lw_heating_rate', vector_space, checkpoint_restart_flag)

    ! Increments on pgrid (not U and V grids) but rho levels
    vector_space => function_space_collection%get_fs(mesh_id, 0, W3)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'du_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'dv_conv', vector_space, checkpoint_restart_flag)

    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3) 
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'lw_down_surf', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'sw_down_surf', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'sw_direct_surf', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'sw_down_blue_surf', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(physics_incs, depository, prognostic_fields, &
      'sw_direct_blue_surf', vector_space, checkpoint_restart_flag, twod=.true.)


    !========================================================================
    ! Surface fields (a temporary treatment for multi-dimensional fields
    !                 using higher-order fields)
    !========================================================================

    ! Jules ancillaries
    jules_ancils = field_collection_type(name='jules_ancils')
    checkpoint_restart_flag = .false. ! #1920 required for this to work

    vector_space => function_space_collection%get_fs(twod_mesh_id, tile_order, W3)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'tile_fraction', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, pft_order, W3)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'leaf_area_index', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'canopy_height', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, 0, W3)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'sd_orog', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'peak_to_trough_orog', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'silhouette_area_orog', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_albedo', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_roughness', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_thermal_cond', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_carbon_content', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'albedo_obs_sw', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'albedo_obs_vis', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'albedo_obs_nir', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, soil_order, W3)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_moist_wilt', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_moist_crit', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_moist_sat', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_cond_sat', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_thermal_cap', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'soil_suction_sat', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_ancils, depository, prognostic_fields, &
      'clapp_horn_b', vector_space, checkpoint_restart_flag, twod=.true.)

    ! Jules prognostics
    jules_prognostics = field_collection_type(name='jules_prognostics')
    checkpoint_restart_flag = .false. ! #1920 required for this to work

    vector_space => function_space_collection%get_fs(twod_mesh_id, tile_order, W3)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'tile_temperature', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'tile_snow_mass', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'tile_snow_rgrain', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'n_snow_layers', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'snow_depth', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'canopy_water', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'lw_up_tile', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'sw_up_tile', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'sw_up_blue_tile', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, 0, W3)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'snow_soot', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'surface_conductance', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'chloro_sea', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, soil_order, W3)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'soil_temperature', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'soil_moisture', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'unfrozen_soil_moisture', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'frozen_soil_moisture', vector_space, checkpoint_restart_flag, twod=.true.)

    vector_space => function_space_collection%get_fs(twod_mesh_id, sice_order, W3)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'sea_ice_thickness', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'sea_ice_pond_frac', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(jules_prognostics, depository, prognostic_fields, &
      'sea_ice_pond_depth', vector_space, checkpoint_restart_flag, twod=.true.)

  end subroutine create_physics_prognostics

  !>@brief Add field to field collection and set its write,
  !>       checkpoint and restart behaviour
  !> @param[in,out] field_collection Field collection that 'name' will be added to
  !> @param[in]     name             Name of field to be added to collection
  !> @param[in]     vector_space     Function space of field to set behaviour for
  !> @param[in]     checkpoint_restart_flag  Flag to allow checkpoint and
  !>                                        restart behaviour of field to be set
  !> @param[in]     twod             Optional flag to determine if this is a
  !>                                        2D field for diagnostic output
  subroutine add_physics_field(field_collection, &
                               depository, prognostic_fields, &
                               name, vector_space, &
                               checkpoint_restart_flag, twod)

    use io_config_mod,      only : use_xios_io, &
                                   write_diag
    use io_mod,             only : xios_write_field_face,        &
                                   xios_write_field_single_face, &
                                   checkpoint_write_xios,        &
                                   checkpoint_write_netcdf,      &
                                   checkpoint_read_netcdf,       &
                                   checkpoint_read_xios

    implicit none

    character(*), intent(in)                   :: name
    type(field_collection_type), intent(inout) :: field_collection
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(function_space_type), intent(in)      :: vector_space
    logical(l_def), intent(in)                 :: checkpoint_restart_flag
    logical, optional, intent(in)              :: twod
    !Local variables
    type(field_type)                           :: new_field

    ! pointers for xios write interface
    procedure(write_interface), pointer   :: write_diag_behaviour => null()
    procedure(checkpoint_write_interface), pointer  :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer   :: checkpoint_read_behaviour => null()

    if ( use_xios_io) then
      checkpoint_write_behaviour => checkpoint_write_xios
      checkpoint_read_behaviour => checkpoint_read_xios
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour => checkpoint_read_netcdf
    endif

    new_field = field_type( vector_space, name=trim(name) )

    if (use_xios_io .and. write_diag) then
      ! All physics fields currently require output on faces...
      write_diag_behaviour => xios_write_field_face
      if (present(twod))then
        if (twod) write_diag_behaviour => xios_write_field_single_face
      end if
      call new_field%set_write_behaviour(write_diag_behaviour)
    end if
    if (checkpoint_restart_flag) then
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( depository%get_field(name) )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if (checkpoint_restart_flag) then
      call prognostic_fields%add_reference_to_field( depository%get_field(name) )
    endif

  end subroutine add_physics_field

end module create_physics_prognostics_mod
