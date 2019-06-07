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
  use section_choice_config_mod,      only : cloud, cloud_none
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
  !> @param[out]   twod_fields Collection of two fields
  !> @param[out]   physics_incs Collection of physics increments
  subroutine create_physics_prognostics( mesh_id, twod_mesh_id, &
                                         prognostic_fields, &
                                         derived_fields, cloud_fields, &
                                         twod_fields, physics_incs )

    implicit none

    integer(i_def), intent(in)               :: mesh_id
    integer(i_def), intent(in)               :: twod_mesh_id

    ! Collections of fields
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(out) :: twod_fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: derived_fields
    type(field_collection_type), intent(out) :: physics_incs

    ! pointers to vector spaces
    type(function_space_type), pointer  :: vector_space => null()

    type( field_type ), pointer   :: theta => null()

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

    call add_physics_field(derived_fields, 'w_physics',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'rho_in_wth',     vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'wetrho_in_wth',  vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'exner_in_wth',   vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'w_physics_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'theta_star',     vector_space, checkpoint_restart_flag)

    ! W3 fields
    vector_space=> function_space_collection%get_fs(mesh_id, 0, W3)

    call add_physics_field(derived_fields, 'u1_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u2_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u3_in_w3',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'theta_in_w3',   vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'wetrho_in_w3',  vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u1_in_w3_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u2_in_w3_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u3_in_w3_star', vector_space, checkpoint_restart_flag)

    ! W2 fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, W2)

    call add_physics_field(derived_fields, 'u_physics',      vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u_physics_star', vector_space, checkpoint_restart_flag)
    call add_physics_field(derived_fields, 'u_star',         vector_space, checkpoint_restart_flag)


    !========================================================================
    ! Here we create some 2d fields for the UM physics
    !========================================================================
    twod_fields = field_collection_type(name='twod_fields')
    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3)
    checkpoint_restart_flag = .true.

    call add_physics_field(twod_fields, 'tstar',   vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'zh',      vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'z0msea',  vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'conv_rain',  vector_space, checkpoint_restart_flag, twod=.true.)

    checkpoint_restart_flag = .false.
    call add_physics_field(twod_fields, 'ntml',    vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'cumulus', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'cos_zenith_angle',   vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'lit_fraction',       vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'stellar_irradiance', vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'ls_rain',  vector_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(twod_fields, 'ls_snow',  vector_space, checkpoint_restart_flag, twod=.true.)

    !========================================================================
    ! Here we create some cloud fields
    !========================================================================

    cloud_fields = field_collection_type(name='cloud_fields')
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)
    checkpoint_restart_flag = .true.

    call add_physics_field(cloud_fields, 'area_fraction',   vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, 'ice_fraction',    vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, 'liquid_fraction', vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, 'bulk_fraction',   vector_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, 'rh_crit_wth',     vector_space, checkpoint_restart_flag)

    !========================================================================
    ! Increment values from individual physics parametrizations
    ! either for use by subsequent parametrizations or as diagnostics
    !========================================================================
    physics_incs = field_collection_type(name='physics_incs')
    vector_space => function_space_collection%get_fs(mesh_id, 0, Wtheta)
    checkpoint_restart_flag = .false. ! no need to dump any of these

    call add_physics_field(physics_incs, 'dt_bl', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'dmv_bl', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'dt_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'dmv_conv', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'dtl_mphys', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'dmt_mphys', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'sw_heating_rate', vector_space, checkpoint_restart_flag)
    call add_physics_field(physics_incs, 'lw_heating_rate', vector_space, checkpoint_restart_flag)

    ! Put references to fields requiring checkpointing into the prognostic fields collection
    call prognostic_fields%add_reference_to_field( twod_fields%get_field('tstar') )
    call prognostic_fields%add_reference_to_field( twod_fields%get_field('zh') )
    call prognostic_fields%add_reference_to_field( twod_fields%get_field('z0msea') )

    if (cloud /= cloud_none)then
      call prognostic_fields%add_reference_to_field( cloud_fields%get_field('area_fraction') )
      call prognostic_fields%add_reference_to_field( cloud_fields%get_field('ice_fraction') )
      call prognostic_fields%add_reference_to_field( cloud_fields%get_field('liquid_fraction') )
      call prognostic_fields%add_reference_to_field( cloud_fields%get_field('bulk_fraction') )
      call prognostic_fields%add_reference_to_field( cloud_fields%get_field('rh_crit_wth') )
    end if

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
  subroutine add_physics_field(field_collection, name, vector_space, &
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
      if (present(twod)) then
        if (twod) write_diag_behaviour => xios_write_field_single_face
      end if
      call new_field%set_write_behaviour(write_diag_behaviour)
    end if
    if (checkpoint_restart_flag) then
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    call field_collection%add_field(new_field)

  end subroutine add_physics_field

end module create_physics_prognostics_mod
