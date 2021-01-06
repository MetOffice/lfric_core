!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_ancils_mod

  use constants_mod,                  only : i_def, l_def, dp_xios, str_def, r_def
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use field_mod,                      only : field_type
  use field_parent_mod,               only : read_interface, &
                                             write_interface
  use io_config_mod,                  only : use_xios_io
  use initialization_config_mod,      only : ancil_option,            &
                                             ancil_option_basic_gal,  &
                                             ancil_option_prototype_gal
  use linked_list_mod,                only : linked_list_type
  use read_methods_mod,               only : read_field_face, &
                                             read_field_single_face, &
                                             read_field_time_var, &
                                             read_time_data
  use write_methods_mod,              only : write_field_face, &
                                             write_field_single_face
  use field_collection_mod,           only : field_collection_type, &
                                             field_collection_real_iterator_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3, WTheta
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use time_axis_mod,                  only : time_axis_type, update_interface
  use jules_control_init_mod,         only: n_land_tile
  use jules_surface_types_mod,        only: npft

  implicit none

  public   :: create_fd_ancils,         &
              setup_ancil_field

contains

  !> @details Organises fields to be read from ancils into ancil_fields
  !           collection then reads them.
  !> @param[in,out] depository The depository field collection
  !> @param[out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh_id The identifier given to the current 3d mesh
  !> @param[in] twod_mesh_id The identifier given to the current 2d mesh
  subroutine create_fd_ancils( depository, ancil_fields, mesh_id, &
                               twod_mesh_id, ancil_times_list )

    implicit none

    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( out )   :: ancil_fields
    integer(i_def), intent(in) :: mesh_id
    integer(i_def), intent(in) :: twod_mesh_id
    type(linked_list_type), intent(out) :: ancil_times_list

    ! Pointer to time-axis update procedure
    procedure(update_interface), pointer :: tmp_update_ptr => null()

    ! Time axis objects for different ancil groups - must be saved to be
    ! available after function call
    type(time_axis_type), save :: sea_time_axis
    type(time_axis_type), save :: sst_time_axis
    type(time_axis_type), save :: sea_ice_time_axis
    type(time_axis_type), save :: aerosol_time_axis
    type(time_axis_type), save :: albedo_vis_time_axis
    type(time_axis_type), save :: albedo_nir_time_axis
    type(time_axis_type), save :: pft_time_axis
    type(time_axis_type), save :: ozone_time_axis

    ! Set pointer to time axis read behaviour
    tmp_update_ptr => read_field_time_var

    ! Set up ancil_fields collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Setting up ancil field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    ancil_fields = field_collection_type(name='ancil_fields')

    ! Here ancil fields are set up with a call to setup_ancil_field. For ancils
    ! that are time-varying, the time-axis (set up with a call to
    ! init_time_axis) is passed to the setup_ancil_field subroutine.

    ! We populate the ancil fields collection based on the model configuration

    ! Set up ancils for the basic GAL configuration (proto GAL also
    ! includes these):
    if ( ancil_option == ancil_option_basic_gal .or. &
         ancil_option == ancil_option_prototype_gal ) then

      !=====  SURFACE ANCILS  =====
      call setup_ancil_field("land_area_fraction", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("land_tile_fraction", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true., ndata=n_land_tile)
      call init_time_axis("pft_time", pft_time_axis)
      call setup_ancil_field("canopy_height_in", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
      call setup_ancil_field("leaf_area_index_in", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
      call pft_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(pft_time_axis)

      call init_time_axis("sea_time", sea_time_axis)
      call setup_ancil_field("chloro_sea", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true., time_axis=sea_time_axis)
      call sea_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(sea_time_axis)

      call init_time_axis("sst_time", sst_time_axis)
      call setup_ancil_field("tstar_sea", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true., time_axis=sst_time_axis)
      call sst_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(sst_time_axis)

      !=====  SEA ICE ANCILS  =====
      call init_time_axis("sea_ice_time", sea_ice_time_axis)
      call setup_ancil_field("sea_ice_thickness", depository, ancil_fields, &
                mesh_id, twod_mesh_id, twod=.true., time_axis=sea_ice_time_axis)
      call setup_ancil_field("sea_ice_fraction", depository, ancil_fields, &
                mesh_id, twod_mesh_id, twod=.true., time_axis=sea_ice_time_axis)
      call sea_ice_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(sea_ice_time_axis)

      !=====  RADIATION ANCILS  =====
      call init_time_axis("albedo_vis_time", albedo_vis_time_axis)
      call setup_ancil_field("albedo_obs_vis", depository, ancil_fields, mesh_id, &
                        twod_mesh_id, twod=.true., time_axis=albedo_vis_time_axis)
      call albedo_vis_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(albedo_vis_time_axis)

      call init_time_axis("albedo_nir_time", albedo_nir_time_axis)
      call setup_ancil_field("albedo_obs_nir", depository, ancil_fields, mesh_id, &
                        twod_mesh_id, twod=.true., time_axis=albedo_nir_time_axis)
      call albedo_nir_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(albedo_nir_time_axis)

      !=====  SOIL ANCILS  =====
      call setup_ancil_field("soil_albedo", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_carbon_content", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_thermal_cond", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_moist_wilt", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_moist_crit", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_moist_sat", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_cond_sat", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_thermal_cap", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("soil_suction_sat", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("clapp_horn_b", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("mean_topog_index", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("stdev_topog_index", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)

      !=====  OROGRAPHY ANCILS  =====
      call setup_ancil_field("sd_orog", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true.)
      call setup_ancil_field("grad_xx_orog", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true.)
      call setup_ancil_field("grad_xy_orog", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true.)
      call setup_ancil_field("grad_yy_orog", depository, ancil_fields, mesh_id, &
                              twod_mesh_id, twod=.true.)
      call setup_ancil_field("peak_to_trough_orog", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)
      call setup_ancil_field("silhouette_area_orog", depository, ancil_fields, &
                              mesh_id, twod_mesh_id, twod=.true.)

      !=====  OZONE ANCIL  =====
      call init_time_axis("ozone_time", ozone_time_axis)
      call setup_ancil_field("ozone", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=ozone_time_axis)
      call ozone_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(ozone_time_axis)

      !=====  AEROSOL ANCILS  =====
      call init_time_axis("aerosols_time", aerosol_time_axis)
      call setup_ancil_field("acc_sol_bc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_oc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_su", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_ss", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("nd_acc_sol", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_bc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_oc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_su", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("nd_ait_sol", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_ins_bc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_ins_oc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("nd_ait_ins", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_bc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_oc", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_su", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_ss", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call setup_ancil_field("nd_cor_sol", depository, ancil_fields, mesh_id, &
                             twod_mesh_id, time_axis=aerosol_time_axis)
      call aerosol_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(aerosol_time_axis)
    end if

    ! Now the field collection is set up, the fields will be initialised in
    ! gungho_model_data_mod

  end subroutine create_fd_ancils

  !> @details Adds fields to the ancil collection, sets up their read and write
  !>      behaviour and creates them in the depository if they do not yet exist
  !> @param[in] name The field name
  !> @param[in, out] depository The depository field collection
  !> @param[in, out] ancil_fields The ancil field collection
  !> @param[in] mesh_id The identifier given to the current 3d mesh
  !> @param[in, optional] twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in, optional] ndata Number of non-spatial dimensions for multi-data field
  !> @param[in, out, optional] time_axis Time axis associated with ancil field
  subroutine setup_ancil_field( name, depository, ancil_fields, mesh_id, &
                                twod_mesh_id, twod, ndata, time_axis )

    implicit none

    character(*),                   intent(in)    :: name
    type( field_collection_type ),  intent(inout) :: depository
    type( field_collection_type ),  intent(inout) :: ancil_fields
    integer(i_def),                 intent(in)    :: mesh_id
    integer(i_def),                 intent(in)    :: twod_mesh_id
    logical(l_def),       optional, intent(in)    :: twod
    integer(i_def),       optional, intent(in)    :: ndata
    type(time_axis_type), optional, intent(inout) :: time_axis

    ! Local variables
    type(field_type)          :: new_field
    integer(i_def)            :: ndat
    integer(i_def), parameter :: fs_order = 0

    ! Pointers
    type(function_space_type),       pointer :: wth_space => null()
    type(function_space_type),       pointer :: twod_space => null()
    procedure(read_interface),       pointer :: tmp_read_ptr => null()
    procedure(write_interface),      pointer :: tmp_write_ptr => null()
    class(field_type),               pointer :: fld_ptr => null()
    class(pure_abstract_field_type), pointer :: abs_fld_ptr => null()

    ! Set field ndata if argument is present, else leave as default value
    if (present(ndata)) then
      ndat = ndata
    else
      ndat = 1
    end if

    ! Set up function spaces for field initialisation
    wth_space   => function_space_collection%get_fs( mesh_id, fs_order, &
                                                     WTheta, ndat )
    twod_space => function_space_collection%get_fs( twod_mesh_id, fs_order, &
                                                    W3, ndat )

    ! If field does not yet exist, then create it
    if ( .not. depository%field_exists( name ) ) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating new field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)
      if (present(twod)) then
        call new_field%initialise( twod_space, name=trim(name) )
      else
        call new_field%initialise( wth_space, name=trim(name) )
      end if
      ! Add the new field to the field depository
      call depository%add_field(new_field)
    end if

    ! If field is time-varying, also create field storing raw data to be
    ! interpolated
    if (present(time_axis)) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating time axis field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)
      if (present(twod)) then
        call new_field%initialise( twod_space, name=trim(name) )
        call time_axis%add_field(new_field)
      else
        call new_field%initialise( wth_space, name=trim(name) )
        call time_axis%add_field(new_field)
      end if
    end if

    ! Get a field pointer from the depository
    fld_ptr => depository%get_field(name)
    ! Set up field read behaviour for 2D and 3D fields
    if (present(twod)) then
      tmp_write_ptr => write_field_single_face
    else
      tmp_write_ptr => write_field_face
    end if
    ! Set field write behaviour for target field
    call fld_ptr%set_write_behaviour(tmp_write_ptr)


    if (.not. present(time_axis)) then
      !Set up field read behaviour for 2D and 3D fields
      if (present(twod))then
        tmp_read_ptr => read_field_single_face
      else
        tmp_read_ptr => read_field_face
      end if

      ! Set field read behaviour for target field
      call fld_ptr%set_read_behaviour(tmp_read_ptr)
    end if

    ! Add the field pointer to the target field collection
    abs_fld_ptr => depository%get_field(name)
    call ancil_fields%add_reference_to_field(abs_fld_ptr)

    ! Nullify pointers
    nullify(wth_space)
    nullify(twod_space)
    nullify(tmp_read_ptr)
    nullify(tmp_write_ptr)
    nullify(fld_ptr)
    nullify(abs_fld_ptr)

  end subroutine setup_ancil_field

  !> @details Initialises a time_axis object for a group of time-varying ancils
  !> @param[in] time_id The name of the time axis (also the XIOS id of the time
  !>                    data)
  !> @param[out] time_axis The resulting time_axis object
  subroutine init_time_axis(time_id, time_axis)

    implicit none

    character(len=*),     intent(in)  :: time_id
    type(time_axis_type), intent(out) :: time_axis

    ! Local/derived variables for time axis initialisation
    real(dp_xios),  allocatable :: time_data_dpxios(:)
    real(r_def),    allocatable :: time_data(:)
    integer(i_def), allocatable :: axis_indices(:)

    ! Other local variables for XIOS interface
    integer(i_def) :: i

    ! Read time data from ancil file, then cast to r_def
    call read_time_data(time_id, time_data_dpxios)
    time_data = real(time_data_dpxios, kind=r_def)
    deallocate(time_data_dpxios)

    ! Set up axis index array
    allocate(axis_indices(size(time_data)))
    axis_indices = (/ (i, i=1,size(time_data),1) /)

    ! Initialise time axis
    call time_axis%initialise( time_data, axis_indices, time_id )

  end subroutine init_time_axis

end module init_ancils_mod
