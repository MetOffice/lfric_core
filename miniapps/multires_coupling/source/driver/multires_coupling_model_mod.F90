!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the multires_coupling miniapp model simulations.
!>
module multires_coupling_model_mod

  use assign_orography_field_mod, only : assign_orography_field
  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use checksum_alg_mod,           only : checksum_alg
  use check_configuration_mod,    only : get_required_stencil_depth, &
                                         check_any_shifted
  use driver_fem_mod,             only : init_fem, final_fem, &
                                         init_function_space_chains
  use driver_mesh_mod,            only : init_mesh
  use driver_io_mod,              only : init_io, final_io, &
                                         filelist_populator
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, i_native,          &
                                         PRECISION_REAL, r_second, &
                                         str_def, r_def, l_def
  use convert_to_upper_mod,       only : convert_to_upper
  use derived_config_mod,         only : set_derived_config
  use extrusion_mod,              only : extrusion_type, TWOD, &
                                         SHIFTED, DOUBLE_LEVEL
  use field_mod,                  only : field_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type
  use formulation_config_mod,     only : l_multigrid,              &
                                         moisture_formulation,     &
                                         moisture_formulation_dry, &
                                         use_physics,              &
                                         use_multires_coupling
  use geometric_constants_mod,    only : get_chi_inventory, get_panel_id_inventory
  use gungho_extrusion_mod,       only : create_extrusion
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use mesh_mod,                   only : mesh_type
  use model_clock_mod,            only : model_clock_type
  use multires_coupling_config_mod, &
                                  only : physics_mesh_name,           &
                                         multires_coupling_mesh_tags, &
                                         orography_mesh_name
  use gungho_model_data_mod,      only : model_data_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use init_altitude_mod,          only : init_altitude
  use io_config_mod,              only : write_conservation_diag, &
                                         write_minmax_tseries
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_DEBUG
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use mpi_mod,                    only : mpi_type
  use multigrid_config_mod,       only : chain_mesh_tags
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_final
  use runtime_constants_mod,      only : create_runtime_constants, &
                                         final_runtime_constants
  use semi_implicit_timestep_alg_mod, &
                                  only : semi_implicit_alg_init, &
                                         semi_implicit_alg_final
  use section_choice_config_mod,  only : radiation,         &
                                         radiation_socrates,&
                                         surface, surface_jules
  use setup_orography_alg_mod,    only : setup_orography_alg
  use timestepping_config_mod,    only : method,               &
                                         method_semi_implicit, &
                                         method_rk,            &
                                         spinup_period

#ifdef UM_PHYSICS
  use jules_control_init_mod,     only : jules_control_init
  use jules_physics_init_mod,     only : jules_physics_init
  use planet_constants_mod,       only : set_planet_constants
  use socrates_init_mod,          only : socrates_init
  use um_clock_init_mod,          only : um_clock_init
  use um_control_init_mod,        only : um_control_init
  use um_physics_init_mod,        only : um_physics_init
  use um_radaer_lut_init_mod,     only : um_radaer_lut_init
  use um_ukca_init_mod,           only : um_ukca_init
#endif

  implicit none

  private
  public initialise_infrastructure, &
         initialise_model,          &
         finalise_infrastructure,   &
         finalise_model

contains

  !> @brief Initialises the infrastructure and sets up constants used by the
  !>        model.
  !>
  !> @param[in]     program_name An identifier given to the model begin run
  !> @param[out]    model_clock  Time within the model
  !> @param[in,out] mpi          Communication object
  !>
  subroutine initialise_infrastructure( program_name, model_clock, calendar, mpi )

    use check_configuration_mod, only: get_required_stencil_depth

    implicit none

    character(*),           intent(in)    :: program_name
    type(model_clock_type), intent(inout) :: model_clock
    class(calendar_type),   intent(in)    :: calendar
    class(mpi_type),        intent(inout) :: mpi

    type(mesh_type), pointer :: mesh => null()
    type(mesh_type), pointer :: shifted_mesh => null()
    type(mesh_type), pointer :: double_level_mesh => null()
    type(mesh_type), pointer :: twod_mesh => null()
    type(mesh_type), pointer :: physics_mesh => null()
    type(mesh_type), pointer :: orography_twod_mesh => null()
    type(mesh_type), pointer :: orography_mesh => null()

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    type(field_type) :: surface_altitude
    class(extrusion_type), allocatable :: extrusion

    type(inventory_by_mesh_type), pointer :: chi_inventory => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()

    logical(l_def)                  :: mesh_already_exists
    integer(i_def)                  :: i, j, mesh_ctr, max_num_meshes
    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: shifted_mesh_names(:)
    character(str_def), allocatable :: double_level_mesh_names(:)
    character(str_def), allocatable :: tmp_mesh_names(:)
    character(str_def), allocatable :: extra_io_mesh_names(:)

#ifdef UM_PHYSICS
    integer(i_def) :: ncells

#endif

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Work out which meshes are required
    !-------------------------------------------------------------------------

    ! Gather together names of meshes
    max_num_meshes = 1
    if ( l_multigrid ) max_num_meshes = max_num_meshes + SIZE(chain_mesh_tags)
    if ( use_multires_coupling ) max_num_meshes = max_num_meshes + SIZE(multires_coupling_mesh_tags)
    ! Don't know the full number of meshes, so allocate maximum for now
    allocate(tmp_mesh_names(max_num_meshes))

    ! Prime extrusions -------------------------------------------------------
    ! Always have the prime mesh
    mesh_ctr = 1
    tmp_mesh_names(1) = prime_mesh_name
    ! Multigrid meshes
    if ( l_multigrid ) then
      do i = 1, SIZE(chain_mesh_tags)
        ! Only add mesh if it has not already been added
        mesh_already_exists = .false.
        do j = 1, mesh_ctr
          if ( tmp_mesh_names(j) == chain_mesh_tags(i) ) then
            mesh_already_exists = .true.
            exit
          end if
        end do
        if ( .not. mesh_already_exists ) then
          mesh_ctr = mesh_ctr + 1
          tmp_mesh_names(mesh_ctr) = chain_mesh_tags(i)
        end if
      end do
    end if
    ! Multires coupling meshes
    if ( use_multires_coupling ) then
      do i = 1, SIZE(multires_coupling_mesh_tags)
        ! Only add mesh if it has not already been added
        mesh_already_exists = .false.
        do j = 1, mesh_ctr
          if ( tmp_mesh_names(j) == multires_coupling_mesh_tags(i) ) then
            mesh_already_exists = .true.
            exit
          end if
        end do
        if ( .not. mesh_already_exists ) then
          mesh_ctr = mesh_ctr + 1
          tmp_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
        end if
      end do
    end if

    ! Transfer mesh names from temporary array to an array of appropriate size
    allocate(base_mesh_names(mesh_ctr))
    do i = 1, mesh_ctr
      base_mesh_names(i) = tmp_mesh_names(i)
    end do
    deallocate(tmp_mesh_names)

    ! Shifted meshes ---------------------------------------------------------
    if ( check_any_shifted() ) then
      allocate(tmp_mesh_names(SIZE(base_mesh_names)))

      mesh_ctr = 1
      tmp_mesh_names(1) = prime_mesh_name

      if ( use_multires_coupling ) then
        do i = 1, SIZE(multires_coupling_mesh_tags)
          ! Only add mesh if it has not already been added
          mesh_already_exists = .false.
          do j = 1, mesh_ctr
            if ( tmp_mesh_names(j) == multires_coupling_mesh_tags(i) ) then
              mesh_already_exists = .true.
              exit
            end if
          end do
          if ( .not. mesh_already_exists ) then
            mesh_ctr = mesh_ctr + 1
            tmp_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          end if
        end do
      end if

      ! Transfer mesh names from temporary array to an array of appropriate size
      allocate(shifted_mesh_names(mesh_ctr))
      do i = 1, mesh_ctr
        shifted_mesh_names(i) = tmp_mesh_names(i)
      end do
      deallocate(tmp_mesh_names)
    end if

    ! Double level meshes ------------------------------------------------------
    if ( check_any_shifted() ) then
      allocate(tmp_mesh_names(SIZE(base_mesh_names)))

      mesh_ctr = 1
      tmp_mesh_names(1) = prime_mesh_name

      if ( use_multires_coupling ) then
        do i = 1, SIZE(multires_coupling_mesh_tags)
          ! Only add mesh if it has not already been added
          mesh_already_exists = .false.
          do j = 1, mesh_ctr
            if ( tmp_mesh_names(j) == multires_coupling_mesh_tags(i) ) then
              mesh_already_exists = .true.
              exit
            end if
          end do
          if ( .not. mesh_already_exists ) then
            mesh_ctr = mesh_ctr + 1
            tmp_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          end if
        end do
      end if

      ! Transfer mesh names from temporary array to an array of appropriate size
      allocate(double_level_mesh_names(mesh_ctr))
      do i = 1, mesh_ctr
        double_level_mesh_names(i) = tmp_mesh_names(i)
      end do
      deallocate(tmp_mesh_names)
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    ! Generate prime mesh extrusion
    allocate( extrusion, source=create_extrusion() )

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),               &
                    base_mesh_names,                                        &
                    shifted_mesh_names      = shifted_mesh_names,           &
                    double_level_mesh_names = double_level_mesh_names,      &
                    required_stencil_depth  = get_required_stencil_depth(), &
                    input_extrusion         = extrusion )

    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )
    if ( l_multigrid ) then
      call init_function_space_chains( mesh_collection, chain_mesh_tags )
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    allocate(extra_io_mesh_names(SIZE(multires_coupling_mesh_tags)-1))
    mesh_ctr = 1
    do i = 1, SIZE(multires_coupling_mesh_tags)
      if (multires_coupling_mesh_tags(i) /= prime_mesh_name) then
        extra_io_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
        mesh_ctr = mesh_ctr + 1
      end if
    end do

    files_init_ptr => init_gungho_files
    call init_io( program_name,                      &
                  mpi%get_comm(),                    &
                  chi_inventory, panel_id_inventory, &
                  model_clock, calendar,             &
                  populate_filelist=files_init_ptr,  &
                  alt_mesh_names=extra_io_mesh_names )

    if ( use_multires_coupling ) then
      orography_mesh => mesh_collection%get_mesh(trim(orography_mesh_name))
      orography_twod_mesh => mesh_collection%get_mesh(orography_mesh, TWOD)
    else
      orography_mesh => mesh_collection%get_mesh(prime_mesh_name)
      orography_twod_mesh => mesh_collection%get_mesh(orography_mesh, TWOD)
    end if

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GAL)
    call init_altitude( orography_twod_mesh, surface_altitude )

    call setup_orography_alg( base_mesh_names,                &
                              orography_mesh%get_mesh_name(), &
                              chi_inventory,                  &
                              panel_id_inventory,             &
                              surface_altitude        )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_collection, chi_inventory, &
                                  panel_id_inventory, model_clock )

    physics_mesh => mesh_collection%get_mesh( physics_mesh_name )

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()

    ! Initialise UM to run in columns
    ncells = 1_i_def

    if ( use_physics ) then
      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call socrates_init()
      end if

      ! Initialisation of UM high-level variables
      call um_control_init( physics_mesh, ncells )

      ! Initialisation of UM clock
      call um_clock_init(clock)
      ! Initialisation of UM physics variables
      call um_physics_init(ncells)
      !Read all the radaer lut namelist files
      call um_radaer_lut_init()
      ! Initialisation of Jules high-level variables
      call jules_control_init()
      if (surface == surface_jules) then
        ! Initialisation of Jules physics variables
        call jules_physics_init()
      end if
      ! Initialisation of UKCA physics variables
      call um_ukca_init()
    end if
#endif

    nullify(physics_mesh, mesh, twod_mesh, shifted_mesh, double_level_mesh, &
            chi_inventory, panel_id_inventory, files_init_ptr)
    deallocate(base_mesh_names)
    if (allocated(shifted_mesh_names)) deallocate(shifted_mesh_names)
    if (allocated(double_level_mesh_names)) deallocate(double_level_mesh_names)

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the multires_coupling application
  !>
  !> @param[in] dynamics_mesh           The dynamics mesh
  !> @param[in,out] dynamics_model_data The dynamics data set
  !> @param[in] physics_mesh            The physics mesh
  !> @param[in,out] physics_model_data  The physics data set
  !>
  subroutine initialise_model( dynamics_mesh,       &
                               dynamics_model_data, &
                               physics_mesh,        &
                               physics_model_data )
    implicit none

    type(mesh_type),       pointer, intent(in)    :: dynamics_mesh
    type(mesh_type),       pointer, intent(in)    :: physics_mesh
    type(model_data_type), target,  intent(inout) :: dynamics_model_data
    type(model_data_type), target,  intent(inout) :: physics_model_data

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_type),            pointer :: dynamics_mr(:) => null()
    type(field_type),            pointer :: physics_mr(:) => null()

    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_exner => null()
    type(field_type), pointer :: physics_theta => null()
    type(field_type), pointer :: physics_u => null()
    type(field_type), pointer :: physics_rho => null()
    type(field_type), pointer :: physics_exner => null()

    ! Get pointers to field collections for use downstream
    dynamics_mr => dynamics_model_data%mr
    physics_mr  => physics_model_data%mr

    prognostic_fields => dynamics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', dynamics_theta)
    call prognostic_fields%get_field('u', dynamics_u)
    call prognostic_fields%get_field('rho', dynamics_rho)
    call prognostic_fields%get_field('exner', dynamics_exner)

    prognostic_fields => physics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', physics_theta)
    call prognostic_fields%get_field('u', physics_u)
    call prognostic_fields%get_field('rho', physics_rho)
    call prognostic_fields%get_field('exner', physics_exner)

    if (write_minmax_tseries) then
      call minmax_tseries_init('u')
      call minmax_tseries(dynamics_u, 'u', dynamics_mesh)
    end if

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        ! Initialise and output initial conditions for first timestep
        call semi_implicit_alg_init(dynamics_mesh, dynamics_u, dynamics_rho, &
                                    dynamics_theta, dynamics_exner, dynamics_mr)
        if ( write_conservation_diag ) then
         call conservation_algorithm( dynamics_rho,     &
                                      dynamics_u,       &
                                      dynamics_theta,   &
                                      dynamics_mr,      &
                                      dynamics_exner )
         if ( moisture_formulation /= moisture_formulation_dry ) &
           call moisture_conservation_alg( dynamics_rho,         &
                                           dynamics_mr,          &
                                           'Before timestep' )
        end if
      case( method_rk )             ! RK
        ! Initialise and output initial conditions for first timestep
        call rk_alg_init(dynamics_mesh, dynamics_u, dynamics_rho, &
                         dynamics_theta, dynamics_exner)
        if ( write_conservation_diag ) then
         call conservation_algorithm( dynamics_rho,     &
                                      dynamics_u,       &
                                      dynamics_theta,   &
                                      dynamics_mr,      &
                                      dynamics_exner )
         if ( moisture_formulation /= moisture_formulation_dry ) &
           call moisture_conservation_alg( dynamics_rho,         &
                                           dynamics_mr,          &
                                           'Before timestep' )
        end if
      case default
        call log_event("Multires Coupling: Incorrect time stepping option chosen, "// &
                        "stopping program! ", LOG_LEVEL_ERROR)
    end select

  end subroutine initialise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalises infrastructure and constants used by the model.
  !>
  subroutine finalise_infrastructure()

    implicit none

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise I/O
    !-------------------------------------------------------------------------

    call final_io()

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------

    call final_fem()

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the Multires Coupling application
  !>
  !> @param[in,out] dynamics_model_data The working data set for the model run
  !> @param[in,out] physics_model_data  The working data set for the model run
  !> @param[in] program_name An identifier given to the model begin run
  !>
  subroutine finalise_model( dynamics_model_data, &
                             physics_model_data,  &
                             program_name )
    implicit none


    type(model_data_type), target, intent(inout) :: dynamics_model_data
    type(model_data_type), target, intent(inout) :: physics_model_data
    character(*),                  intent(in)    :: program_name

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_type),            pointer :: dynamics_mr(:) => null()
    type(field_type),            pointer :: physics_mr(:) => null()

    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_exner => null()
    type(field_type), pointer :: physics_theta => null()
    type(field_type), pointer :: physics_u => null()
    type(field_type), pointer :: physics_rho => null()
    type(field_type), pointer :: physics_exner => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to field collections for use downstream
    dynamics_mr => dynamics_model_data%mr
    physics_mr => physics_model_data%mr

    prognostic_fields => dynamics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', dynamics_theta)
    call prognostic_fields%get_field('u', dynamics_u)
    call prognostic_fields%get_field('rho', dynamics_rho)
    call prognostic_fields%get_field('exner', dynamics_exner)

    prognostic_fields => physics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', physics_theta)
    call prognostic_fields%get_field('u', physics_u)
    call prognostic_fields%get_field('rho', physics_rho)
    call prognostic_fields%get_field('exner', physics_exner)

    ! Log fields
    call dynamics_rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call dynamics_theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call dynamics_exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call dynamics_u%log_field(     LOG_LEVEL_DEBUG, 'u' )
    call physics_rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call physics_theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call physics_exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call physics_u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if ( moisture_formulation /= moisture_formulation_dry ) then
      call checksum_alg(program_name,             &
                        dynamics_rho, 'rho',      &
                        dynamics_theta, 'theta',  &
                        dynamics_u, 'u',          &
                        field_bundle=dynamics_mr, &
                        bundle_name='mr')
    else
      call checksum_alg(program_name,             &
                        dynamics_rho, 'rho',      &
                        dynamics_theta, 'theta',  &
                        dynamics_u, 'u')
    end if

    if (write_minmax_tseries) call minmax_tseries_final()

    if ( method == method_semi_implicit ) call semi_implicit_alg_final()
    if ( method == method_rk )            call rk_alg_final()

  end subroutine finalise_model

end module multires_coupling_model_mod
