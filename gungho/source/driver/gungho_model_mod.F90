!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the gungho model simulations.
!>
module gungho_model_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use checksum_alg_mod,           only : checksum_alg
  use driver_fem_mod,             only : init_fem, final_fem, &
                                         init_function_space_chains
  use driver_io_mod,              only : init_io, final_io,  &
                                         filelist_populator
  use driver_mesh_mod,            only : init_mesh
  use check_configuration_mod,    only : get_required_stencil_depth, &
                                         check_any_shifted
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, i_native, r_def, l_def, &
                                         PRECISION_REAL, r_second, str_def
  use convert_to_upper_mod,       only : convert_to_upper
  use derived_config_mod,         only : set_derived_config
  use extrusion_mod,              only : extrusion_type, TWOD, &
                                         SHIFTED, DOUBLE_LEVEL
  use field_mod,                  only : field_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type
  use field_spec_mod,             only : field_spec_type, processor_type
  use create_gungho_prognostics_mod, &
                                  only : enable_gungho_prognostics, &
                                         enable_checkpointing
  use boundaries_config_mod,      only : limited_area
  use create_lbcs_mod,            only : enable_lbc_fields
  use create_physics_prognostics_mod, &
                                  only : process_physics_prognostics
  use multigrid_config_mod,       only : chain_mesh_tags
  use multires_coupling_config_mod, &
                                  only : multires_coupling_mesh_tags, &
                                         orography_mesh_name
  use formulation_config_mod,     only : l_multigrid,              &
                                         moisture_formulation,     &
                                         moisture_formulation_dry, &
                                         use_physics,              &
                                         use_multires_coupling
  use geometric_constants_mod,    only : get_chi_inventory, get_panel_id_inventory
  use gungho_extrusion_mod,       only : create_extrusion
  use gungho_model_data_mod,      only : model_data_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use gungho_transport_control_alg_mod, &
                                  only : gungho_transport_control_alg_final
  use init_altitude_mod,          only : init_altitude
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use io_config_mod,              only : use_xios_io,             &
                                         write_conservation_diag, &
                                         write_dump,              &
                                         write_minmax_tseries
  use lfric_xios_context_mod,     only : lfric_xios_context_type
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_ALWAYS
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use clock_mod,                  only : clock_type
  use model_clock_mod,            only : model_clock_type
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use mpi_mod,                    only : mpi_type
  use mr_indices_mod,             only : nummr
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
  use time_config_mod,            only : timestep_end, timestep_start
  use timestepping_config_mod,    only : dt,                     &
                                         method,                 &
                                         method_semi_implicit,   &
                                         method_rk,              &
                                         method_no_timestepping, &
                                         spinup_period
  use derived_config_mod,         only : l_esm_couple
#ifdef COUPLED
  use coupler_mod,                only : cpl_define, cpl_fields
#endif

#ifdef UM_PHYSICS
  use gas_calc_all_mod,            only : gas_calc_all
  use jules_control_init_mod,      only : jules_control_init
  use jules_physics_init_mod,      only : jules_physics_init
  use planet_constants_mod,        only : set_planet_constants
  use socrates_init_mod,           only : socrates_init
  use illuminate_alg_mod,          only : illuminate_alg
  use um_clock_init_mod,           only : um_clock_init
  use um_control_init_mod,         only : um_control_init
  use um_domain_init_mod,          only : um_domain_init
  use um_domain_init_mod,          only : um_domain_init
  use um_sizes_init_mod,           only : um_sizes_init
  use um_physics_init_mod,         only : um_physics_init
  use um_radaer_lut_init_mod,      only : um_radaer_lut_init
  use um_ukca_init_mod,            only : um_ukca_init
  use random_seed_gen_alg_mod,     only : random_seed_gen_alg
  use stochastic_physics_config_mod, only : use_spt, use_skeb
#endif

  implicit none

  private

  logical(l_def) :: use_moisture

  !> @brief Processor class for enabling checkpointing fields
  type, extends(processor_type) :: enabler_type
  contains
    private
    ! main interface
    procedure, public :: init => enabler_init
    procedure, public :: apply => enabler_apply

    ! destructor - here to avoid gnu compiler bug
    final :: enabler_dtor
  end type enabler_type

  public initialise_infrastructure, &
         initialise_model,          &
         finalise_infrastructure,   &
         finalise_model
contains
  !> @brief  Initialise processor object for enabling checkpointing
  !> @param[in]   self      Enabler object
  !> @param[in]   clock     Model clock
  subroutine enabler_init(self, clock)
    implicit none
    class(enabler_type),       intent(inout) :: self
    class(clock_type), target, intent(in)    :: clock

    call self%set_clock(clock)
  end subroutine enabler_init

  !> @brief     Enabler's apply method
  !> @details   Enable a field given by field specifier for checkpointing.
  !> @param[in]   self      Enabler object
  !> @param[in]   spec      Field specifier
  subroutine enabler_apply(self, spec)
    implicit none
    class(enabler_type),   intent(in) :: self
    type(field_spec_type), intent(in) :: spec

    if (spec%ckp) then
      call enable_checkpointing(spec%name)
    end if
  end subroutine enabler_apply

  !> @brief Destructor of enabler object
  !> @param[inout] self  Enabler object
  subroutine enabler_dtor(self)
    type(enabler_type), intent(inout) :: self
    ! empty
  end subroutine enabler_dtor

  !> @brief Enable active state fields for checkpointing; sync xios axis dimensions
  !> @param[in] clock        The clock providing access to time information
  subroutine before_context_close(clock)
    use multidata_field_dimensions_mod, only: sync_multidata_field_dimensions
    use time_dimensions_mod, only: sync_time_dimensions

    implicit none
    class(clock_type), intent(in) :: clock

    type(enabler_type) :: enabler

    call enable_gungho_prognostics()
    if (limited_area) call enable_lbc_fields()
    if (use_physics) then
      call enabler%init(clock)
      call process_physics_prognostics(enabler)
      call sync_multidata_field_dimensions()
      call sync_time_dimensions()
    end if
  end subroutine before_context_close

  !> @brief Initialisations that depend only on loaded namelists.
  !> @details To be called before the IO context opens so that the
  !> multidata dimensions are completely defined when the XIOS
  !> axis dimensions are being synched.
  !> @param[in] model_clock     The clock providing access to time information
  subroutine basic_initialisations(model_clock)
    implicit none

    class(model_clock_type), intent(inout) :: model_clock

#ifdef UM_PHYSICS
    integer(i_def) :: ncells

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
      call um_control_init()
      call um_sizes_init(ncells)
      ! Initialisation of UM clock
      call um_clock_init(model_clock)

      ! Initialisation of UM physics variables
      call um_physics_init()
      ! Read all the radaer lut namelist files
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
  end subroutine basic_initialisations

  !> @brief Initialises the infrastructure and sets up constants used by the
  !>        model.
  !>
  !> @param [in,out] model_data   The working data set for the model run
  !> @param [out]    model_clock  Time within the model
  !> @param [in]     mpi          Communication object
  !>
  subroutine initialise_infrastructure( model_data, &
                                        model_clock, calendar, mpi )

    use logging_config_mod, only: key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    ! @todo I think the communication object aught to work as intent(in) but
    !       there seems to be an issue with Intel 19 which means (inout) is
    !       needed.
    !
    class(mpi_type), intent(inout) :: mpi

    type(model_data_type),   intent(inout) :: model_data
    class(model_clock_type), intent(inout) :: model_clock
    class(calendar_type),    intent(in)    :: calendar

    character(len=*), parameter :: io_context_name = "gungho_atm"

    type(mesh_type), pointer :: mesh => null()
    type(mesh_type), pointer :: shifted_mesh => null()
    type(mesh_type), pointer :: double_level_mesh => null()
    type(mesh_type), pointer :: twod_mesh => null()
    type(mesh_type), pointer :: orography_twod_mesh => null()
    type(mesh_type), pointer :: orography_mesh => null()

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    type(field_type) :: surface_altitude

    class(extrusion_type), allocatable :: extrusion

    type(field_type),     pointer :: chi(:) => null()

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

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)

    !-------------------------------------------------------------------------
    ! initialize coupling
    !-------------------------------------------------------------------------
!> @todo this must be done in infrastructure for now (before XIOS context
!>       initialization). With XIOS 3 it will be possible to move it outside
!>       infrastructure and remove change in gungho_prognostics_mod.f90
#ifdef COUPLED
    if( l_esm_couple ) then
       call log_event("Initialising coupler", LOG_LEVEL_INFO)
       ! Add fields used in coupling
       call cpl_fields( mesh, twod_mesh, model_data%depository, &
                        model_data%prognostic_fields )
       ! Define coupling interface
       call model_data%cpl_snd%initialise(name="cpl_snd")
       call model_data%cpl_rcv%initialise(name="cpl_rcv")
       call cpl_define( twod_mesh, chi, model_data%depository, &
                        model_data%cpl_snd, model_data%cpl_rcv )

    endif
#endif

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    call basic_initialisations( model_clock )

    call log_event("Initialising I/O context", LOG_LEVEL_INFO)

    files_init_ptr => init_gungho_files

    if ( use_multires_coupling ) then
      ! Compose list of meshes for I/O
      ! This is list of namelist mesh names minus the prime mesh name
      allocate(extra_io_mesh_names(SIZE(multires_coupling_mesh_tags)-1))
      mesh_ctr = 1
      do i = 1, SIZE(multires_coupling_mesh_tags)
        if (multires_coupling_mesh_tags(i) /= prime_mesh_name) then
          extra_io_mesh_names(mesh_ctr) = multires_coupling_mesh_tags(i)
          mesh_ctr = mesh_ctr + 1
        end if
      end do

      call init_io( io_context_name, mpi%get_comm(),   &
                    chi_inventory, panel_id_inventory, &
                    model_clock, calendar,             &
                    populate_filelist=files_init_ptr,  &
                    alt_mesh_names=extra_io_mesh_names,&
                    before_close=before_context_close )

    else
      call init_io( io_context_name, mpi%get_comm(),   &
                    chi_inventory, panel_id_inventory, &
                    model_clock, calendar,             &
                    populate_filelist=files_init_ptr,  &
                    before_close=before_context_close )
    end if

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
    call create_runtime_constants( mesh_collection,    &
                                   chi_inventory,      &
                                   panel_id_inventory, &
                                   model_clock )
#ifdef UM_PHYSICS
    if ( use_physics ) then

      ! Initialise time-varying trace gases
      call gas_calc_all()

      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call illuminate_alg( model_data%radiation_fields, &
                             model_clock%get_step(),      &
                             model_clock%get_seconds_per_step())
      end if
      ! Initialisation of UM variables related to the mesh
      call um_domain_init(mesh)

      ! Initialise random seed for stochastic physics
      if ( use_spt .or. use_skeb ) then
        call random_seed_gen_alg()
      end if

    end if
#endif
    nullify(mesh, twod_mesh, shifted_mesh, double_level_mesh, chi, &
            chi_inventory, panel_id_inventory, files_init_ptr,     &
            orography_mesh, orography_twod_mesh)
    deallocate(base_mesh_names)
    if (allocated(shifted_mesh_names)) deallocate(shifted_mesh_names)
    if (allocated(double_level_mesh_names)) deallocate(double_level_mesh_names)

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the gungho application
  !>
  !> @param[in] mesh  The primary mesh
  !> @param[in,out] model_data The working data set for the model run
  !>
  subroutine initialise_model( mesh, model_data )
    implicit none

    type( mesh_type ),       intent(in),    pointer :: mesh
    type( model_data_type ), intent(inout), target  :: model_data

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    use_moisture = ( moisture_formulation /= moisture_formulation_dry )

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    mr => model_data%mr

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    if (write_minmax_tseries) then
      call minmax_tseries_init('u')
      call minmax_tseries(u, 'u', mesh)
    end if

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        ! Initialise and output initial conditions for first timestep
        call semi_implicit_alg_init(mesh, u, rho, theta, exner, mr)

        if ( write_conservation_diag ) then
         call conservation_algorithm( rho,              &
                                      u,                &
                                      theta,            &
                                      mr,               &
                                      exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_rk )             ! RK
        ! Initialise and output initial conditions for first timestep
        call rk_alg_init(mesh, u, rho, theta, exner)
        if ( write_conservation_diag ) then
         call conservation_algorithm( rho,              &
                                      u,                &
                                      theta,            &
                                      mr,               &
                                      exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_no_timestepping )
        write( log_scratch_space, &
           '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                      ' Prognostic fields not evolved'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

      case default
        call log_event("Gungho: Incorrect time stepping option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
    end select

  end subroutine initialise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalises infrastructure and constants used by the model.
  !>
  subroutine finalise_infrastructure()

    implicit none

    !-------------------------------------------------------------------------
    ! Finalise I/O
    !-------------------------------------------------------------------------

    call final_io()

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------
    call final_fem()

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the gungho application
  !>
  !> @param[in,out] model_data The working data set for the model run
  !> @param[in] program_name   An identifier given to the model begin run
  !>
  subroutine finalise_model( model_data, &
                             program_name )

    implicit none

    type( model_data_type ), target,  intent(inout) :: model_data
    character(*),                     intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: fd_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr
    fd_fields => model_data%fd_fields

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    ! Log fields
    call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if (use_moisture) then
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u', &
         field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u')
    end if

    if (write_minmax_tseries) call minmax_tseries_final()

    if ( method == method_semi_implicit ) call semi_implicit_alg_final()
    if ( method == method_rk )            call rk_alg_final()
    call gungho_transport_control_alg_final()

  end subroutine finalise_model

end module gungho_model_mod
