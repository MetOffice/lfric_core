!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the transport miniapp.
!>
module transport_driver_mod

  use base_mesh_config_mod,             only: prime_mesh_name
  use calendar_mod,                     only: calendar_type
  use checksum_alg_mod,                 only: checksum_alg
  use check_configuration_mod,          only: get_required_stencil_depth
  use configuration_mod,                only: final_configuration
  use constants_mod,                    only: i_def, l_def, &
                                              r_def, r_second, str_def
  use driver_fem_mod,                   only: init_fem
  use driver_io_mod,                    only: init_io, final_io
  use driver_mesh_mod,                  only: init_mesh
  use driver_time_mod,                  only: init_time, get_calendar
  use derived_config_mod,               only: set_derived_config
  use diagnostics_io_mod,               only: write_scalar_diagnostic, &
                                              write_vector_diagnostic
  use divergence_alg_mod,               only: divergence_alg
  use field_mod,                        only: field_type
  use formulation_config_mod,           only: use_multires_coupling
  use geometric_constants_mod,          only: get_chi_inventory, &
                                              get_panel_id_inventory
  use io_context_mod,                   only: io_context_type
  use io_config_mod,                    only: diagnostic_frequency, &
                                              nodal_output_on_w3,   &
                                              subroutine_timers,    &
                                              write_diag,           &
                                              use_xios_io
  use inventory_by_mesh_mod,            only: inventory_by_mesh_type
  use local_mesh_mod,                   only: local_mesh_type
  use log_mod,                          only: log_event,         &
                                              log_scratch_space, &
                                              LOG_LEVEL_ALWAYS,  &
                                              LOG_LEVEL_INFO,    &
                                              LOG_LEVEL_TRACE
  use mass_conservation_alg_mod,        only: mass_conservation
  use mesh_mod,                         only: mesh_type
  use mesh_collection_mod,              only: mesh_collection
  use model_clock_mod,                  only: model_clock_type
  use mpi_mod,                          only: mpi_type
  use mr_indices_mod,                   only: nummr
  use multires_coupling_config_mod,     only: aerosol_mesh_name,               &
                                              multires_coupling_mesh_tags
  use runtime_constants_mod,            only: create_runtime_constants
  use step_calendar_mod,                only: step_calendar_type
  use timer_mod,                        only: timer
  use timestepping_config_mod,          only: dt
  use transport_init_fields_alg_mod,    only: transport_init_fields_alg
  use transport_control_alg_mod,        only: transport_prerun_setup,          &
                                              transport_init, transport_step,  &
                                              transport_final, use_w2_vector,  &
                                              use_aerosols
  use transport_runtime_collection_mod, only: init_transport_runtime_collection, &
                                              transport_runtime_collection_final

  implicit none

  private

  public :: initialise_transport, step_transport, finalise_transport

  ! Prognostic fields
  type(field_type) :: wind
  type(field_type) :: density
  type(field_type) :: theta
  type(field_type) :: tracer_con
  type(field_type) :: tracer_adv
  type(field_type) :: constant
  type(field_type) :: mr(nummr)
  type(field_type) :: w2_vector
  type(field_type) :: divergence
  type(field_type) :: w3_aerosol
  type(field_type) :: wt_aerosol
  type(field_type) :: aerosol_wind

  ! Number of moisutre species to transport
  integer(kind=i_def) :: nummr_to_transport

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in,out] mpi          The structure that holds comms details
  !> @param [in]     program_name An identifier given to the model being run
  subroutine initialise_transport( mpi, model_clock, program_name, calendar )

    implicit none

    class(mpi_type),         intent(inout) :: mpi
    class(model_clock_type), intent(inout) :: model_clock
    character(*),            intent(in)    :: program_name
    class(calendar_type),    intent(in)    :: calendar

    character(len=*), parameter :: xios_ctx  = "transport"

    integer(kind=i_def)                   :: num_base_meshes
    integer(kind=i_def),      allocatable :: local_mesh_ids(:)
    type(local_mesh_type),        pointer :: local_mesh => null()
    type(mesh_type),              pointer :: mesh => null()
    type(mesh_type),              pointer :: aerosol_mesh => null()
    type(inventory_by_mesh_type), pointer :: chi_inventory => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()
    character(str_def),       allocatable :: base_mesh_names(:)
    character(str_def),       allocatable :: shifted_mesh_names(:)
    character(str_def),       allocatable :: double_level_mesh_names(:)
    character(str_def),       allocatable :: extra_io_mesh_names(:)
    logical(l_def)                        :: create_rdef_div_operators

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Work out which meshes are required
    !-------------------------------------------------------------------------

    if ( use_multires_coupling ) then
      num_base_meshes = 2
    else
      num_base_meshes = 1
    end if

    ! Always have just one base mesh with shifted version
    allocate(base_mesh_names(num_base_meshes))
    allocate(shifted_mesh_names(num_base_meshes))
    allocate(double_level_mesh_names(num_base_meshes))
    base_mesh_names(1) = prime_mesh_name
    shifted_mesh_names(1) = prime_mesh_name
    double_level_mesh_names(1) = prime_mesh_name

    if ( use_multires_coupling ) then
      base_mesh_names(2) = aerosol_mesh_name
      shifted_mesh_names(2) = aerosol_mesh_name
      double_level_mesh_names(2) = aerosol_mesh_name
    end if

    !-------------------------------------------------------------------------
    ! Initialise mesh and FEM infrastructure
    !-------------------------------------------------------------------------

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),           &
                    base_mesh_names,                                    &
                    shifted_mesh_names=shifted_mesh_names,              &
                    double_level_mesh_names=double_level_mesh_names,    &
                    required_stencil_depth=get_required_stencil_depth() )

    ! FEM initialisation
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Create runtime_constants object.
    create_rdef_div_operators = .true.
    call create_runtime_constants( mesh_collection,    &
                                   chi_inventory,      &
                                   panel_id_inventory, &
                                   model_clock,        &
                                   create_rdef_div_operators )

    ! Set up transport runtime collection type
    ! Transport on only one horizontal local mesh
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    local_mesh => mesh%get_local_mesh()

    if ( use_multires_coupling ) then
      allocate(local_mesh_ids(2))
      local_mesh_ids(1) = local_mesh%get_id()
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      local_mesh => aerosol_mesh%get_local_mesh()
      local_mesh_ids(2) = local_mesh%get_id()
    else
      allocate(local_mesh_ids(1))
      local_mesh_ids(1) = local_mesh%get_id()
      aerosol_mesh => mesh_collection%get_mesh(prime_mesh_name)
    end if

    call init_transport_runtime_collection(local_mesh_ids)

    ! Set transport metadata for primal mesh
    call transport_prerun_setup( num_base_meshes )

    ! Initialise prognostic variables
    call transport_init_fields_alg( mesh, wind, density, theta, &
                                    tracer_con, tracer_adv,     &
                                    constant, mr, w2_vector,    &
                                    aerosol_mesh, aerosol_wind, &
                                    w3_aerosol,  wt_aerosol,    &
                                    divergence )

    ! Initialise all transport-only control algorithm
    call transport_init( density, theta, tracer_con, tracer_adv,          &
                         constant, mr, w2_vector, w3_aerosol, wt_aerosol  )

    nummr_to_transport = 1_i_def

    ! I/O initialisation
    if (use_multires_coupling) then
      allocate(extra_io_mesh_names(1))
      extra_io_mesh_names(1) = aerosol_mesh%get_mesh_name()
      call init_io( xios_ctx,           &
                    mpi%get_comm(),     &
                    chi_inventory,      &
                    panel_id_inventory, &
                    model_clock,        &
                    get_calendar(),     &
                    alt_mesh_names=extra_io_mesh_names )
    else
      call init_io( xios_ctx,           &
                    mpi%get_comm(),     &
                    chi_inventory,      &
                    panel_id_inventory, &
                    model_clock,        &
                    get_calendar() )
    end if

    ! Output initial conditions
    if (model_clock%is_initialisation() .and. write_diag) then

      call write_vector_diagnostic( 'u', wind, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'rho', density, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'theta', theta, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_con', tracer_con, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_adv', tracer_adv, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'constant', constant, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'm_v', mr(1), model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'divergence', divergence, model_clock, &
                                    mesh, nodal_output_on_w3 )
      if (use_w2_vector) then
        call write_vector_diagnostic( 'w2_vector', w2_vector, model_clock, &
                                      mesh, nodal_output_on_w3 )
      end if
      if (use_aerosols) then
        call write_vector_diagnostic( 'aerosol_wind', aerosol_wind, model_clock, &
                                      aerosol_mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'w3_aerosol', w3_aerosol, model_clock, &
                                      aerosol_mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'wt_aerosol', wt_aerosol, model_clock, &
                                      aerosol_mesh, nodal_output_on_w3 )
      end if
    end if

    deallocate(base_mesh_names)
    deallocate(shifted_mesh_names)
    deallocate(double_level_mesh_names)
    if (allocated(extra_io_mesh_names)) deallocate(extra_io_mesh_names)
    deallocate(local_mesh_ids)
    nullify(chi_inventory, panel_id_inventory, mesh, local_mesh, aerosol_mesh)

  end subroutine initialise_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs a time step.
  !>
  subroutine step_transport( model_clock )

    implicit none

    class(model_clock_type), intent(in) :: model_clock

    type(mesh_type), pointer :: mesh => null()
    type(mesh_type), pointer :: aerosol_mesh => null()

    call log_event( 'Miniapp will run with default precision set as:', LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( '        r_def kind = '//log_scratch_space , LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( '        i_def kind = '//log_scratch_space , LOG_LEVEL_INFO )

    call mass_conservation( model_clock%get_step(), density, mr )
    call density%log_minmax( LOG_LEVEL_INFO, 'rho' )
    call theta%log_minmax( LOG_LEVEL_INFO, 'theta' )
    call tracer_con%log_minmax( LOG_LEVEL_INFO, 'tracer_con' )
    call tracer_adv%log_minmax( LOG_LEVEL_INFO, 'tracer_adv' )
    call constant%log_minmax( LOG_LEVEL_INFO, 'constant' )
    call mr(1)%log_minmax( LOG_LEVEL_INFO, 'm_v' )

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    if (use_multires_coupling) then
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
    else
      aerosol_mesh => mesh_collection%get_mesh(prime_mesh_name)
    end if

    write(log_scratch_space, '("/", A, "\ ")') repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) &
      'Start of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if ( subroutine_timers ) call timer( 'transport step' )

    call transport_step( model_clock,                          &
                         wind, density, theta, tracer_con,     &
                         tracer_adv, constant, mr, w2_vector,  &
                         w3_aerosol, wt_aerosol, aerosol_wind, &
                         nummr_to_transport )

    if ( subroutine_timers ) call timer( 'transport step' )

    ! Write out conservation diagnostics
    call mass_conservation( model_clock%get_step(), density, mr )
    call density%log_minmax( LOG_LEVEL_INFO, 'rho' )
    call theta%log_minmax( LOG_LEVEL_INFO, 'theta' )
    call tracer_con%log_minmax( LOG_LEVEL_INFO, 'tracer_con' )
    call tracer_adv%log_minmax( LOG_LEVEL_INFO, 'tracer_adv' )
    call constant%log_minmax( LOG_LEVEL_INFO, 'constant' )
    call mr(1)%log_minmax( LOG_LEVEL_INFO, 'm_v' )
    if (use_aerosols) then
      call w3_aerosol%log_minmax( LOG_LEVEL_INFO, 'w3_aerosol' )
      call wt_aerosol%log_minmax( LOG_LEVEL_INFO, 'wt_aerosol' )
    end if

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write(log_scratch_space, '("\", A, "/ ")') repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Output wind and density values.
    if ( (mod( model_clock%get_step(), diagnostic_frequency ) == 0) &
         .and. write_diag ) then

      ! Compute divergence
      call divergence_alg( divergence, wind )

      call write_vector_diagnostic( 'u', wind,                &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'rho', density,           &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'theta', theta,           &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_con', tracer_con, &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_adv', tracer_adv, &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'constant', constant,     &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'm_v', mr(1),             &
                                    model_clock, mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'divergence', divergence, &
                                    model_clock, mesh, nodal_output_on_w3 )
      if (use_w2_vector) then
        call write_vector_diagnostic( 'w2_vector', w2_vector,   &
                                      model_clock, mesh, nodal_output_on_w3 )
      end if
      if (use_aerosols) then
        call write_vector_diagnostic( 'aerosol_wind', aerosol_wind, model_clock, &
                                      aerosol_mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'w3_aerosol', w3_aerosol,   &
                                      model_clock, aerosol_mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'wt_aerosol', wt_aerosol,   &
                                      model_clock, aerosol_mesh, nodal_output_on_w3 )
      end if
    end if

    nullify(mesh)

  end subroutine step_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_transport( program_name )

    implicit none

    character(*), intent(in) :: program_name

    call transport_final( density, theta, tracer_con, tracer_adv, &
                          constant, mr, w2_vector, w3_aerosol, wt_aerosol )

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    ! Write checksums to file
    call checksum_alg( program_name, density, 'rho',  wind, 'u',  &
                       theta, 'theta', tracer_adv, 'tracer',      &
                       field_bundle=mr, bundle_name='mr' )

    call final_io()

    call transport_runtime_collection_final()

  end subroutine finalise_transport

end module transport_driver_mod
