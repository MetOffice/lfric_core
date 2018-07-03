!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use checksum_alg_mod,           only : checksum_alg
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, imdi
  use yz_bip_cosmic_alg_mod,      only : yz_bip_cosmic_step
  use cosmic_threed_alg_mod,      only : cosmic_threed_transport_step
  use cusph_cosmic_transport_alg_mod, &
                                  only : cusph_cosmic_transport_init, &
                                         cusph_cosmic_transport_step
  use derived_config_mod,         only : set_derived_config
  use diagnostic_alg_mod,         only : divergence_diagnostic_alg, &
                                         density_diagnostic_alg,    &
                                         hydbal_diagnostic_alg
  use yaxt,                       only : xt_initialize, xt_finalize
  use field_mod,                  only : field_type
  use formulation_config_mod,     only : transport_only, &
                                         use_moisture,   &
                                         use_physics
  use function_space_collection_mod, &
                                  only : function_space_collection
  use field_collection_mod,       only : field_collection_type
  use global_mesh_collection_mod, only : global_mesh_collection, &
                                         global_mesh_collection_type
  use gungho_mod,                 only : load_configuration, final_configuration
  use init_fem_mod,               only : init_fem
  use init_gungho_mod,            only : init_gungho
  use init_mesh_mod,              only : init_mesh
  use init_physics_mod,           only : init_physics
  use io_mod,                     only : output_nodal, &
                                         output_xios_nodal, &
                                         xios_domain_init
  use iter_timestep_alg_mod,      only : iter_alg_init, &
                                         iter_alg_step, &
                                         iter_alg_final
  use log_mod,                    only : log_event,         &
                                         log_set_level,     &
                                         log_scratch_space, &
                                         initialise_logging, &
                                         finalise_logging, &
                                         LOG_LEVEL_ERROR,   &
                                         LOG_LEVEL_INFO,    &
                                         LOG_LEVEL_DEBUG,   &
                                         LOG_LEVEL_TRACE
  use mesh_collection_mod,        only : mesh_collection
  use mod_wait
  use minmax_tseries_mod,         only : minmax_tseries, &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mr_indices_mod,             only : imr_v, imr_c, imr_r, imr_nc, &
                                         imr_nr, nummr
  use output_config_mod,          only : diagnostic_frequency, &
                                         subroutine_timers, &
                                         write_minmax_tseries, &
                                         subroutine_counters, &
                                         write_nodal_output, &
                                         write_xios_output

  use restart_config_mod,         only : restart_filename => filename
  use restart_control_mod,        only : restart_type
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_step, &
                                         rk_alg_final
  use rk_transport_mod,           only : rk_transport_init, &
                                         rk_transport_step, &
                                         rk_transport_final
  use runge_kutta_init_mod,       only : runge_kutta_init, runge_kutta_final
  use runtime_constants_mod,      only : final_runtime_constants
  use timer_mod,                  only : timer, output_timer
  use timestepping_config_mod,    only : method, dt, &
                                         timestepping_method_semi_implicit, &
                                         timestepping_method_rk
  use transport_config_mod,       only : scheme, &
                                         transport_scheme_method_of_lines, &
                                         transport_scheme_yz_bip_cosmic,   &
                                         transport_scheme_cosmic_3D,      &
                                         transport_scheme_horz_cosmic
  use xios
  use count_mod,                  only : count_type, halo_calls
  use mpi_mod,                    only : initialise_comm, store_comm, &
                                         finalise_comm, &
                                         get_comm_size, get_comm_rank

  implicit none

  private
  public initialise, run, finalise

  type(restart_type) :: restart

  ! Prognostic fields
  type( field_type ) :: u,         &
                        rho,       &
                        theta,     &
                        exner,     &
                        xi,        &
                        mr(nummr)

  ! Field collections
  type( field_collection_type ) :: derived_fields
  type( field_collection_type ) :: cloud_fields
  type( field_collection_type ) :: twod_fields

  ! Coordinate field
  type(field_type), target :: chi(3)

  integer(i_def) :: mesh_id      = imdi
  integer(i_def) :: twod_mesh_id = imdi

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration 
  subroutine initialise( filename )

    implicit none

    character(:), intent(in), allocatable :: filename

    character(len=*), parameter :: xios_id   = "lfric_client"
    character(len=*), parameter :: xios_ctx  = "gungho_atm"

    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: comm = -999
    integer(i_def) :: timestep, ts_init, dtime

    ! Initialse mpi and create the default communicator: mpi_comm_world
    call initialise_comm(comm)

    ! Initialise XIOS and get back the split mpi communicator
    call init_wait()
    call xios_initialize(xios_id, return_comm = comm)

    ! Save lfric's part of the split communicator for later use
    call store_comm(comm)

    ! Initialise YAXT
    call xt_initialize(comm)

    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging(local_rank, total_ranks, 'gungho')

    ! First call to log event must occur after the call to initialise_logging
    call log_event( 'gungho running...', LOG_LEVEL_INFO )

    call load_configuration( filename )
    call set_derived_config( .true. )

    restart = restart_type( restart_filename, local_rank, total_ranks )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) call timer('gungho')

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter('gungho')
    end if

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Create the mesh
    call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    write(log_scratch_space,'(A)') &
        "Purging global mesh collection."
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call global_mesh_collection%clear()
    deallocate(global_mesh_collection)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi)

    !-------------------------------------------------------------------------
    ! IO init
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up
    ! XIOS domain and context

    if ( (write_xios_output) .or. (restart%use_xios()) ) then

      dtime = int(dt)

      call xios_domain_init( xios_ctx,   &
                             comm,       &
                             dtime,      &
                             restart,    &
                             mesh_id,    &
                             chi)

    end if


    ! Create and initialise prognostic fields
    timestep = 0
    call init_gungho( mesh_id, chi, u, rho, theta, exner, mr, xi, restart )

    if (use_physics)then
      call init_physics(mesh_id, twod_mesh_id,                      &
                        u, exner, rho, theta,                       &
                        derived_fields, cloud_fields, twod_fields)

    end if

    ! Initial output
    ! We only want these once at the beginning of a run
    ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then

      ! Original nodal output
      if ( write_nodal_output)  then

        call output_nodal('theta', ts_init, theta, mesh_id)
        call output_nodal('xi',    ts_init, xi,    mesh_id)
        call output_nodal('u',     ts_init, u,     mesh_id)
        call output_nodal('rho',   ts_init, rho,   mesh_id)
        call output_nodal('exner', ts_init, exner, mesh_id)

        if (use_moisture)then
        call output_nodal('m_v',   ts_init, mr(imr_v),   mesh_id)
        call output_nodal('m_c',   ts_init, mr(imr_c),   mesh_id)
        call output_nodal('m_r',   ts_init, mr(imr_r),   mesh_id)
        call output_nodal('m_nc',  ts_init, mr(imr_nc),   mesh_id)
        call output_nodal('m_nr',  ts_init, mr(imr_nr),   mesh_id)
        end if

      end if

      ! XIOS output
      if (write_xios_output) then

        ! Make sure XIOS calendar is set to first timestep to be computed
        call xios_update_calendar(ts_init+1)

        ! Output scalar fields
        call output_xios_nodal("init_theta", theta, mesh_id)
        call output_xios_nodal("init_rho", rho, mesh_id)
        call output_xios_nodal("init_exner", exner, mesh_id)

        ! Output vector fields
        call output_xios_nodal("init_u", u, mesh_id)
        call output_xios_nodal("init_xi", xi, mesh_id)

        if (use_moisture)then
        call output_xios_nodal('init_m_v',  mr(imr_v),   mesh_id)
        call output_xios_nodal('init_m_c',  mr(imr_c),   mesh_id)
        call output_xios_nodal('init_m_r',  mr(imr_r),   mesh_id)
        call output_xios_nodal('init_m_nc', mr(imr_nc),   mesh_id)
        call output_xios_nodal('init_m_nr', mr(imr_nr),   mesh_id)
        end if

      end if

      if(write_minmax_tseries .and. ts_init == 0) then

         call minmax_tseries_init('u', mesh_id)
         call minmax_tseries(u, 'u', mesh_id)

      end if

      call divergence_diagnostic_alg(u, ts_init, mesh_id)
      call hydbal_diagnostic_alg(theta, exner, mesh_id)

    end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run()

    implicit none

    integer(i_def) :: timestep

    do timestep = restart%ts_start(),restart%ts_end()

      ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
      if ( (write_xios_output) .or. (restart%use_xios()) ) then
        call log_event( "Gungho: Updating XIOS timestep", LOG_LEVEL_INFO )
        call xios_update_calendar(timestep)
      end if

      write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if ( transport_only ) then

        select case( scheme )
          case ( transport_scheme_method_of_lines )
            if (timestep == restart%ts_start()) then
              ! Initialise and output initial conditions on first timestep
              call runge_kutta_init()
              call rk_transport_init( mesh_id, u, rho, theta)
            end if
            call rk_transport_step( u, rho, theta)
          case ( transport_scheme_yz_bip_cosmic )
            call cusph_cosmic_transport_init(mesh_id,u,timestep)
            call yz_bip_cosmic_step(rho,u,mesh_id,timestep)
          case ( transport_scheme_horz_cosmic )
            call cusph_cosmic_transport_init(mesh_id, u, timestep)
            call cusph_cosmic_transport_step( mesh_id, rho)
          case ( transport_scheme_cosmic_3D )
            call cusph_cosmic_transport_init(mesh_id, u, timestep)
            call cosmic_threed_transport_step( rho, u, mesh_id, timestep )
          case default
          call log_event("Dynamo: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
          stop
        end select
        call density_diagnostic_alg(rho, timestep)
        call conservation_algorithm(timestep, rho, u, theta, exner, xi)
      else
        select case( method )
          case( timestepping_method_semi_implicit )  ! Semi-Implicit
            ! Initialise and output initial conditions on first timestep
            if (timestep == restart%ts_start()) then
              call runge_kutta_init()
              call iter_alg_init(mesh_id, u, rho, theta, exner, mr, &
                                 twod_fields)
              call conservation_algorithm(timestep, rho, u, theta, exner, xi)
            end if
            call iter_alg_step(u, rho, theta, exner, mr, xi, &
                               derived_fields, cloud_fields, twod_fields,    &
                               timestep)

          case( timestepping_method_rk )             ! RK
            ! Initialise and output initial conditions on first timestep
            if (timestep == restart%ts_start()) then
              call runge_kutta_init()
              call rk_alg_init(mesh_id, u, rho, theta, exner)
              call conservation_algorithm(timestep, rho, u, theta, exner, xi)
            end if
            call rk_alg_step(u, rho, theta, exner, xi)
          case default
            call log_event("Dynamo: Incorrect time stepping option chosen, "// &
                            "stopping program! ",LOG_LEVEL_ERROR)
            stop
        end select

        call conservation_algorithm(timestep, rho, u, theta, exner, xi)

        if(write_minmax_tseries) call minmax_tseries(u, 'u', mesh_id)

       call u%log_minmax(LOG_LEVEL_INFO, ' u')

      end if

      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( mod(timestep, diagnostic_frequency) == 0 ) then

        ! Original nodal output
        if ( write_nodal_output)  then

          call log_event("Gungho: writing nodal output", LOG_LEVEL_INFO)
          call output_nodal("theta", timestep, theta, mesh_id)
          call output_nodal("rho",   timestep, rho,   mesh_id)
          call output_nodal("exner", timestep, exner, mesh_id)
          call output_nodal("u",     timestep, u,     mesh_id)
          call output_nodal("xi",    timestep, xi,    mesh_id)

          if (use_moisture)then
            call output_nodal('m_v',    timestep, mr(imr_v),   mesh_id)
            call output_nodal('m_c',    timestep, mr(imr_c),   mesh_id)
            call output_nodal('m_r',    timestep, mr(imr_r),   mesh_id)
            call output_nodal('m_nc',   timestep, mr(imr_nc),   mesh_id)
            call output_nodal('m_nr',   timestep, mr(imr_nr),   mesh_id)
          end if

        end if

        ! XIOS output
        if (write_xios_output) then

          ! Output scalar fields
          call log_event("Gungho: writing xios output", LOG_LEVEL_INFO)
          call output_xios_nodal("theta", theta, mesh_id)
          call output_xios_nodal("rho", rho, mesh_id)
          call output_xios_nodal("exner", exner, mesh_id)

          if (use_moisture)then
            call output_xios_nodal('m_v', mr(imr_v),   mesh_id)
            call output_xios_nodal('m_c', mr(imr_c),   mesh_id)
            call output_xios_nodal('m_r', mr(imr_r),   mesh_id)
            call output_xios_nodal('m_nc', mr(imr_nc),   mesh_id)
            call output_xios_nodal('m_nr', mr(imr_nr),   mesh_id)
          end if

          ! Output vector fields
          call output_xios_nodal("u", u, mesh_id)
          call output_xios_nodal("xi", xi, mesh_id)

        end if

        call divergence_diagnostic_alg(u, timestep, mesh_id)
        call hydbal_diagnostic_alg(theta, exner, mesh_id)
      end if

    end do ! end ts loop

    call runge_kutta_final()

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none

    integer(i_def) :: rc
    integer(i_def) :: ierr
    integer(i_def) :: i 

    character(5) :: name

    ! Log fields
    call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if (use_moisture)then
      call checksum_alg('gungho', rho, 'rho', theta, 'theta', u, 'u', &
         field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg('gungho', rho, 'rho', theta, 'theta', u, 'u')
    end if

    ! Write checkpoint files if required
    if( restart%checkpoint() ) then 

       write(log_scratch_space,'(A,I6)') "Checkpointing rho at ts ",restart%ts_end() 
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call rho%write_checkpoint("checkpoint_rho", trim(restart%endfname("rho")))

       write(log_scratch_space,'(A,I6)') "Checkpointing u at ts ",restart%ts_end() 
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call u%write_checkpoint("checkpoint_u", trim(restart%endfname("u")))

       write(log_scratch_space,'(A,I6)') "Checkpointing theta at ts ",restart%ts_end() 
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call theta%write_checkpoint("checkpoint_theta", trim(restart%endfname("theta")))

       write(log_scratch_space,'(A,I6)') "Checkpointing exner at ts ",restart%ts_end() 
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call exner%write_checkpoint("checkpoint_exner", trim(restart%endfname("exner")))
  
       write(log_scratch_space,'(A,I6)') "Checkpointing xi at ts ",restart%ts_end() 
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call xi%write_checkpoint("checkpoint_xi", trim(restart%endfname("xi")))

       if (use_moisture)then
         do i=1,nummr
           write(name, '(A,I2.2)') 'mr_', i
           write(log_scratch_space,'(A,A,A,I6)') "Checkpointing ",  trim(name), " at ts ",restart%ts_end() 
           call log_event(log_scratch_space,LOG_LEVEL_INFO)
           call mr(i)%write_checkpoint(name, trim(restart%endfname(name)))
         end do
       end if

    end if

    ! Call timestep finalizers
    if ( transport_only .and. scheme == transport_scheme_method_of_lines) then
      call rk_transport_final( rho, theta)
    end if

    call theta%field_final()
    call rho%field_final()
    call exner%field_final()
    call u%field_final()
    call xi%field_final()
    do i=1, nummr
      call mr(i)%field_final()
    end do

    if(write_minmax_tseries) call minmax_tseries_final(mesh_id)

    call iter_alg_final()
    call rk_alg_final()
    call final_runtime_constants()
    call final_configuration()

    if ( subroutine_timers ) then
      call timer('gungho')
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter('gungho')
      call halo_calls%output_counters()
    end if

    call log_event( 'gungho completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( (write_xios_output) .or. (restart%use_xios()) ) then
      call xios_context_finalize()
    end if

    if (allocated(mesh_collection)) then
      call mesh_collection%clear()
      deallocate(mesh_collection)
    end if

    if (allocated(function_space_collection)) then
      call function_space_collection%clear()
      deallocate(function_space_collection)
    end if

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise mpi and release the communicator
    call finalise_comm()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise

end module gungho_driver_mod
