!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls various init subroutines to create a mesh, function spaces
!> and then prognostic fields on those function spaces.
!> Controls the main timestepping loop, initialises and calls
!> timestepping algorithms and also handles output of diagnostics at
!> specified intervals as well as checkpoint restart files.

program dynamo

  use cli_mod,                        only : get_initial_filename
  use constants_mod,                  only : i_def
  use dynamo_mod,                     only : load_configuration
  use ESMF
  use field_io_mod,                   only : write_state_netcdf
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use formulation_config_mod,         only : transport_only, use_moisture
  use init_gungho_mod,                only : init_gungho
  use init_dynamo_mod,                only : init_dynamo
  use io_utility_mod,                 only : open_file, close_file
  use iter_timestep_alg_mod,          only : iter_alg_init, &
                                             iter_alg_step
  use runge_kutta_init_mod,           only : runge_kutta_init
  use rk_alg_timestep_mod,            only : rk_alg_init, &
                                             rk_alg_step
  use transport_config_mod,           only : scheme, &
                                             transport_scheme_method_of_lines, &
                                             transport_scheme_bip_cosmic, &
                                             transport_scheme_cusph_cosmic
  use rk_transport_mod,               only : rk_transport_init, &
                                             rk_transport_step, &
                                             rk_transport_final
  use cosmic_transport_alg_mod,       only : cosmic_transport_init, &
                                             cosmic_transport_step
  use conservation_algorithm_mod,     only : conservation_algorithm
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG,   &
                                             LOG_LEVEL_TRACE
  use restart_config_mod,             only : restart_filename => filename
  use restart_control_mod,            only : restart_type
  use output_config_mod,              only : diagnostic_frequency
  use output_alg_mod,                 only : output_alg
  use timestepping_config_mod,        only : method, &
                                           timestepping_method_semi_implicit, &
                                           timestepping_method_rk
  use derived_config_mod,             only : set_derived_config
  use runtime_constants_mod,          only : create_runtime_constants, &
                                             get_geopotential
  use checksum_alg_mod,               only : checksum_alg
  use diagnostic_alg_mod,             only : divergence_diagnostic_alg
  use mr_indices_mod,                 only : imr_v, imr_c, imr_r, imr_nc, &
                                             imr_nr, nummr
  use set_rho_alg_mod,                only : set_rho_alg

  implicit none

  character(:), allocatable :: filename
  integer                   :: namelist_unit

  type(ESMF_VM)      :: vm
  integer            :: rc
  integer            :: total_ranks, local_rank
  integer            :: petCount, localPET

  type(restart_type) :: restart

  integer(i_def)     :: mesh_id

  ! coordinate fields
  type( field_type ) :: chi(3)

  ! prognostic fields
  type( field_type ) :: u, rho, theta, xi, mr(nummr), rho_in_wth

  ! Array to hold fields for checkpoint output
  type( field_type ), allocatable   :: checkpoint_output(:)

  ! temps to hold things retrieved from runtime_constants
  ! that are needed for output
  type( field_type )               :: geopotential


  integer                          :: timestep, ts_init

  integer :: i
  character(5) :: name

  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise ESMF and get the rank information from the virtual machine
  CALL ESMF_Initialize(vm=vm, defaultlogfilename="dynamo.Log", &
                  logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) &
    call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  total_ranks = petCount
  local_rank  = localPET

  ! Currently log_event can only use ESMF so it cannot be used before ESMF is
  ! initialised.
  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  allocate( filename, source='dynamo_configuration.nml' )
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call load_configuration( namelist_unit )
  call set_derived_config()
  call close_file( namelist_unit )
  deallocate( filename )

  restart = restart_type( restart_filename, local_rank, total_ranks )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------

  ! Create the mesh and function space collection
  call init_gungho(mesh_id, local_rank, total_ranks)

  ! Create and initialise prognostic fields
  timestep = 0
  call init_dynamo(mesh_id, chi, u, rho, theta, rho_in_wth, mr, xi, restart)

  ! Create runtime_constants object. This in turn creates various things
  ! needed by the timestepping algorithms such as mass matrix operators, mass
  ! matrix diagonal fields and the geopotential field

  call create_runtime_constants(mesh_id, chi)

  geopotential = get_geopotential()

  ! Initial output
  ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.
  call output_alg('theta', ts_init, theta, mesh_id)
  call output_alg('xi',    ts_init, xi,    mesh_id)
  call output_alg('u',     ts_init, u,     mesh_id)
  call output_alg('rho',   ts_init, rho,   mesh_id)
  call divergence_diagnostic_alg(u, ts_init, mesh_id)

  if (use_moisture)then
    call output_alg('m_v',   ts_init, mr(imr_v),   mesh_id)
    call output_alg('m_c',   ts_init, mr(imr_c),   mesh_id)
    call output_alg('m_r',   ts_init, mr(imr_r),   mesh_id)
    call output_alg('m_nc',   ts_init, mr(imr_nc),   mesh_id)
    call output_alg('m_nr',   ts_init, mr(imr_nr),   mesh_id)
  end if
  
  !-----------------------------------------------------------------------------
  ! model step 
  !-----------------------------------------------------------------------------
  do timestep = restart%ts_start(),restart%ts_end()
    call log_event( &
    "/****************************************************************************\ ", &
     LOG_LEVEL_TRACE )
     write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
     call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if ( transport_only ) then

      select case( scheme )
        case ( transport_scheme_method_of_lines)
          if (timestep == restart%ts_start()) then 
            ! Initialise and output initial conditions on first timestep
            call runge_kutta_init()
            call rk_transport_init( mesh_id, u, rho, theta)
            call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
          end if
          call rk_transport_step( mesh_id, chi, u, rho, theta)
        case ( transport_scheme_bip_cosmic)
          if (timestep == restart%ts_start()) then 
            ! Initialise and output initial conditions on first timestep
            call cosmic_transport_init(mesh_id, u)
            call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
          end if
          call cosmic_transport_step( mesh_id, chi, rho)
        case ( transport_scheme_cusph_cosmic)
          call set_rho_alg(mesh_id, chi, rho, timestep)
        case default
         call log_event("Dynamo: Incorrect transport option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
         stop
      end select
    else
       select case( method )
         case( timestepping_method_semi_implicit )  ! Semi-Implicit 
           ! Initialise and output initial conditions on first timestep
           if (timestep == restart%ts_start()) then 
             call runge_kutta_init()
             call iter_alg_init(mesh_id, u, rho, theta)
             call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
             call conservation_algorithm(timestep, mesh_id, rho, u, theta, xi, geopotential, chi)
           end if
           call iter_alg_step(chi, u, rho, theta, mr, xi, timestep)

         case( timestepping_method_rk )             ! RK
           ! Initialise and output initial conditions on first timestep
           if (timestep == restart%ts_start()) then
             call runge_kutta_init()
             call rk_alg_init( mesh_id, u, rho, theta)
             call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
             call conservation_algorithm(timestep, mesh_id, rho, u, theta, xi, geopotential, chi)
           end if
           call rk_alg_step( mesh_id, chi, u, rho, theta, xi)
         case default
           call log_event("Dynamo: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
           stop
       end select

       call conservation_algorithm(timestep, mesh_id, rho, u, theta, xi, geopotential, chi)

    end if

    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
     LOG_LEVEL_INFO)


    ! Use diagnostic output frequency to determine whether to write diagnostics
    ! on this timestep

    if ( mod(timestep, diagnostic_frequency) == 0 ) then
      call log_event("Dynamo: writing diagnostic output", LOG_LEVEL_INFO)
      call output_alg('theta', timestep, theta, mesh_id)
      call output_alg('xi',    timestep, xi,    mesh_id)
      call output_alg('u',     timestep, u,     mesh_id)
      call output_alg('rho',   timestep, rho,   mesh_id)
      call divergence_diagnostic_alg(u, timestep, mesh_id)
      if (use_moisture)then
        call output_alg('m_v',    timestep, mr(imr_v),   mesh_id)
        call output_alg('m_c',    timestep, mr(imr_c),   mesh_id)
        call output_alg('m_r',    timestep, mr(imr_r),   mesh_id)
        call output_alg('m_nc',    timestep, mr(imr_nc),   mesh_id)
        call output_alg('m_nr',    timestep, mr(imr_nr),   mesh_id)
      end if
    end if


  end do ! end ts loop

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Log fields
  call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
  call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
  call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

  ! Write checksums to file
  call checksum_alg('dynamo', rho, 'rho', theta, 'theta', u, 'u')

  ! Write checkpoint/restart files if required
  if( restart%write_file() ) then 

    allocate(checkpoint_output(4))
    write(log_scratch_space,'(A,A)') "writing file:",  &
                            trim(restart%endfname("rho"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     checkpoint_output(1) = rho
     call write_state_netcdf( 1, checkpoint_output(1), trim(restart%endfname("rho")) )

     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("u"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     checkpoint_output(2) = u
     call write_state_netcdf( 1, checkpoint_output(2), trim(restart%endfname("u")) )

     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("theta"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     checkpoint_output(3) = theta
     call write_state_netcdf( 1, checkpoint_output(3), trim(restart%endfname("theta")) )
  
     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("xi"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     checkpoint_output(4) = xi
     call write_state_netcdf( 1, checkpoint_output(4), trim(restart%endfname("xi")) )

     if (use_moisture)then
       do i=1,nummr
         write(name, '(A,I2.2)') 'mr_', i
         write(log_scratch_space,'(A,A)') "writing file:",  &
            trim(restart%endfname(name))
         call log_event(log_scratch_space,LOG_LEVEL_INFO)
         checkpoint_output(1) = mr(i)
         call write_state_netcdf( 1, checkpoint_output(1), trim(restart%endfname(name)) )
       end do
     end if
  end if

  ! Call timestep finalizers
  if ( transport_only .and. scheme == transport_scheme_method_of_lines) then
    call rk_transport_final( mesh_id, rho, theta)
  end if

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Close down ESMF
  call ESMF_Finalize(rc=rc)

end program dynamo
