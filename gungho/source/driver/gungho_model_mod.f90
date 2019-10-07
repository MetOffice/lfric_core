!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief initialise and finalise functionality for gungho simulations

!> @details Handles initialise and finalise of prognostic and coordinate fields

module gungho_model_mod

  use checksum_alg_mod,               only : checksum_alg
  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def
  use field_collection_mod,           only : field_collection_type, &
                                             field_collection_iterator_type
  use field_mod,                      only : field_type, &
                                             write_interface
  use formulation_config_mod,         only : transport_only, &
                                             use_moisture,   &
                                             use_physics
  use io_config_mod,                  only : write_diag, &
                                             write_dump, &
                                             write_minmax_tseries
  use io_mod,                         only : write_state, &
                                             xios_write_field_single_face
  use iter_timestep_alg_mod,          only : iter_alg_init, &
                                             iter_alg_final
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_DEBUG, &
                                             LOG_LEVEL_INFO,  &
                                             LOG_LEVEL_ERROR
  use minmax_tseries_mod,             only : minmax_tseries,      &
                                             minmax_tseries_init, &
                                             minmax_tseries_final
  use mr_indices_mod,                 only : nummr
  use rk_alg_timestep_mod,            only : rk_alg_init, &
                                             rk_alg_final
  use rk_transport_mod,               only : rk_transport_init, &
                                             rk_transport_final
  use runge_kutta_init_mod,           only : runge_kutta_init, &
                                             runge_kutta_final
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk
  use transport_config_mod,           only : scheme, &
                                             scheme_method_of_lines
  implicit none

  private 
  public initialise_model, finalise_model

  contains

  !> @brief Initialises the gungho app
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] prognostic_fields A collection of all the prognostic fields
  !> @param[inout] diagnostic_fields A collection of all the diagnostic fields
  !> @param[inout] mr Array of fields containing the mixing ratios
  !> @param[inout] twod_fields 2D field collection for physics
  !> @param[in] timestep_start number of timestep at which this run started
  subroutine initialise_model( mesh_id,           &
                               prognostic_fields, &
                               diagnostic_fields, &
                               mr,                &
                               twod_fields,       &
                               timestep_start )
    implicit none

    integer(i_def),                intent(in)    :: mesh_id
    type( field_collection_type ), intent(inout) :: prognostic_fields
    type( field_collection_type ), intent(inout) :: diagnostic_fields
    type( field_type ),            intent(inout) :: mr(nummr)
    type( field_collection_type ), intent(inout) :: twod_fields
    integer(i_def),                intent(in)    :: timestep_start

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    if (write_minmax_tseries) then
      call minmax_tseries_init('u', mesh_id)
      call minmax_tseries(u, 'u', mesh_id)
    end if

    if ( transport_only ) then

      select case( scheme )
        case ( scheme_method_of_lines )
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_transport_init( mesh_id, u, rho, theta)
        case default
          call log_event("Gungho: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select
    else
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call iter_alg_init(mesh_id, u, rho, theta, exner, mr, &
                             twod_fields)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner)
        case( method_rk )             ! RK
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_alg_init(mesh_id, u, rho, theta, exner)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner)
        case default
          call log_event("Gungho: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

  end subroutine initialise_model


  !> @brief Finalise the gungho application
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] prognostic_fields A collection containing
  !>                  the prognostic fields
  !> @param[inout] diagnostic_fields A collection containing
  !>                  the diagnostic fields
  !> @param[inout] mr Array of fields containing the mixing ratios
  !> @param[inout] twod_fields Collection of two dimensional physics fields
  !> @param[inout] fd_fields Collection of a copy of the prognostic field
  !>                       in finite-difference form
  !> @param [in] program_name An identifier given to the model begin run
  subroutine finalise_model( mesh_id,           &
                             prognostic_fields, &
                             diagnostic_fields, &
                             mr,                &
                             twod_fields,       &
                             fd_fields,         &
                             program_name )

    implicit none

    integer(i_def),                intent(in)    :: mesh_id
    type( field_collection_type ), intent(inout) :: prognostic_fields
    type( field_collection_type ), intent(inout) :: diagnostic_fields
    type( field_type ),            intent(inout) :: mr(nummr)
    type( field_collection_type ), intent(inout) :: twod_fields
    type( field_collection_type ), intent(inout) :: fd_fields
    character(*),                  intent(in)    :: program_name

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Pointer for tstar_2d to allow write to dump
    type( field_type ), pointer   :: tstar_2d => null()
    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

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

    !===================== Write fields to dump ======================!

    if( write_dump ) then

      ! Current dump writing is only relevant for physics runs at the moment
      if (use_physics) then

        ! For the purposes of dumping from one collection, we add a pointer
        ! to tstar to the fd_prognostics collection

        tmp_write_ptr => xios_write_field_single_face
        tstar_2d => twod_fields%get_field('tstar')
        call tstar_2d%set_write_behaviour(tmp_write_ptr)

        call fd_fields%add_field(tstar_2d)

        call log_event("Gungho: writing FD fields to dump", LOG_LEVEL_INFO)

        ! Write prognostic fields to dump
        call write_state(fd_fields)
      
        nullify(tmp_write_ptr, tstar_2d)

      end if

    end if

    ! Call timestep finalizers
    if ( transport_only .and. scheme == scheme_method_of_lines) then
      call rk_transport_final( rho, theta)
    end if

    if(write_minmax_tseries) call minmax_tseries_final(mesh_id)

    call runge_kutta_final()

    if ( .not. transport_only ) then
      if ( method == method_semi_implicit ) call iter_alg_final()
      if ( method == method_rk )            call rk_alg_final()
    end if

  end subroutine finalise_model

end module gungho_model_mod
