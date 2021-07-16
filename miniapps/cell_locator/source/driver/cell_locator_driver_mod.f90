!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the cell_locator miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module cell_locator_driver_mod

  use constants_mod,                  only : i_def, i_native, r_def
  use cli_mod,                        only : get_initial_filename
  use cell_locator_mod,               only : load_configuration
  use create_mesh_mod,                only : init_mesh
  use create_fem_mod,                 only : init_fem
  use yaxt,                           only : xt_initialize, xt_finalize
  use field_mod,                      only : field_type
  use runtime_constants_mod,          only : create_runtime_constants
  use derived_config_mod,             only : set_derived_config
  use local_mesh_collection_mod,      only : local_mesh_collection, &
                                             local_mesh_collection_type
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             initialise_logging, &
                                             finalise_logging,  &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG,   &
                                             LOG_LEVEL_TRACE
  use mesh_collection_mod,            only : mesh_collection, &
                                             mesh_collection_type
  use mpi_mod,                        only : store_comm, &
                                             get_comm_size, get_comm_rank

  use cell_locator_api_mod,           only : cell_locator_api_type, &
                                             cell_locator_api_clear, &
                                             cell_locator_api_build, &
                                             cell_locator_api_find

  use cell_locator_config_mod,        only : verbose, output_filename

  use timer_mod,                      only : init_timer, timer, output_timer

  use, intrinsic :: iso_c_binding,    only : c_long_long

  implicit none

  private
  public initialise, run, finalise

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi_xyz
  type(field_type), target, dimension(3) :: chi_sph
  type(field_type), target               :: panel_id

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  type(cell_locator_api_type) :: cell_locator_obj

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

    implicit none

    character(:),      intent(in), allocatable :: filename
    integer(i_native), intent(in)              :: model_communicator

    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: ier

    ! Store the MPI communicator for later use.
    call store_comm( model_communicator )

    ! Initialise YAXT
    call xt_initialize( model_communicator )

    ! and get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging( local_rank, total_ranks, 'cell_locator' )

    call log_event( 'cell_locator: Running miniapp ...', LOG_LEVEL_INFO )

    call load_configuration( filename )
    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------

    call init_timer()
    call timer( 'cell_locator' )

    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )

    allocate( mesh_collection, &
              source=mesh_collection_type() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id, &
                    twod_mesh_id=twod_mesh_id )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_id, chi_xyz, chi_sph, panel_id )

    ! Sets up chi. Create runtime_constants object. This creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field.
    call create_runtime_constants( mesh_id, twod_mesh_id, chi_xyz, chi_sph, panel_id )

    ! Construct cell locator
    call timer( 'constructor' )
    cell_locator_obj = cell_locator_api_type( ier )
    if ( ier /= 0 ) then
      call log_event( 'cell_locator: Error after constructing cell locator', &
                      LOG_LEVEL_ERROR )
    endif
    call timer( 'constructor' )

    call timer( 'build' )
    call cell_locator_obj%build( mesh_id, ier )
    if ( ier /= 0 ) then
      call log_event( 'cell_locator: Error after building cell locator', &
                      LOG_LEVEL_ERROR )
    endif
    call timer( 'build' )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Find all the cells containing target points and write the results
  !>
  subroutine run()

    implicit none
    real(r_def)            :: target_point(3)
    integer(c_long_long)   :: cell_id_0 ! index of cell, zero based
    real(r_def)            :: pcoords(3)
    real(r_def)            :: dist_error_square
    character(len=1024)    :: msg
    integer(i_def)         :: ier, ipoint, npts, num_not_found, counter
    real(r_def)            :: total_error

    ! initially set to invalid
    cell_id_0 = -1

    npts = cell_locator_obj%get_number_of_target_points()

    counter = 1
    total_error = 0
    num_not_found = 0
    call timer( 'run' )
    do ipoint = 1, npts
      target_point = cell_locator_obj%get_target_point( ipoint )
      call cell_locator_obj%find( target_point, cell_id_0, pcoords, &
                                  dist_error_square, ier )
      if ( ier /= 0 ) then
        write ( msg, &
               '(A, I4, A, I8, A, 1PE11.4, A, 1PE11.4, A, 1PE11.4)' ) &
               'cell_locator: Error ier = ', ier, &
               ' after calling find in cell locator for ', &
               ipoint, '-th target_point = ', &
               target_point(1), ', ', target_point(2), ', ', target_point(3)
        call log_event( msg, LOG_LEVEL_ERROR )
      endif
      if (cell_id_0 >= 0) then
        if ( verbose > 1 ) then
          write( msg, &
            '(A, I8, A, F5.3, A, F5.3, A, F5.3, A, 1PE11.4, A,' &
            // ' 1PE11.4, A, 1PE11.4, A, 1PE9.2)' ) &
            'cell_locator: Found cell id ', cell_id_0, &
            ' pcoords = ', pcoords(1), ', ', pcoords(2), ', ', pcoords(3), &
            ' for target point = ', target_point(1), ', ', target_point(2), &
            ', ', target_point(3), &
            ' interp error sq = ', dist_error_square
          call log_event( msg, LOG_LEVEL_DEBUG )
        endif
        ! Store
        call cell_locator_obj%set_result(counter, ipoint, cell_id_0, &
                                         pcoords, dist_error_square)
        counter = counter + 1

        total_error = total_error + dist_error_square
      else
        if ( verbose > 1 ) then
          write( msg, '(A, 1PE11.4, A, 1PE11.4, A, 1PE11.4)' ) &
            'cell_locator: ~ Failed to find cell for target point = ', &
            target_point(1), ', ', target_point(2), ', ', target_point(3)
          call log_event( msg, LOG_LEVEL_TRACE )
        end if
        num_not_found = num_not_found + 1
      endif

    enddo
    call timer( 'run' )

    write( msg, '(A, I8, A, I8, A, F7.2, A, 1PE9.2)') &
           'cell_locator: Ratio of points not found = ', num_not_found, &
           '/', npts, ' = ', &
           100*real(num_not_found)/real(npts), &
           ' % total error = ', total_error
    call log_event( msg, LOG_LEVEL_INFO )

    call cell_locator_obj%write_results( output_filename, ier )
    if ( ier /= 0 ) then
      call log_event( &
        'cell_locator: Error occurred when writing results in file "' &
        //output_filename//'"', LOG_LEVEL_ERROR )
    endif

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    integer(i_def) :: ier

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    call log_event( 'cell_locator: Miniapp completed', LOG_LEVEL_INFO )


    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Write results and free cell locator object
    call timer( 'reclaiming memory' )
    call cell_locator_obj%clear( ier )
    if ( ier /= 0 ) then
      call log_event( &
        'cell_locator: Error occurred when calling clear', LOG_LEVEL_ERROR )
    endif
    call timer( 'reclaiming memory' )

    ! End of timing so write out timer output file
    call timer( 'cell_locator' )
    call output_timer()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise

end module cell_locator_driver_mod
