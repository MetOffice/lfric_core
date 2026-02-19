!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module containing an object containing all data structures needed to test the
!> LFRic-XIOS component
module test_db_mod

  use calendar_mod,                   only: calendar_type
  use cli_mod,                        only: parse_command_line
  use config_loader_mod,              only: read_configuration
  use constants_mod,                  only: i_def, r_def, str_def, imdi, &
                                            r_second, i_timestep
  use extrusion_mod,                  only: TWOD
  use field_collection_mod,           only: field_collection_type
  use field_parent_mod,               only: read_interface, write_interface
  use field_mod,                      only: field_type, field_proxy_type
  use halo_comms_mod,                 only: initialise_halo_comms, &
                                            finalise_halo_comms
  use halo_routing_collection_mod,    only: halo_routing_collection_type, &
                                            halo_routing_collection
  use halo_comms_mod,                 only: halo_routing_type
  use function_space_collection_mod,  only: function_space_collection_type, &
                                            function_space_collection
  use lfric_mpi_mod,                  only: lfric_mpi_type, global_mpi, &
                                            create_comm, destroy_comm, &
                                            lfric_comm_type
  use local_mesh_mod,                 only: local_mesh_type
  use log_mod,                        only: initialise_logging, &
                                            finalise_logging,   &
                                            log_set_level, log_event, &
                                            LOG_LEVEL_TRACE, LOG_LEVEL_ERROR
  use config_mod,                     only: config_type
  use lfric_xios_read_mod,            only: read_field_generic
  use lfric_xios_write_mod,           only: write_field_generic
  use local_mesh_collection_mod,      only: local_mesh_collection_type, &
                                            local_mesh_collection
  use mesh_collection_mod,            only: mesh_collection_type, &
                                            mesh_collection
  use mesh_mod,                       only: mesh_type, PLANE, PLANE_TWOD
  use model_clock_mod,                only: model_clock_type
  use function_space_mod,             only: function_space_type
  use fs_continuity_mod,              only: Wchi, W0, W2H, W3
  use step_calendar_mod,              only: step_calendar_type

  implicit none

  !> Object containing infrastructure for testing LFRic-XIOS
  type, public :: test_db_type
    private
    type(lfric_comm_type), public :: comm
    type(config_type),     public :: config
    type(field_type),      public :: chi(3)
    type(field_type),      public :: panel_id
    type(model_clock_type), public, allocatable :: clock
    class(calendar_type),   public, allocatable :: calendar
    type(field_collection_type), public :: temporal_fields

  contains
    procedure initialise
    procedure finalise
  end type

contains

  !> Initialiser for test db object
  subroutine initialise(self)

    implicit none

    class(test_db_type), target, intent(inout) :: self

    type(local_mesh_type), target  :: local_mesh
    type(local_mesh_type), pointer :: local_mesh_ptr
    type(mesh_type), target :: mesh, twod_mesh
    type(mesh_type), pointer :: mesh_ptr
    type(mesh_type), pointer :: twod_mesh_ptr
    type(function_space_type), pointer :: wchi_fs
    type(function_space_type), pointer :: tmp_fs
    type(field_proxy_type) :: chi_p(3), pid_p, rproxy
    type(field_type) :: temporal_field

    procedure(read_interface),  pointer :: read_ptr
    procedure(write_interface), pointer :: write_ptr

    integer(i_def), parameter :: n_x = 3
    integer(i_def), parameter :: n_y = 3
    integer(i_def), parameter :: n_z = 5

    integer(i_def) :: local_mesh_id, mesh_id, twod_mesh_id, coord, i, x_start, y_start, rc, cells_per_layer

    character(:), allocatable :: filename

    integer(i_timestep) :: first
    integer(i_timestep) :: last
    character(str_def) :: start_date, calendar_start, calendar_origin
    character(str_def) :: timestep_start
    character(str_def) :: timestep_end

    real(r_second) :: timestep_length

    call parse_command_line( filename )

    ! Initialise MPI & logging
    call create_comm(self%comm)
    call global_mpi%initialise(self%comm)
    call initialise_halo_comms(self%comm)
    call initialise_logging(self%comm%get_comm_mpi_val(), 'lfric_xios_context_test')
    call log_set_level(LOG_LEVEL_TRACE)

    call self%config%initialise("lfric_xios_integration_tests")
    call read_configuration(trim(adjustl(filename)), config=self%config)

    deallocate(filename)

    start_date      = self%config%time%calendar_start()
    timestep_start  = self%config%time%timestep_start()
    timestep_end    = self%config%time%timestep_end()
    timestep_length = self%config%timestepping%dt()

    ! Create top level mesh collection, function spaces & routing tables
    local_mesh_collection = local_mesh_collection_type()
    mesh_collection = mesh_collection_type()
    function_space_collection = function_space_collection_type()
    halo_routing_collection = halo_routing_collection_type()

    ! Dummy mesh mod has 9 cells, 3 layers and is uniform in vertical
    call local_mesh%initialise()
    local_mesh_id = local_mesh_collection%add_new_local_mesh( local_mesh )
    local_mesh_ptr => local_mesh_collection%get_mesh_by_id(local_mesh_id)

    mesh = mesh_type(PLANE, local_mesh_ptr)
    mesh_id = mesh_collection%add_new_mesh(mesh)
    mesh_ptr => mesh_collection%get_mesh_by_id(mesh_id)

    twod_mesh = mesh_type(PLANE_TWOD, local_mesh_ptr)
    twod_mesh_id = mesh_collection%add_new_mesh(twod_mesh)
    twod_mesh_ptr => mesh_collection%get_mesh_by_id(twod_mesh_id)

    ! Create function spaces
    wchi_fs => function_space_collection%get_fs(mesh_ptr, 0, 0, WChi)

    tmp_fs => function_space_collection%get_fs(mesh_ptr, 0, 0, W0)
    tmp_fs => function_space_collection%get_fs(mesh_ptr, 0, 0, W2H)
    tmp_fs => function_space_collection%get_fs(mesh_ptr, 0, 0, W3)
    tmp_fs => function_space_collection%get_fs(twod_mesh_ptr, 0, 0, W0)
    tmp_fs => function_space_collection%get_fs(twod_mesh_ptr, 0, 0, W2H)
    tmp_fs => function_space_collection%get_fs(twod_mesh_ptr, 0, 0, W3)

    ! Create coordinate fields
    do coord = 1, size(self%chi)
      call self%chi(coord)%initialise(vector_space = wchi_fs, halo_depth = 1 )
      chi_p(coord) = self%chi(coord)%get_proxy()
    end do

    call self%panel_id%initialise(vector_space = wchi_fs, halo_depth = 1 )
    pid_p = self%panel_id%get_proxy()
    pid_p%data(:) = 1.0_r_def

    ! Assign coordinate values
    cells_per_layer = size(chi_p(1)%data) / mesh%get_nlayers()

    do coord = 1, n_x
      do i = 1, mesh%get_nlayers()
        x_start =  1 + (coord - 1) * size(self%chi) + (i - 1) * cells_per_layer
        chi_p(1)%data(x_start:x_start + (size(self%chi) - 1)) = coord
      end do
    end do

    do coord = 1, n_y
      y_start = 1 + (coord - 1) * n_x * n_z
      chi_p(2)%data(y_start:y_start + (n_x * n_z) - 1) = coord
    end do

    do coord = 1, n_z
      chi_p(3)%data(coord::n_z) = real(coord, r_def)
    end do

    ! Set halos to clean
    do coord = 1, size(self%chi)
      call chi_p(coord)%set_clean(depth = 1)
    end do

    ! Init clock and calendar
    calendar_origin = start_date
    calendar_start = start_date
    allocate( self%calendar,                                                &
              source = step_calendar_type( trim(adjustl(calendar_origin)),  &
                                           trim(adjustl(calendar_start)) ), &
                                           stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate calendar", LOG_LEVEL_ERROR )
    end if


    first = self%calendar%parse_instance(timestep_start)
    last  = self%calendar%parse_instance(timestep_end)

    if ( allocated(self%clock) ) deallocate (self%clock)
    allocate( self%clock,                                 &
              source = model_clock_type( first, last,     &
                                         timestep_length, &
                                         0.0_r_second ),  &
                                         stat=rc )

    if (rc /= 0) then
      call log_event( "Unable to allocate model clock", LOG_LEVEL_ERROR )
    end if

    ! Create field for reading
    call self%temporal_fields%initialise(name="temporal_fields", table_len=1)
    tmp_fs => function_space_collection%get_fs(mesh_ptr, 0, 0, W3)
    call temporal_field%initialise(vector_space = tmp_fs, name="temporal_field" )
    rproxy = temporal_field%get_proxy()
    rproxy%data(:) = 0.0_r_def
    read_ptr => read_field_generic
    call temporal_field%set_read_behaviour(read_ptr)
    write_ptr => write_field_generic
    call temporal_field%set_write_behaviour(write_ptr)
    call self%temporal_fields%add_field(temporal_field)


    nullify(local_mesh_ptr)
    nullify(mesh_ptr)
    nullify(twod_mesh_ptr)
    nullify(wchi_fs)
    nullify(tmp_fs)

  end subroutine initialise

  !> Destroy the test infrastructure
  subroutine finalise( self )

    implicit none

    class(test_db_type), intent(inout) :: self

    call finalise_logging()
    call global_mpi%finalise()
    call destroy_comm()

  end subroutine finalise

end module test_db_mod
