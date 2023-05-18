!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Create fields for holding diagnostics
!> @details Initialise a diagnostic field

module initialise_diagnostics_mod

  use log_mod,                         only:                                  &
    log_event,                                                                &
    log_scratch_space,                                                        &
    log_level_info,                                                           &
    log_level_debug,                                                          &
    log_level_error
  use constants_mod,                   only: r_def, i_def, l_def, str_def
  use field_mod,                       only: field_type
  use mesh_mod,                        only: mesh_type
  use fs_continuity_mod,               only:                                  &
    W3,                                                                       &
    Wtheta,                                                                   &
    W2H,                                                                      &
    W0,                                                                       &
    name_from_functionspace
  use function_space_collection_mod,   only: function_space_collection
  use function_space_mod,              only: function_space_type
  use field_parent_mod,                only: write_interface
  use lfric_xios_write_mod,            only: write_field_generic
  use lfric_xios_mock_mod,             only: lfric_xios_mock_pull_in
  use lfric_xios_diag_mod,             only:                                  &
    field_is_active,                                                          &
    get_field_order,                                                          &
    get_field_grid_ref,                                                       &
    get_field_domain_ref,                                                     &
    get_field_axis_ref
  use multidata_field_dimensions_mod,  only:                                  &
    get_ndata                                                                 &
    => get_multidata_field_dimension
  use empty_data_mod,                  only: empty_real_data
  use extrusion_mod,                   only: TWOD
  use mesh_collection_mod,             only:                                  &
       mesh_collection,                                                       &
       mesh_collection_type
  use base_mesh_config_mod,            only: prime_mesh_name
  use io_config_mod,                   only: diag_always_on_sampling, use_xios_io

  implicit none

  private

  type(mesh_type), pointer :: diag_mesh_3d => null()
  type(mesh_type), pointer :: diag_mesh_2d => null()
  character(str_def), parameter :: diag_file_name = 'lfric_diag'

  ! field flavours
  character(str_def), parameter :: vanilla &
    = 'VanillaField'                                 ! scalar field on 3D mesh
  character(str_def), parameter :: vanilla_multi &
    = 'VanillaMulti'                              ! multidata field on 3D mesh
  character(str_def), parameter :: planar = 'PlanarField'
  character(str_def), parameter :: tile = 'TileField'
  character(str_def), parameter :: radiation = 'RadiationField'

  ! field status indicators
  character(str_def), parameter :: activated                                  &
    = 'Activated'     ! needed as a dependency
  character(str_def), parameter :: enabled                                    &
    = 'Enabled'       ! dynamically enabled, will be sampled
  character(str_def), parameter :: disabled                                   &
    = 'Disabled'      ! dynamically disabledm will not be sampled

  public :: init_diagnostic_field, diagnostic_to_be_sampled

contains

  ! derive function space enumerator from xios metadata
  function get_field_fsenum(field_name, grid_ref, domain_ref)                 &
    result(fsenum)
    implicit none
    character(*), intent(in) :: field_name
    character(*), intent(in) :: grid_ref
    character(*), intent(in) :: domain_ref
    integer(i_def) :: fsenum

    ! from RB's python metadata generator
    if (grid_ref == 'full_level_face_grid') then
      fsenum = Wtheta
    else if (grid_ref == 'half_level_face_grid'                               &
      .or. domain_ref == 'face') then
      fsenum = W3
    else if (grid_ref == 'half_level_edge_grid') then
      fsenum = W2H
    else if (grid_ref == 'node_grid') then
      fsenum = W0
    else
      fsenum = 0 ! silence compiler warning
      write(log_scratch_space, *)                                             &
        'cannot derive function space enumerator for field: ' //              &
        trim(field_name) //                                                   &
        ', grid_ref: ' // trim(grid_ref) //                                   &
        ', domain_ref: ' // trim(domain_ref)
      call log_event(log_scratch_space, log_level_error)
    end if
  end function get_field_fsenum

  ! derive field flavour (vanilla/planar/tile/radiation) from xios metadata
  function get_field_flavour(field_name, grid_ref, domain_ref, axis_ref)      &
    result(flavour)
    implicit none
    character(*), intent(in) :: field_name
    character(*), intent(in) :: grid_ref
    character(*), intent(in) :: domain_ref
    character(*), intent(in) :: axis_ref
    character(str_def) :: flavour

    if (grid_ref /= "") then
      if (domain_ref /= "") then
        write(log_scratch_space, *)                                           &
        'field ' // trim(field_name) //                                       &
        'with grid_ref and domain_ref : ' //                                  &
        grid_ref // ' ' // domain_ref
        call log_event(log_scratch_space, log_level_error)
      end if
      if (axis_ref /= "") then
        flavour = vanilla_multi
      else
        flavour = vanilla
      end if
    else
      if (domain_ref /= "") then
        if (axis_ref /= "") then
          if (axis_ref == 'radiation_levels') then
            flavour = radiation
          else
            flavour = tile
          end if
        else
          ! only domain - must not happen except for face domain
          if (domain_ref == 'face') then
            flavour = planar
          else
            write(log_scratch_space, *)                                       &
            'field ' // trim(field_name) //                                   &
            ' with only domain_ref: ' // domain_ref
            call log_event(log_scratch_space, log_level_error)
          end if
        end if
      else
        ! only axis - must not happen
        write(log_scratch_space, *)                                           &
        'field ' // trim(field_name) //                                       &
        ' with only an axis_ref: ' // axis_ref
        call log_event(log_scratch_space, log_level_error)
      end if
    end if
  end function get_field_flavour

  ! derive multidata item name from xios metadata
  function get_field_tile_id(field_name, axis_ref) result(tile_id)
    implicit none
    character(*), intent(in) :: field_name
    character(*), intent(in) :: axis_ref
    character(str_def) :: tile_id
    tile_id = axis_ref
  end function get_field_tile_id

  !> @brief Return true if and only if XIOS will sample the field.
  !> @param[in]   unique_id   XIOS id of field
  !> @return                  Sampling on/off status of the field
  function diagnostic_to_be_sampled(unique_id) result(sampling_on)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: sampling_on
    if (diag_always_on_sampling .or. .not. use_xios_io) then
      sampling_on = .true. ! for testing
      ! This is used when xios is off to ensure existing behvaiour of
      ! old nodal output is preserved
    else
      sampling_on = field_is_active(unique_id, .true.) ! derived from metadata
    end if
  end function diagnostic_to_be_sampled

  !> @brief Initialise a diagnostic field.
  !> @details If the field was requested, or if it is needed as a dependency,
  !> it will be created as an active field, otherwise with an empty data
  !> array (to save memory).
  !> Pass activate=.true. for fields needed as dependencies. If activate
  !> is not passed, the field's status will be derived from the XIOS metadata.
  !> Passing activate=.false. is equivalent to not passing activate.
  !> @post  The field name will be set equal to the XIOS id passed in.
  !> @param[out]          field            Field to initialise
  !> @param[in]           unique_id        XIOS id of field
  !> @param[in, optional] activate         Force-activate field
  !> @param[in, optional] force_mesh       Override derived mesh
  !> @param[in, optional] force_rad_levels Override derived radiation levels
  !> @return                               Sampling status of the field
  function init_diagnostic_field(field, unique_id,                            &
    activate, force_mesh, force_rad_levels) result (sampling_on)

    implicit none

    type(field_type),                   intent(out) :: field
    character(*),                       intent(in)  :: unique_id
    logical(l_def), optional,           intent(in)  :: activate
    type(mesh_type), pointer, optional, intent(in)  :: force_mesh
    integer(kind=i_def), optional,      intent(in)  :: force_rad_levels

    character(str_def), parameter :: routine_name = 'init_diagnostic_field'

    logical(kind=l_def) :: is_activated
    logical(kind=l_def) :: sampling_on
    logical(kind=l_def) :: active
    character(str_def)  :: field_name
    character(str_def)  :: grid_ref
    character(str_def)  :: domain_ref
    character(str_def)  :: axis_ref
    integer(kind=i_def) :: order
    integer(kind=i_def) :: fsenum
    character(str_def)  :: flavour
    integer(kind=i_def) :: ndata
    character(str_def)  :: status

    type(mesh_type), pointer            :: this_mesh    => null()
    type(function_space_type), pointer  :: vector_space => null()
    procedure(write_interface), pointer :: write_behaviour => null()

    ! dynamic initialisation
    if (.not. associated(diag_mesh_3d)) then
#ifdef UNIT_TEST
      diag_mesh_3d => mesh_collection%get_mesh('test mesh: planar bi-periodic')
      diag_mesh_2d => mesh_collection%get_mesh('test mesh: planar bi-periodic')
#else
      diag_mesh_3d => mesh_collection%get_mesh(prime_mesh_name)
      diag_mesh_2d &
        => mesh_collection%get_mesh_variant(diag_mesh_3d, extrusion_id=TWOD)
#endif
  end if

    ! the field_name is used only for logging and error reporting
    field_name = unique_id

    ! field sampling status
    sampling_on = diagnostic_to_be_sampled(unique_id)

    ! field activation status
    is_activated = .false.
    if (present(activate)) then
      if (activate) then
        is_activated = .true.
      end if
    end if
    if (is_activated) then
      active = .true.
      status = activated
    else
      active = sampling_on
      if (active) then
        status = enabled
      else
        status = disabled
      end if
    end if

    ! metadata lookup
    grid_ref = get_field_grid_ref(unique_id)
    domain_ref = get_field_domain_ref(unique_id)
    axis_ref = get_field_axis_ref(unique_id)
    order = get_field_order(unique_id )

    ! derive function space and flavour from metadata
    fsenum = get_field_fsenum(field_name, grid_ref, domain_ref)
    flavour = get_field_flavour(field_name, grid_ref, domain_ref, axis_ref)
    write_behaviour => write_field_generic

    ! select mesh
    if (present(force_mesh)) then
      this_mesh => force_mesh
    else if (                                                                 &
      flavour == radiation .or.                                               &
      flavour == tile .or.                                                    &
      flavour == planar) then
      this_mesh => diag_mesh_2d
    else if (flavour == vanilla .or. flavour == vanilla_multi) then
      this_mesh => diag_mesh_3d
    else
      call log_event( "unexpected flavour: " // flavour, log_level_error )
    end if

    ! derive ndata (multidata dimension)
    select case (flavour)
    case (radiation)
     ! radiation fields are treated as special multidata fields
      if (present(force_rad_levels)) then
        ndata = force_rad_levels
      else
        ndata = diag_mesh_3d%get_nlayers() + 1
      end if
    case (tile)
      ! genuine multidata field - derive dimension from axis metadata
      if (axis_ref == "") then
        ndata = 1
      else
        ndata = get_ndata(get_field_tile_id(field_name, axis_ref))
      end if
    case (planar)
      ! scalar field on 2d mesh
      ndata = 1
    case (vanilla )
      ! scalar field on 3d mesh
      ndata = 1
    case (vanilla_multi)
      ! multidata field on 3d mesh
      ndata = get_ndata(get_field_tile_id(field_name, axis_ref))
    case default
      call log_event("unexpected flavour: " // flavour, log_level_error)
    end select

    ! set up functions space - needed by psyclone even for inactive fields
    vector_space => function_space_collection%get_fs(                         &
      this_mesh,                                                              &
      order,                                                                  &
      fsenum,                                                                 &
      ndata)
#ifndef UNIT_TEST
    write(log_scratch_space,'(A, A, A, A, A, A, I2, A, A, A, A, A, I5)')      &
      'field: ', trim(field_name),                                            &
      ', ' // trim(status),                                                   &
      ', mesh: ', trim(this_mesh%get_mesh_name()),                            &
      ', order: ', order,                                                     &
      ', fs: ', trim(name_from_functionspace(fsenum)),                        &
      ', flavour: ', trim(flavour),                                           &
      ', ndata: ', ndata
    call log_event(log_scratch_space, log_level_debug)
#endif
    if (active) then ! field was requested or is needed as a dependency
      call field%initialise(                                                  &
        vector_space = vector_space,                                          &
        name = field_name)
      call field%set_write_behaviour(write_behaviour)
    else ! field is not needed
      ! ceate a field that points to the empty data array
      call field%initialise(vector_space = vector_space,                      &
       override_data = empty_real_data,                                       &
       name = field_name)
    end if

    ! check post-condition
    if (field%get_name() /= unique_id) then
      write(log_scratch_space,*) trim(routine_name) // ': ' //                &
       'name mismatch: ' // trim(field%get_name()) //' vs '// trim(unique_id)
       call log_event(log_scratch_space, log_level_error)
    end if

    ! paranoia
    nullify(this_mesh, vector_space, write_behaviour)
  end function init_diagnostic_field

end module initialise_diagnostics_mod
