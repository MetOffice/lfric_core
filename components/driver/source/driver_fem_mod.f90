!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief    Initialisation and finalisation for FEM-specific choices for model.
!> @details  Contains routines related to FEM choices:
!>           * Populates the global function space collection with function spaces
!>             required by the model. Corresponding coordinate (chi) and panel_id
!>             inventories are also captured.
module driver_fem_mod

  use constants_mod,                 only: i_def, r_def, l_def, &
                                           str_def, imdi, cmdi
  use extrusion_mod,                 only: twod, prime_extrusion
  use fs_continuity_mod,             only: W0, W3, Wtheta, Wchi
  use function_space_mod,            only: function_space_type
  use function_space_collection_mod, only: function_space_collection
  use driver_coordinates_mod,        only: assign_coordinate_field
  use log_mod,                       only: log_event,       &
                                           log_level_info,  &
                                           log_level_error, &
                                           log_scratch_space
  use mesh_collection_mod,           only: mesh_collection
  use sci_chi_transform_mod,         only: init_chi_transforms, &
                                           final_chi_transforms

  ! Object types
  use config_mod, only: config_type
  use field_mod,  only: field_type
  use mesh_mod,   only: mesh_type
  use inventory_by_mesh_mod, only: inventory_by_mesh_type

  ! Configuration modules
  use base_mesh_config_mod,      only: geometry_spherical,    &
                                       geometry_planar,       &
                                       topology_non_periodic, &
                                       topology_fully_periodic
  use finite_element_config_mod, only: coord_system_xyz, &
                                       coord_space_W0,   &
                                       coord_space_Wchi, &
                                       coord_space_Wtheta

  implicit none

  private
  public :: init_fem, final_fem

contains

  !> @brief  Initialises the coordinate fields (chi) and FEM components.
  !>
  !> @param[in]      config               Application namelist configuration object
  !> @param[in,out]  chi_inventory        Inventory object, containing all of
  !!                                      the chi fields indexed by mesh
  !> @param[in,out]  panel_id_inventory   Inventory object, containing all of
  !!                                      the fields with the ID of mesh panels
  subroutine init_fem(config, chi_inventory, panel_id_inventory)

    implicit none

    ! Coordinate field
    type(config_type), intent(in) :: config

    type(inventory_by_mesh_type), intent(inout) :: chi_inventory
    type(inventory_by_mesh_type), intent(inout) :: panel_id_inventory

    character(str_def),    allocatable :: all_mesh_names(:)
    type(mesh_type),           pointer :: mesh
    type(mesh_type),           pointer :: twod_mesh
    type(field_type)                   :: chi(3)
    type(field_type)                   :: panel_id
    type(function_space_type), pointer :: fs

    integer(i_def) :: chi_space, coord, i
    integer(i_def) :: coord_order_h, coord_order_v
    integer(i_def) :: this_coord_order
    integer(i_def) :: halo_depth
    logical(l_def) :: is_valid

    character(str_def) :: mesh_name, prime_mesh_name
    integer(i_def)     :: geometry, topology, coord_system
    integer(i_def)     :: coord_space, coord_order, coord_order_nonprime
    real(r_def)        :: scaled_radius

    call log_event( 'FEM specifics: creating function spaces...', &
                    log_level_info )

    nullify(mesh, twod_mesh, fs)

    prime_mesh_name = cmdi
    if (config%namelist_exists('base_mesh')) then
      prime_mesh_name = config%base_mesh%prime_mesh_name()
    end if

    coord_system         = config%finite_element%coord_system()
    coord_order          = config%finite_element%coord_order()
    coord_space          = config%finite_element%coord_space()
    coord_order_nonprime = config%finite_element%coord_order_nonprime()
    scaled_radius        = config%planet%scaled_radius()

    ! ======================================================================== !
    ! Initialise coordinates
    ! ======================================================================== !

    ! To loop through mesh collection, get all mesh names
    ! Then get mesh from collection using these names
    all_mesh_names = mesh_collection%get_mesh_names()

    call chi_inventory%initialise(name="chi", table_len=size(all_mesh_names))
    call panel_id_inventory%initialise(name="panel_id", &
                                       table_len=size(all_mesh_names))

    ! ======================================================================== !
    ! Loop through all 3D meshes
    ! ======================================================================== !

    do i = 1, size(all_mesh_names)

      mesh => mesh_collection%get_mesh(all_mesh_names(i))
      mesh_name = mesh%get_mesh_name()

      if (mesh%is_geometry_spherical()) then
        geometry = geometry_spherical
      else
        geometry = geometry_planar
      end if

      if (mesh%is_topology_periodic()) then
        topology = topology_fully_periodic
      else
        topology = topology_non_periodic
      end if

      ! Initialise coordinate transformations
      call init_chi_transforms( geometry, topology, &
                                mesh_collection=mesh_collection )

      ! Only create coordinates for 3D meshes
      if (mesh%get_extrusion_id() /= twod) then

        ! Initialise panel ID field object -------------------------------------
        twod_mesh => mesh_collection%get_mesh(mesh, twod)
        fs => function_space_collection%get_fs(twod_mesh, 0, 0, W3)
        halo_depth = twod_mesh%get_halo_depth()
        call panel_id%initialise(fs, halo_depth=halo_depth)

        ! Initialise chi field object ------------------------------------------
        ! Set coordinate order for this mesh
        if (all_mesh_names(i) == prime_mesh_name) then
          this_coord_order = coord_order
        else
          this_coord_order = coord_order_nonprime
        end if

        ! Determine coordinate space
        select case (coord_space)
        case (coord_space_W0)
          ! Check domain/topology is valid
          is_valid = ( geometry     == geometry_spherical .and. &
                       coord_system == coord_system_xyz )       &
                .or. ( geometry     == geometry_planar .and.    &
                       topology     == topology_non_periodic )
          if (.not. is_valid) then
            write(log_scratch_space,'(A)')                               &
                'Coordinate space W0 is only valid for non-periodic ' // &
                'planar domains or when using the xyz coordinate system.'
            call log_event(log_scratch_space, log_level_error)
          end if

          ! Correct the coord_order for W0 polynomials being 1 order above W3
          this_coord_order = this_coord_order - 1
          chi_space = W0

        case (coord_space_Wchi)
          chi_space = Wchi

        case (coord_space_Wtheta)
          chi_space = Wtheta

        case default
          call log_event('Invalid value for coord_space', log_level_error)
        end select

        ! Set horizontal and vertical coordinate orders separately
        if (coord_system == coord_system_xyz) then
          ! Geocentric Cartesian coordinates - same order in all directions
          coord_order_h = this_coord_order
          coord_order_v = this_coord_order

        else
          ! For native coordinates, we separate horizontal and vertical coords
          ! and can still accurately represent space with linear vertical coords
          coord_order_h = this_coord_order

          if (coord_space == coord_space_Wchi) then
            coord_order_v = 1  ! Linear vertical coords for Wchi
          else
            coord_order_v = 0  ! Linear vertical coords for Wtheta
          end if
        end if

        ! Create coordinate space
        fs => function_space_collection%get_fs( mesh,          &
                                                coord_order_h, &
                                                coord_order_v, &
                                                chi_space )

        do coord = 1, size(chi)
          call chi(coord)%initialise(fs, halo_depth=halo_depth)
        end do

        ! Set coordinate fields --------------------------------------------------
        call assign_coordinate_field( config, mesh, chi, panel_id)

        ! Add fields to inventory
        call chi_inventory%copy_field_array(chi, mesh)
        call panel_id_inventory%copy_field(panel_id, mesh)

        nullify(mesh, fs)
      end if
    end do

    call log_event( 'FEM specifics created', log_level_info )

  end subroutine init_fem

  !> @brief  Finalises the function_space_collection.
  subroutine final_fem()

    implicit none

    call final_chi_transforms()

  end subroutine final_fem

end module driver_fem_mod
