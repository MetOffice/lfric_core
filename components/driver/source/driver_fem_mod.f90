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
!>           * Initialises function space chains for use by the model.
module driver_fem_mod

  use config_mod,                     only: config_type
  use constants_mod,                  only: i_def, r_def, l_def, str_def, imdi
  use extrusion_mod,                  only: TWOD, PRIME_EXTRUSION
  use field_mod,                      only: field_type
  use fs_continuity_mod,              only: W0, W2, W3, Wtheta, Wchi, W2v, W2h
  use function_space_mod,             only: function_space_type
  use function_space_collection_mod,  only: function_space_collection
  use function_space_chain_mod,       only: function_space_chain_type,          &
                                            single_layer_function_space_chain,  &
                                            multigrid_function_space_chain,     &
                                            W2_multigrid_function_space_chain,  &
                                            W2v_multigrid_function_space_chain, &
                                            W2h_multigrid_function_space_chain, &
                                            wtheta_multigrid_function_space_chain
  use driver_coordinates_mod,         only: assign_coordinate_field
  use inventory_by_mesh_mod,          only: inventory_by_mesh_type
  use log_mod,                        only: log_event,         &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_ERROR,   &
                                            log_scratch_space
  use mesh_mod,                       only: mesh_type
  use mesh_collection_mod,            only: mesh_collection

  use sci_chi_transform_mod, only: init_chi_transforms, &
                                   final_chi_transforms

  implicit none

  private
  public :: init_fem, init_function_space_chains, final_fem

contains

  !> @brief  Initialises the coordinate fields (chi) and FEM components.
  !>
  !> @param[in]      config               Configuration object
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
    integer(i_def)                     :: chi_space, coord, i

    character(str_def) :: mesh_name
    integer(i_def)     :: coord_order, geometry, topology, coord_system
    real(r_def)        :: scaled_radius

    call log_event( 'FEM specifics: creating function spaces...', &
                    log_level_info )

    nullify(mesh, twod_mesh, fs)

    if (config%namelist_exists('base_mesh')) then
      geometry = config%base_mesh%geometry()
      topology = config%base_mesh%topology()
    else
      geometry = imdi
      topology = imdi
    end if

    coord_system  = config%finite_element%coord_system()
    coord_order   = config%finite_element%coord_order()
    scaled_radius = config%planet%scaled_radius()

    ! ======================================================================== !
    ! Initialise coordinates
    ! ======================================================================== !

    ! Initialise coordinate transformations
    call init_chi_transforms( geometry, topology, mesh_collection=mesh_collection )

    ! To loop through mesh collection, get all mesh names
    ! Then get mesh from collection using these names
    all_mesh_names = mesh_collection%get_mesh_names()

    call chi_inventory%initialise(name="chi", table_len=size(all_mesh_names))
    call panel_id_inventory%initialise(name="panel_id", &
                                       table_len=size(all_mesh_names))

    ! ======================================================================== !
    ! Loop through all 3D meshes
    ! ======================================================================== !

    do i = 1, SIZE(all_mesh_names)
      mesh => mesh_collection%get_mesh(all_mesh_names(i))
      mesh_name = mesh%get_mesh_name()

      ! Only create coordinates for 3D meshes
      if (mesh%get_extrusion_id() /= TWOD) then

        ! Initialise panel ID field object ---------------------------------------
        twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
        fs => function_space_collection%get_fs(twod_mesh, 0, 0, W3)
        call panel_id%initialise( vector_space=fs, &
                                  halo_depth=twod_mesh%get_halo_depth() )

        ! Initialise chi field object --------------------------------------------
        if ( coord_order == 0 ) then
          chi_space = W0
          write(log_scratch_space,'(A)') &
              'Computing W0 coordinate fields for ' // &
              trim(mesh_name) // 'mesh'
          call log_event( log_scratch_space, log_level_info )
        else
          chi_space = Wchi
          write(log_scratch_space,'(A)') &
              'Computing Wchi coordinate fields for ' // &
              trim(mesh_name) // 'mesh'
          call log_event( log_scratch_space, log_level_info )
        end if
        fs => function_space_collection%get_fs( mesh, coord_order, &
                                                coord_order, chi_space )

        do coord = 1, size(chi)
          call chi(coord)%initialise( vector_space=fs, &
                                      halo_depth=twod_mesh%get_halo_depth() )
        end do

        ! Set coordinate fields --------------------------------------------------
        call assign_coordinate_field(chi, panel_id, mesh, &
                                     geometry, topology, &
                                     coord_system, scaled_radius)

        ! Add fields to inventory
        call chi_inventory%copy_field_array(chi, mesh)
        call panel_id_inventory%copy_field(panel_id, mesh)

        nullify(mesh, fs)
      end if
    end do

    call log_event( 'FEM specifics created', log_level_info )

  end subroutine init_fem

  !> @brief  Initialises the function space chains used in multigrid.
  !> @param[in] multigrid_mesh_names  Names of the multigrid meshes
  subroutine init_function_space_chains(multigrid_mesh_names)

    implicit none

    character(str_def),  intent(in) :: multigrid_mesh_names(:)

    type(mesh_type), pointer :: mesh
    type(mesh_type), pointer :: twod_mesh

    type(function_space_type), pointer :: fs

    integer(i_def) :: i

    nullify(mesh, twod_mesh, fs)

    call log_event( 'FEM specifics: creating function space chains...', &
                    log_level_info )

    ! ======================================================================== !
    ! Create function space chains
    ! ======================================================================== !

    multigrid_function_space_chain        = function_space_chain_type()
    w2_multigrid_function_space_chain     = function_space_chain_type()
    w2v_multigrid_function_space_chain    = function_space_chain_type()
    w2h_multigrid_function_space_chain    = function_space_chain_type()
    wtheta_multigrid_function_space_chain = function_space_chain_type()

    write(log_scratch_space,'(A,I1,A)') &
        'Initialising MultiGrid ', size(multigrid_mesh_names), &
        '-level function space chain.'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    do i = 1, size(multigrid_mesh_names)

      mesh => mesh_collection%get_mesh( multigrid_mesh_names(i) )

      ! Make sure this function_space is in the collection
      fs => function_space_collection%get_fs( mesh, 0, 0, W3 )
      call multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, W2 )
      call w2_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, W2v )
      call w2v_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, W2h )
      call w2h_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, Wtheta )
      call wtheta_multigrid_function_space_chain%add( fs )
    end do

    single_layer_function_space_chain = function_space_chain_type()
    do i = 1, size(multigrid_mesh_names)
      mesh => mesh_collection%get_mesh( multigrid_mesh_names(i) )
      twod_mesh => mesh_collection%get_mesh( mesh, TWOD )
      fs => function_space_collection%get_fs( twod_mesh, 0, 0, W3 )
      call single_layer_function_space_chain%add( fs )
    end do

    nullify(mesh, twod_mesh, fs)

    call log_event( 'Function space chains created', LOG_LEVEL_INFO )

  end subroutine init_function_space_chains

  !> @brief  Finalises the function_space_collection.
  subroutine final_fem()

    implicit none

    call final_chi_transforms()

  end subroutine final_fem

end module driver_fem_mod
