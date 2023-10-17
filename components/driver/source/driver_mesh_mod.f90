!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief    Set up and destroy partitioned 3D mesh(es).
!> @details  Contains routines to:
!!            i) Read global UGRID meshes and set up partitioned 3D mesh(es),
!!           ii) Destroy partitioned mesh(es).
module driver_mesh_mod

  use constants_mod,              only: i_def, l_def, r_def, str_def, imdi, &
                                        str_max_filename
  use extrusion_mod,              only: extrusion_type,              &
                                        uniform_extrusion_type,      &
                                        geometric_extrusion_type,    &
                                        quadratic_extrusion_type,    &
                                        shifted_extrusion_type,      &
                                        double_level_extrusion_type, &
                                        PRIME_EXTRUSION
  use global_mesh_collection_mod, only: global_mesh_collection, &
                                        global_mesh_collection_type
  use global_mesh_mod,            only: global_mesh_type
  use local_mesh_collection_mod,  only: local_mesh_collection, &
                                        local_mesh_collection_type
  use local_mesh_mod,             only: local_mesh_type
  use log_mod,                    only: log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_INFO,    &
                                        LOG_LEVEL_ERROR
  use mesh_collection_mod,        only: mesh_collection, mesh_collection_type
  use mesh_mod,                   only: mesh_type
  use ncdf_quad_mod,              only: ncdf_quad_type
  use partition_mod,              only: partition_type, &
                                        partitioner_interface

  use ugrid_2d_mod,               only: ugrid_2d_type
  use ugrid_mesh_data_mod,        only: ugrid_mesh_data_type
  use ugrid_file_mod,             only: ugrid_file_type


  use base_mesh_config_mod,       only: file_prefix,             &
                                        prime_mesh_name,         &
                                        prepartitioned,          &
                                        key_from_geometry,       &
                                        key_from_topology,       &
                                        geometry,                &
                                        geometry_spherical,      &
                                        geometry_planar,         &
                                        topology,                &
                                        topology_fully_periodic, &
                                        topology_non_periodic
  use extrusion_config_mod,       only: method,           &
                                        method_uniform,   &
                                        method_geometric, &
                                        method_quadratic, &
                                        domain_top, number_of_layers
  use planet_config_mod,          only: scaled_radius

  implicit none

  private
  public  :: init_mesh

contains

!> @brief  Generates a mesh and determines the basis functions and dofmaps.
!>
!> @param[in]   local_rank                Number of the MPI rank of this process
!> @param[in]   total_ranks               Total number of MPI ranks in this job
!> @param[in]   base_mesh_names           Array of base mesh names. For each of
!!                                        these, a 3D "prime" extrusion will be
!!                                        constructed as well as its 2D version
!> @param[in]   shifted_mesh_names        Optional, array of names of meshes for
!!                                        which to make a shifted extrusion
!> @param[in]   double_level_mesh_names   Optional, array of names of meshes for
!!                                        which to make a double level extrusion
!> @param[in]   required_stencil_depth    Optional, stencil depth that local
!!                                        meshes should support
!> @param[in]   input_extrusion           Optional, mesh extrusion object
subroutine init_mesh( local_rank, total_ranks,        &
                      base_mesh_names,                &
                      shifted_mesh_names,             &
                      double_level_mesh_names,        &
                      required_stencil_depth,         &
                      input_extrusion )

  use finite_element_config_mod,  only: cellshape,          &
                                        key_from_cellshape, &
                                        cellshape_triangle, &
                                        cellshape_quadrilateral

  implicit none

  integer(kind=i_def),    intent(in) :: local_rank
  integer(kind=i_def),    intent(in) :: total_ranks
  character(len=str_def), intent(in) :: base_mesh_names(:)

  ! Optional lists of meshes to create shifted or double level extrusions of
  character(len=str_def), optional, allocatable, intent(in) :: shifted_mesh_names(:)
  character(len=str_def), optional, allocatable, intent(in) :: double_level_mesh_names(:)

  integer(kind=i_def),    optional, intent(in) :: required_stencil_depth

  ! Optional prime extrusion can be passed in for model-specific cases
  class(extrusion_type),  optional, intent(in) :: input_extrusion

  ! Parameters
  integer(kind=i_def), parameter :: max_factor_iters = 10000

  ! Local variables
  integer(kind=i_def) :: xproc  ! Processor ranks in mesh panel x-direction
  integer(kind=i_def) :: yproc  ! Processor ranks in mesh panel y-direction

  procedure(partitioner_interface), pointer :: partitioner_ptr => null()

  logical(kind=l_def) :: create_shifted_meshes      = .false.
  logical(kind=l_def) :: create_double_level_meshes = .false.

  integer(kind=i_def) :: i, j, k, stencil_depth

  character(len=str_def) :: mesh_name_A
  character(len=str_def) :: mesh_name_B
  character(len=str_def) :: twod_name_A
  character(len=str_def) :: twod_name_B

  class(extrusion_type), allocatable :: prime_extrusion

  character(str_max_filename)     :: input_mesh_file
  character(str_def), allocatable :: tmp_mesh_names(:)
  character(str_def), allocatable :: target_mesh_names(:)
  character(str_def), allocatable :: shifted_meshes_to_create(:)
  character(str_def), allocatable :: double_level_meshes_to_create(:)

  character(len=9), parameter :: routine_name = 'init_mesh'

  type(mesh_type),       pointer :: mesh => null()
  type(local_mesh_type), pointer :: local_mesh => null()

  ! Set up stencil depth
  if (present(required_stencil_depth)) then
    stencil_depth = required_stencil_depth
  else
    stencil_depth = 1
  end if

  ! Sort out prime mesh extrusion
  if (present(input_extrusion)) then
    allocate( prime_extrusion, source=input_extrusion )
  else
    allocate( prime_extrusion, source=create_prime_extrusion() )
  end if

  ! Currently only quad elements are fully functional
  if (cellshape /= cellshape_quadrilateral) then
    call log_event( "Reference_element must be QUAD for now...", &
                    LOG_LEVEL_ERROR )
  end if

  !=================================================================
  ! 1.0 Use input args to determine which meshes to create
  !=================================================================
  if (present(shifted_mesh_names)) then
    if (allocated(shifted_mesh_names)) then
      create_shifted_meshes = .true.
      allocate(shifted_meshes_to_create(SIZE(shifted_mesh_names)))
      shifted_meshes_to_create(:) = shifted_mesh_names(:)
    end if
  end if
  if (present(double_level_mesh_names)) then
    if (allocated(double_level_mesh_names)) then
      create_double_level_meshes = .true.
      allocate(double_level_meshes_to_create(SIZE(double_level_mesh_names)))
      double_level_meshes_to_create(:) = double_level_mesh_names(:)
    end if
  end if

  if (prepartitioned) then

    !=================================================================
    ! 2.0 Read in local meshes / partition information / mesh maps
    !     direct from file.
    !=================================================================
    !
    ! For this local rank, a mesh input file with a common base name
    ! with the following form should exist.
    !
    !   <input_basename>_<local_rank>_<total_ranks>.nc
    !
    ! Where 1 rank is assigned to each mesh partition.
    write(input_mesh_file,'(A,2(I0,A))') &
        trim(file_prefix) // '_', local_rank, '-', total_ranks, '.nc'

    call log_event( 'Using pre-partitioned mesh file:', LOG_LEVEL_INFO )
    call log_event( '   '//trim(input_mesh_file), LOG_LEVEL_INFO )
    call log_event( "Loading local mesh(es)", LOG_LEVEL_INFO )

    ! 2.1 Read in all local meshes for this rank and
    !     initialise local meshes from them.
    !=================================================================
    ! Each partitioned mesh file will contain meshes of the
    ! same name as all other partitions.
    call load_all_local_meshes( input_mesh_file, base_mesh_names )



    ! 2.2 Check loaded local meshes have the required stencil depth.
    !=================================================================
    call check_stencil_depths( stencil_depth )


    ! 2.3 Assign load mesh maps.
    !=================================================================
    ! Mesh map identifiers are determined by the source/target mesh
    ! IDs they relate to. As a result integrid mesh maps need to be
    ! loaded after all local meshes have been loaded.
    tmp_mesh_names = local_mesh_collection%get_mesh_names()

    do i=1, size(tmp_mesh_names)
      call load_local_mesh_maps( input_mesh_file, tmp_mesh_names(i) )
    end do

    if (allocated(tmp_mesh_names)) deallocate(tmp_mesh_names)

  else

    write(input_mesh_file,'(A)') trim(file_prefix) // '.nc'

    !=================================================================
    ! 3.0 Load and partition from global 2D meshes.
    !=================================================================
    call log_event( "Setting up partition mesh(es)", LOG_LEVEL_INFO )

    ! 3.1 Set constants that will control partitioning.
    !=================================================================
    call set_partition_parameters( total_ranks,       &
                                   xproc, yproc,      &
                                   partitioner_ptr )

    ! 3.2 Read in all global meshes and create local meshes from them.
    !=================================================================
    call create_all_base_meshes( input_mesh_file,                 &
                                 local_rank, total_ranks,         &
                                 xproc, yproc,                    &
                                 stencil_depth,                   &
                                 partitioner_ptr,                 &
                                 base_mesh_names )

    ! 3.3 Check loaded local meshes have the required stencil depth.
    !=================================================================
    call check_stencil_depths( stencil_depth )

    ! 3.4 Read in the global intergrid mesh mappings, then create the
    !     associated local mesh maps
    !=================================================================
    call create_mesh_maps(input_mesh_file)

  end if

  !=================================================================
  ! 4.0 Extrude all meshes into 3D local meshes
  !=================================================================
  call create_all_3d_meshes( prime_extrusion,              &
                             base_mesh_names,              &
                             create_shifted_meshes,        &
                             create_double_level_meshes,   &
                             shifted_meshes_to_create,     &
                             double_level_meshes_to_create )

  !=================================================================
  ! 5.0 Assign maps to local 3D meshes
  !=================================================================
  ! Add mesh maps to prime extrusion and 2D meshes
  if (SIZE(base_mesh_names) > 1) then
    do i = 1, SIZE(base_mesh_names)
      mesh_name_A = base_mesh_names(i)
      twod_name_A = trim(mesh_name_A)//'_2d'
      ! Find all of the target meshes associated with this mesh
      mesh => mesh_collection%get_mesh(mesh_name_A)
      local_mesh => mesh%get_local_mesh()
      call local_mesh%get_target_mesh_names(target_mesh_names)

      if (allocated(target_mesh_names)) then
        do j = 1, SIZE(target_mesh_names)
          mesh_name_B = target_mesh_names(j)
          call add_mesh_maps( mesh_name_A, mesh_name_B )
          twod_name_B = trim(mesh_name_B)//'_2d'
          call add_mesh_maps( twod_name_A, twod_name_B )
        end do
        deallocate(target_mesh_names)
      end if
    end do
  end if

  if (create_shifted_meshes) then
    if (SIZE(shifted_mesh_names) > 1) then
      do i = 1, SIZE(shifted_mesh_names)
        mesh_name_A = trim(shifted_mesh_names(i))//'_shifted'
        ! Find all of the target meshes associated with this mesh
        mesh => mesh_collection%get_mesh(mesh_name_A)
        local_mesh => mesh%get_local_mesh()
        call local_mesh%get_target_mesh_names(target_mesh_names)

        if (allocated(target_mesh_names)) then
          do j = 1, SIZE(target_mesh_names)
            mesh_name_B = trim(target_mesh_names(j))//'_shifted'
            do k = 1, SIZE(shifted_mesh_names)
              if ( mesh_name_B == trim(shifted_mesh_names(k))//'_shifted' ) then
                call add_mesh_maps( mesh_name_A, mesh_name_B )
              end if
            end do
          end do
          deallocate(target_mesh_names)
        end if
      end do
    end if
  end if

  if (create_double_level_meshes) then
    if (SIZE(double_level_mesh_names) > 1) then
      do i = 1, SIZE(double_level_mesh_names)
        mesh_name_A = trim(double_level_mesh_names(i))//'_double'
        ! Find all of the target meshes associated with this mesh
        mesh => mesh_collection%get_mesh(mesh_name_A)
        local_mesh => mesh%get_local_mesh()
        call local_mesh%get_target_mesh_names(target_mesh_names)

        if (allocated(target_mesh_names)) then
          do j = 1, SIZE(target_mesh_names)
            mesh_name_B = trim(target_mesh_names(j))//'_double'
            do k = 1, SIZE(double_level_mesh_names)
              if ( mesh_name_B == trim(double_level_mesh_names(k))//'_double' ) then
                call add_mesh_maps( mesh_name_A, mesh_name_B )
              end if
            end do
          end do
          deallocate(target_mesh_names)
        end if
      end do
    end if
  end if

  nullify(mesh, local_mesh)

  !=================================================================
  ! 6.0 Discard global mesh collection and all meshes in it.
  !     (They should not be used past this point in the code)
  !=================================================================
  if (allocated(prime_extrusion)) deallocate(prime_extrusion)

end subroutine init_mesh


!> @brief Sets common partition parameters to be applied to global meshes.
!>
!> @param[in]   total_ranks      Total number of MPI ranks in this job
!> @param[out]  xproc            Number of ranks in mesh panel x-direction
!> @param[out]  yproc            Number of ranks in mesh panel y-direction
!> @param[out]  partitioner_ptr  Mesh partitioning strategy
subroutine set_partition_parameters( total_ranks,       &
                                     xproc, yproc,      &
                                     partitioner_ptr )

  use partitioning_config_mod,    only: panel_decomposition,        &
                                        panel_xproc, panel_yproc,   &
                                        PANEL_DECOMPOSITION_AUTO,   &
                                        PANEL_DECOMPOSITION_ROW,    &
                                        PANEL_DECOMPOSITION_COLUMN, &
                                        PANEL_DECOMPOSITION_CUSTOM

  use partition_mod,              only: partitioner_cubedsphere_serial, &
                                        partitioner_cubedsphere,        &
                                        partitioner_planar

  implicit none

  integer(kind=i_def), intent(in)  :: total_ranks
  integer(kind=i_def), intent(out) :: xproc
  integer(kind=i_def), intent(out) :: yproc

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(kind=i_def) :: ranks_per_panel
  integer(kind=i_def) :: start_factor
  integer(kind=i_def) :: end_factor
  integer(kind=i_def) :: fact_count
  logical(kind=l_def) :: found_factors

  character(len=str_def) :: domain_desc

  integer(kind=i_def), parameter :: max_factor_iters = 10000

  partitioner_ptr => null()

  ! 1.0 Setup the partitioning strategy
  !===================================================================
  if (geometry == geometry_spherical  .and. &
      topology == topology_fully_periodic ) then

    ! Assume that we have a cubed sphere (and not a global lon-lat mesh)
    if (total_ranks == 1 .or. mod(total_ranks,6) == 0) then

      ranks_per_panel = total_ranks/6
      domain_desc = "6x"

      if (total_ranks == 1) then
        ! Serial run job
        ranks_per_panel = 1
        partitioner_ptr => partitioner_cubedsphere_serial
        call log_event( "Using serial cubed sphere partitioner", &
                        LOG_LEVEL_INFO )

      else
        ! Paralled run job
        partitioner_ptr => partitioner_cubedsphere
        call log_event( "Using parallel cubed sphere partitioner", &
                        LOG_LEVEL_INFO )
      end if

    else
      call log_event( "Total number of processors must be 1 (serial) "// &
                      "or a multiple of 6 for a cubed-sphere domain.",   &
                      LOG_LEVEL_ERROR )
    end if

  else ! Planar/LAM mesh

    ranks_per_panel = total_ranks
    domain_desc = ""

    partitioner_ptr => partitioner_planar
    call log_event( "Using planar mesh partitioner ", &
                    LOG_LEVEL_INFO )
  end if

  ! 2.0 Setup Panel decomposition
  !===================================================================
  select case(panel_decomposition)

  case( PANEL_DECOMPOSITION_AUTO )

    ! For automatic partitioning, try to partition into the squarest
    ! possible partitions by finding the two factors of ranks_per_panel
    ! that are closest to sqrt(ranks_per_panel). If two factors can't
    ! be found after max_factor_iters attempts, they would provide
    ! partitions that are too un-square, so an error is produced.
    start_factor  = nint(sqrt(real(ranks_per_panel, kind=r_def)), kind=i_def)
    end_factor    = max(1,(start_factor-max_factor_iters))
    found_factors = .false.
    do fact_count = start_factor, end_factor, -1
      if (mod(ranks_per_panel,fact_count) == 0) then
        found_factors = .true.
        exit
      end if
    end do

    if (found_factors) then
      xproc = fact_count
      yproc = ranks_per_panel/fact_count
    else
      call log_event( "Could not automatically partition domain.", &
                      LOG_LEVEL_ERROR )
    end if

  case( PANEL_DECOMPOSITION_ROW )
    xproc = ranks_per_panel
    yproc = 1

  case( PANEL_DECOMPOSITION_COLUMN )
    xproc = 1
    yproc = ranks_per_panel

  case( PANEL_DECOMPOSITION_CUSTOM )
    ! Use the values provided from the partitioning namelist
    xproc = panel_xproc
    yproc = panel_yproc

    if (xproc*yproc /= ranks_per_panel) then
      call log_event( "The values of panel_xproc and panel_yproc "// &
                      "are inconsistent with the total number of "// &
                      "processors available.", LOG_LEVEL_ERROR )
    end if

  case default

    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select

  if (total_ranks > 1) then
    write(log_scratch_space, '(2(A,I0))' ) &
        'Panel decomposition: ', xproc,'x', yproc
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end if

end subroutine set_partition_parameters

!> @brief  Reads in global meshes from UGRID file, partitions them
!!         and creates local meshes.
!>
!> @param[in]  input_mesh_file         Input file to load meshes from.
!> @param[in]  local_rank              Number of the local MPI rank
!> @param[in]  total_ranks             Total number of MPI ranks in this job
!> @param[in]  xproc                   Number of ranks in mesh panel x-direction
!> @param[in]  yproc                   Number of ranks in mesh panel y-direction
!> @param[in]  stencil_depth           Depth of cells outside the base cell
!!                                     of stencil
!> @param[in]  partitioner_ptr         Mesh partitioning strategy
!> @param[in]  base_mesh_names         Array of base mesh names to create
subroutine create_all_base_meshes( input_mesh_file,                 &
                                   local_rank, total_ranks,         &
                                   xproc, yproc,                    &
                                   stencil_depth,                   &
                                   partitioner_ptr,                 &
                                   base_mesh_names )

  use multigrid_config_mod, only: chain_mesh_tags, multigrid_chain_nitems

  implicit none

  character(len=str_max_filename),           intent(in) :: input_mesh_file
  integer(kind=i_def),                       intent(in) :: local_rank
  integer(kind=i_def),                       intent(in) :: total_ranks
  integer(kind=i_def),                       intent(in) :: xproc
  integer(kind=i_def),                       intent(in) :: yproc
  integer(kind=i_def),                       intent(in) :: stencil_depth
  procedure(partitioner_interface), pointer, intent(in) :: partitioner_ptr
  character(len=str_def),                    intent(in) :: base_mesh_names(:)

  integer(kind=i_def) :: n_multigrid_levels

  call create_base_meshes( input_mesh_file,         &
                           base_mesh_names,         &
                           local_rank, total_ranks, &
                           xproc, yproc,            &
                           stencil_depth,           &
                           partitioner_ptr )

  ! Check that the partitioning strategy will result in partitions that
  ! align for all levels of multigrid
  if ( SIZE(base_mesh_names) > 1 .and. multigrid_chain_nitems > 1 ) then
    n_multigrid_levels = SIZE(chain_mesh_tags)
    call check_multigrid_partitioning(prime_mesh_name, &
                                      xproc, yproc,    &
                                      n_multigrid_levels)
  end if

end subroutine create_all_base_meshes

!> @brief  Loads the given list of global meshes, partitions them
!!         and creates local meshes from them.
!>
!> @param[in]  input_mesh_file    Input file to load meshes from.
!> @param[in]  mesh_names[:]      Array of requested mesh names to load
!!                                from the mesh input file
!> @param[in]  local_rank         Number of the local MPI rank
!> @param[in]  total_ranks        Total number of MPI ranks in this job
!> @param[in]  xproc              Number of ranks in mesh panel x-direction
!> @param[in]  yproc              Number of ranks in mesh panel y-direction
!> @param[in]  stencil_depth      Depth of cells outside the base cell
!!                                of stencil.
!> @param[in]  partitioner_ptr    Mesh partitioning strategy
subroutine create_base_meshes( input_mesh_file,         &
                               mesh_names,              &
                               local_rank, total_ranks, &
                               xproc, yproc,            &
                               stencil_depth,           &
                               partitioner_ptr )

  implicit none

  character(len=str_max_filename), intent(in) :: input_mesh_file

  character(len=str_def), intent(in) :: mesh_names(:)
  integer(kind=i_def),    intent(in) :: local_rank
  integer(kind=i_def),    intent(in) :: total_ranks
  integer(kind=i_def),    intent(in) :: xproc
  integer(kind=i_def),    intent(in) :: yproc
  integer(kind=i_def),    intent(in) :: stencil_depth

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  type(ugrid_mesh_data_type)      :: ugrid_mesh_data
  type(global_mesh_type)          :: global_mesh
  type(global_mesh_type), pointer :: global_mesh_ptr
  type(partition_type)            :: partition
  type(local_mesh_type)           :: local_mesh
  integer(kind=i_def)             :: local_mesh_id, i
  logical(kind=l_def)             :: valid_geometry, valid_topology

  do i=1, size(mesh_names)
    if (.not. global_mesh_collection%check_for(mesh_names(i))) then

      ! Load mesh data into global_mesh
      call ugrid_mesh_data%read_from_file(trim(input_mesh_file), mesh_names(i))

      global_mesh = global_mesh_type( ugrid_mesh_data )
      call ugrid_mesh_data%clear()


      ! Check mesh has valid domain geometry
      valid_geometry = .false.
      select case(geometry)

      case(geometry_spherical)
        if ( global_mesh%is_geometry_spherical() ) valid_geometry = .true.

      case(geometry_planar)
        if ( global_mesh%is_geometry_planar() ) valid_geometry = .true.

      end select

      if ( .not. valid_geometry) then
        write(log_scratch_space, '(A)')        &
            'Mesh (' // trim(mesh_names(i)) // &
            ') in file is not valid as a '  // &
             trim(key_from_geometry(geometry))//' domain geometry'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR )
      end if


      ! Check mesh has valid topology
      valid_topology = .false.
      select case(topology)

      case(topology_fully_periodic)
        if ( global_mesh%is_topology_periodic() ) valid_topology = .true.

      case(topology_non_periodic)
        if ( global_mesh%is_topology_non_periodic() ) valid_topology = .true.

      end select

      if ( .not. valid_topology) then
        write(log_scratch_space, '(A)')           &
            'Mesh (' // trim(mesh_names(i)) //    &
            ') in file does not have a valid ' // &
             trim(key_from_topology(topology))//' topology'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR )
      end if

      call global_mesh_collection%add_new_global_mesh ( global_mesh )
      global_mesh_ptr => global_mesh_collection%get_global_mesh( mesh_names(i) )

      ! Create partition
      partition = partition_type( global_mesh_ptr, &
                                  partitioner_ptr, &
                                  xproc, yproc,    &
                                  stencil_depth,   &
                                  local_rank, total_ranks )
      ! Create local_mesh
      call local_mesh%initialise( global_mesh_ptr, partition )

      ! Make sure the local_mesh cell owner lookup is correct
      ! (Can only be done when the code is running on its full set of MPI tasks)
      call local_mesh%init_cell_owner()
      local_mesh_id = local_mesh_collection%add_new_local_mesh(local_mesh)

    end if

  end do

end subroutine create_base_meshes

!> @brief    Reads in and assigns available global intergrid maps from file.
!> @details  Global meshes which have been read into the model's global mesh
!!           collection will have a list of target mesh names. These target mesh
!!           names (if any) indicate the valid intergrid maps avaiable in the
!!           mesh file. This routine will read in the appropriate intergrid
!!           maps and assign them to the correct global mesh object.
!!
!> @param[in]  input_mesh_file  Input file to load mesh maps from.
subroutine create_mesh_maps( input_mesh_file )

  implicit none

  character(len=str_max_filename) :: input_mesh_file

  type(ncdf_quad_type) :: file_handler

  character(len=str_def), allocatable :: source_mesh_names(:)
  character(len=str_def), allocatable :: target_mesh_names(:)
  integer(kind=i_def),    allocatable :: gid_mesh_map(:,:,:)
  integer(kind=i_def),    allocatable :: lid_mesh_map(:,:,:)

  integer(kind=i_def) :: i, j, n, x, y
  integer(kind=i_def) :: n_meshes

  type(global_mesh_type), pointer :: source_global_mesh => null()

  type(local_mesh_type), pointer :: source_local_mesh => null()
  type(local_mesh_type), pointer :: target_local_mesh => null()

  integer(kind=i_def) :: ntarget_per_source_cell_x, ntarget_per_source_cell_y
  integer(kind=i_def) :: ncells
  integer(kind=i_def) :: target_local_mesh_id

  ! Read in the maps for each global mesh
  !=================================================================
  call file_handler%file_open(trim(input_mesh_file))

  source_mesh_names = global_mesh_collection%get_mesh_names()
  n_meshes = global_mesh_collection%n_meshes()

  ! Loop over every source mesh
  do i=1, n_meshes
    ! Get the global and local source mesh
    source_global_mesh => &
        global_mesh_collection%get_global_mesh( source_mesh_names(i) )
    source_local_mesh => &
        local_mesh_collection%get_local_mesh( source_mesh_names(i) )
    call source_global_mesh%get_target_mesh_names( target_mesh_names )
    if (allocated(target_mesh_names)) then
      ! Loop over each target mesh
      do j=1, size(target_mesh_names)
        target_local_mesh => &
           local_mesh_collection%get_local_mesh( target_mesh_names(j) )

        if ( associated(target_local_mesh) ) then
          ! Read in the global mesh map
          call file_handler%read_map( source_mesh_names(i), &
                                      target_mesh_names(j), &
                                      gid_mesh_map )

          ! Create the local mesh map
          ntarget_per_source_cell_x = size(gid_mesh_map, 1)
          ntarget_per_source_cell_y = size(gid_mesh_map, 2)
          ncells = source_local_mesh%get_num_cells_in_layer()
          allocate( lid_mesh_map( ntarget_per_source_cell_x, &
                                  ntarget_per_source_cell_y, &
                                  ncells ) )
          ! Convert global cell IDs in the global mesh map
          ! into local cell IDs in a local mesh map
          do x=1, ntarget_per_source_cell_x
            do y=1, ntarget_per_source_cell_y
              do n=1, ncells
                lid_mesh_map( x,y, n ) = target_local_mesh%get_lid_from_gid( &
                    gid_mesh_map( x,y, source_local_mesh%get_gid_from_lid(n) ) )
              end do
            end do
          end do

          ! Put the local mesh map in the local mesh
          target_local_mesh_id = target_local_mesh%get_id()
          call source_local_mesh%add_local_mesh_map( target_local_mesh_id, &
                                                     lid_mesh_map )

          if(allocated( gid_mesh_map )) deallocate( gid_mesh_map )
          if(allocated( lid_mesh_map )) deallocate( lid_mesh_map )
        end if

      end do
      if(allocated( target_mesh_names )) &
                                      deallocate( target_mesh_names )
    end if
  end do

  if(allocated( source_mesh_names ))  deallocate( source_mesh_names)
  call file_handler%file_close()

  return
end subroutine create_mesh_maps

!> @brief    Generates the 3D-meshes required by the model configuration.
!> @details  The extrusion types are set up for the required configuration
!!           before 3D-meshes are instantiated.
!>
!> @param[in]  extrusion                  The prime vertical mesh extrusion
!> @param[in]  base_mesh_names            Array of base mesh names. For each of
!!                                        these, a 3D "prime" extrusion will be
!!                                        constructed as well as its 2D version
!> @param[in]  create_shifted_meshes      Create shifted meshes
!> @param[in]  create_double_level_meshes Create double-level meshes
!> @param[in]  shifted_mesh_names         Optional, array of names of meshes for
!!                                        which to make a shifted extrusion
!> @param[in]  double_level_mesh_names    Optional, array of names of meshes for
!!                                        which to make a double level extrusion
subroutine create_all_3D_meshes( extrusion,                  &
                                 base_mesh_names,            &
                                 create_shifted_meshes,      &
                                 create_double_level_meshes, &
                                 shifted_mesh_names,         &
                                 double_level_mesh_names )

  use extrusion_mod,           only: TWOD, SHIFTED, DOUBLE_LEVEL
  use extrusion_config_mod,    only: domain_top

  implicit none

  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in) :: base_mesh_names(:)
  logical(kind=l_def),   intent(in) :: create_shifted_meshes
  logical(kind=l_def),   intent(in) :: create_double_level_meshes
  character(str_def),    intent(in) :: shifted_mesh_names(:)
  character(str_def),    intent(in) :: double_level_mesh_names(:)

  class(extrusion_type), allocatable :: extrusion_shifted
  class(extrusion_type), allocatable :: extrusion_double
  type(uniform_extrusion_type)       :: extrusion_2d

  character(len=str_def) :: mesh_name

  integer(i_def) :: i

  integer(i_def), parameter :: one_layer = 1_i_def
  real(r_def),    parameter :: atmos_bottom = 0.0_r_def

  ! 1.0 Prime extrusions and 2D meshes
  !===================================================================
  do i = 1, SIZE(base_mesh_names)
    write(log_scratch_space,'(A,I4)') 'Creating prime mesh ', i
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    call create_3d_mesh( base_mesh_names(i), extrusion )

  ! 2.0 Generate addition 3d mesh partitions based on the prime mesh
  ! NOTE: This includes 2D meshes as they are currently implemented
  !       as a 3D mesh of 1-layer thick.
    extrusion_2d = uniform_extrusion_type( atmos_bottom, &
                                           domain_top,   &
                                           one_layer,    &
                                           TWOD )

    mesh_name = trim(base_mesh_names(i))//'_2d'

    call create_3d_mesh( base_mesh_names(i), &
                         extrusion_2d,    &
                         mesh_name=mesh_name )
  end do

  ! 3.0 Generate shifted 3d mesh partition
  if (create_shifted_meshes) then

    if (allocated(extrusion_shifted)) deallocate(extrusion_shifted)
    allocate( extrusion_shifted, &
              source=shifted_extrusion_type(extrusion) )

    do i = 1, SIZE(shifted_mesh_names)
      mesh_name = trim(shifted_mesh_names(i))//'_shifted'

      call create_3d_mesh( shifted_mesh_names(i),   &
                           extrusion_shifted,       &
                           mesh_name=mesh_name )
    end do
  end if

  ! 4.0 Generate double-level 3d mesh partition
  if (create_double_level_meshes) then

    allocate( extrusion_double, &
              source=double_level_extrusion_type(extrusion) )

    do i = 1, SIZE(double_level_mesh_names)
      mesh_name = trim(double_level_mesh_names(i))//'_double'

      call create_3d_mesh( double_level_mesh_names(i),  &
                           extrusion_double,            &
                           mesh_name=mesh_name )
    end do
  end if

  if (allocated(extrusion_shifted)) deallocate( extrusion_shifted )
  if (allocated(extrusion_double))  deallocate( extrusion_double )

  return
end subroutine create_all_3D_meshes

!> @brief    Generates a single (partitioned) 3D-mesh.
!> @details  Instantiates a 3D-mesh partition and adds it to the model's
!!           mesh collection. Multiple meshes may be generated in the model
!!           based on the same global mesh but with differing extrusions.
!>
!> @param[in]  base_mesh_name  Name of base 2D-mesh
!> @param[in]  extrusion       Extrusion type for this 3D-mesh
!> @param[in]  mesh_name       Optional, Name of local 3D-mesh,
!>                             defaults to base_mesh_name
subroutine create_3d_mesh( base_mesh_name, &
                           extrusion,      &
                           mesh_name )

  use mesh_collection_mod, only: mesh_collection

  implicit none

  character(len=str_def),    intent(in) :: base_mesh_name
  class(extrusion_type),     intent(in) :: extrusion
  character(len=str_def),    intent(in), &
                             optional   :: mesh_name

  type(local_mesh_type), pointer :: local_mesh_ptr => null()

  type(mesh_type)        :: mesh
  integer(kind=i_def)    :: mesh_id
  character(len=str_def) :: name

  if (.not. present(mesh_name)) then
    name = base_mesh_name
  else
    name = mesh_name
  end if

  ! 1.0 Check 3D-mesh hasn't already been created
  !===============================================
  if ( mesh_collection%check_for(name) ) return

  ! 2.0 Create the 3D-mesh
  !===============================================
  local_mesh_ptr => local_mesh_collection%get_local_mesh(base_mesh_name)

  mesh = mesh_type( local_mesh_ptr, extrusion, mesh_name=name )

  mesh_id = mesh_collection%add_new_mesh( mesh )
  call mesh%clear()

  ! 3.0 Report on mesh creation
  !===============================================
  write(log_scratch_space,'(A,I0,A)')                 &
      '   ... "'//trim(name)//'"(id:', mesh_id,') '// &
      'based on mesh "'//trim(base_mesh_name)//'"'

  if (mesh_id /= imdi) then
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else
    write(log_scratch_space,'(A,I0,A)') &
        trim(log_scratch_space)//' (FAILED)'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  return
end subroutine create_3d_mesh


!> @brief  Creates the prime vertical mesh extrusion.
!> @return  Resulting extrusion object
function create_prime_extrusion() result(new)

  implicit none

  class(extrusion_type), allocatable :: new

  real(kind=r_def) :: domain_bottom

  if (allocated(new)) deallocate(new)

  select case (geometry)
    case (geometry_planar)
      domain_bottom = 0.0_r_def
    case (geometry_spherical)
      domain_bottom = scaled_radius
    case default
      call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
  end select

  select case (method)
    case (method_uniform)
      allocate( new, source=uniform_extrusion_type( domain_bottom,    &
                                                    domain_top,       &
                                                    number_of_layers, &
                                                    PRIME_EXTRUSION ) )
    case (method_quadratic)
      allocate( new, source=quadratic_extrusion_type( domain_bottom,    &
                                                      domain_top,       &
                                                      number_of_layers, &
                                                      PRIME_EXTRUSION ) )
    case (method_geometric)
      allocate( new, source=geometric_extrusion_type( domain_bottom,    &
                                                      domain_top,       &
                                                      number_of_layers, &
                                                      PRIME_EXTRUSION ) )
    case default
      call log_event("Invalid method for simple extrusion", LOG_LEVEL_ERROR)
  end select

end function create_prime_extrusion


!> @brief    Assigns intergrid maps to 3D-mesh partitions.
!> @details  Adds local ID intergrid mappings 3D-mesh partitions.
!!           Local intergrid maps are assigned to source (source -> target)
!!           and target (target -> source) meshes. A number of
!!           of assumptions are made.
!!
!!           *  The base global meshes of source/target 3D-mesh
!!              are present in the global_mesh_collection.
!!           *  Intergrid maps for the base global meshes were read in
!!              and assigned.
!!           *  Partitioned 3D-mesh is present in the
!!              mesh_collection.
!>
!> @param[in]  source_mesh_name  Name of source 3D mesh partition
!> @param[in]  target_mesh_name  Name of target 3D mesh partition
subroutine add_mesh_maps( source_mesh_name, &
                          target_mesh_name )

  use mesh_collection_mod, only: mesh_collection

  implicit none

  character(len=str_def), intent(in) :: source_mesh_name
  character(len=str_def), intent(in) :: target_mesh_name

  type(mesh_type), pointer :: source_mesh => null()
  type(mesh_type), pointer :: target_mesh => null()


  ! Now add in any mesh maps required by multigrid or
  ! the multires_coupling miniapp
  source_mesh => mesh_collection % get_mesh( source_mesh_name )
  target_mesh => mesh_collection % get_mesh( target_mesh_name )

  if ( associated(source_mesh) .and. &
       associated(target_mesh) ) then

    ! Mesh tag names may be different but "could point to the same mesh
    ! So check the IDs are not the same
    if (source_mesh%get_id() == target_mesh%get_id()) then
      write(log_scratch_space,'(A)')                  &
          'Unable to create intergrid map: Source('// &
          trim(source_mesh_name)//' and target('//    &
          trim(target_mesh_name)//') mesh IDs are the same'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    call source_mesh % add_mesh_map (target_mesh)
    call target_mesh % add_mesh_map (source_mesh)
    write(log_scratch_space,'(A,I0,A)')     &
        'Adding intergrid map "'//          &
         trim(source_mesh_name)//'"<-->"'// &
         trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_INFO)
  else
    write(log_scratch_space,'(A,I0,A)')          &
        'Unable to create mesh map between "'//  &
        trim(source_mesh_name)//'"-"'//          &
        trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR)
  end if

  nullify(source_mesh)
  nullify(target_mesh)

  return
end subroutine add_mesh_maps


!> @brief    Load all the local mesh(es) required.
!> @details  Load the local mesh(es) UGRID mesh input file.
!>
!> @param[in]  input_mesh_file            UGRID file containing data to
!>                                        populate <local_mesh_type> objects.
!> @param[in]  base_mesh_names            Array of base mesh names to read in
subroutine load_all_local_meshes( input_mesh_file, base_mesh_names )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: base_mesh_names(:)

  integer(i_def) :: i

  do i = 1, SIZE(base_mesh_names)
    write(log_scratch_space,'(A)') &
        'Reading mesh: "'//trim(base_mesh_names(i))//'"'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    call load_local_mesh(input_mesh_file, base_mesh_names(i))
  end do

end subroutine load_all_local_meshes


!> @brief    Loads a specified local mesh from a UGRID file.
!> @details  Loads a single local mesh, specified by name from the
!>           UGRID mesh input file.
!>
!> @param[in]  input_mesh_file    UGRID file containing data to
!>                                populate <local_mesh_type> object.
!> @param[in]  mesh_name          The name of the local mesh topology to load
!>                                from the <input_mesh_file>.
subroutine load_local_mesh( input_mesh_file, mesh_name )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: mesh_name

  type(ugrid_mesh_data_type) :: ugrid_mesh_data
  type(local_mesh_type)      :: local_mesh

  integer(i_def) :: local_mesh_id

  if (.not. local_mesh_collection%check_for(mesh_name)) then

    ! Load mesh data into local_mesh

    call ugrid_mesh_data%read_from_file( input_mesh_file, &
                                         mesh_name )
    call local_mesh%initialise_from_ugrid_data( ugrid_mesh_data )

    ! Assign cell ownership for lookup.
    call local_mesh%init_cell_owner()
    local_mesh_id = local_mesh_collection%add_new_local_mesh( local_mesh )

  end if

end subroutine load_local_mesh


!> @brief    Loads local intergrid cell maps from a UGRID file.
!> @details  For a given mesh, (that has been previosly loaded), the
!>           the intergrid cell maps for target local mesh(es) are
!>           loaded the <input_mesh_file> and assigned to the specfied
!>           local mesh.
!>
!> @param[in]  input_mesh_file    UGRID file containing data to
!>                                populate <local_mesh_type> object.
!> @param[in]  mesh_name          The name of the local mesh topology to load
!>                                from the <input_mesh_file>.
subroutine load_local_mesh_maps( input_mesh_file, mesh_name )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: mesh_name

  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: lid_mesh_map(:,:,:)

  integer(i_def) :: i

  type(ncdf_quad_type) :: file_handler


  type(local_mesh_type), pointer :: source_mesh => null()
  type(local_mesh_type), pointer :: target_mesh => null()

  integer(i_def) :: target_mesh_id


  ! Read in the maps for each local mesh
  !=================================================================
  source_mesh => local_mesh_collection%get_local_mesh( mesh_name )
  if (.not. associated(source_mesh)) then
    write(log_scratch_space,'(A)') ' Mesh "'//trim(mesh_name)// &
                                   '" not found in collection'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  call source_mesh%get_target_mesh_names(target_mesh_names)

  if (allocated(target_mesh_names)) then

    call file_handler%file_open(trim(input_mesh_file))

    do i=1, size(target_mesh_names)
      if ( local_mesh_collection%check_for( target_mesh_names(i) ) ) then

        ! Read in the local mesh map.
        call file_handler%read_map( mesh_name,            &
                                    target_mesh_names(i), &
                                    lid_mesh_map )

        target_mesh => local_mesh_collection%get_local_mesh(target_mesh_names(i))
        target_mesh_id =  target_mesh%get_id()

        ! Assign local mesh map to the local mesh object.
        call source_mesh%add_local_mesh_map( target_mesh_id, &
                                             lid_mesh_map )

        if (allocated( lid_mesh_map )) deallocate( lid_mesh_map )

      end if

    end do

    call file_handler%file_close()

    if (allocated( target_mesh_names )) deallocate( target_mesh_names )
  end if

end subroutine load_local_mesh_maps


!> @brief  Checks that the local meshes generated meet the stencil depth
!>         requirements.
!> @param[in]  stencil_depth  The required stencil depth that the local
!>                            meshes should support.
!==========================================================================
subroutine check_stencil_depths( stencil_depth )

  implicit none

  integer(i_def), intent(in) :: stencil_depth

  character(str_def), allocatable :: mesh_names(:)
  type(local_mesh_type), pointer  :: local_mesh_ptr => null()

  integer(i_def) :: max_stencil_depth
  integer(i_def) :: i

  mesh_names = local_mesh_collection%get_mesh_names()

  do i=1, size(mesh_names)

    local_mesh_ptr => local_mesh_collection%get_local_mesh(mesh_names(i))
    max_stencil_depth = local_mesh_ptr%get_max_stencil_depth()

    if ( max_stencil_depth < stencil_depth ) then

      write(log_scratch_space,'(2(A,I0),A)')                      &
         'Insufficient stencil depth, [', max_stencil_depth, '<', &
         stencil_depth, '] for (trim(mesh_names(i))) local mesh.'

      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end do

  local_mesh_ptr => null()
  if ( allocated(mesh_names) ) deallocate(mesh_names)

end subroutine check_stencil_depths

!> @brief  Checks that the partitioning strategy will work on the
!>         lowest resolution multigrid level
!> @param[in]  mesh_name          The name of the mesh on which the
!>                                multigrid is being applied
!> @param[in]  xproc              The number of partitions across a
!>                                panel in the mesh in the x-direction
!> @param[in]  yproc              The number of partitions across a
!>                                panel in the mesh in the y-direction
!> @param[in]  n_multigrid_levels The number of levels of multigrid
!>                                being applied to the mesh
!==========================================================================
subroutine check_multigrid_partitioning( mesh_name, &
                                         xproc, yproc,    &
                                         n_multigrid_levels )
  use reference_element_mod, only : W, E

  implicit none

  character(str_def), intent(in) :: mesh_name
  integer(i_def), intent(in) :: xproc, yproc
  integer(i_def), intent(in) :: n_multigrid_levels

  type(global_mesh_type), pointer :: global_mesh
  integer(i_def) :: nn  ! "C" number of input mesh (as in Cnn)
  integer(i_def) :: xx  ! number of cells across a panel of lowest res multigrid mesh in x-direction
  integer(i_def) :: yy  ! number of cells across a panel of lowest res multigrid mesh in y-direction
  integer(i_def) :: void_cell    ! Cell id that marks the cell as a cell outside of the partition.
  integer(i_def) :: w_cell       ! The id of a cell on the western edge of the domain
  integer(i_def) :: cell_next(4) ! The cells around the cell being queried
  integer(i_def) :: cell_next_e  ! The cell to the east of the cell being queried
  integer(i_def) :: num_cells_x  ! number of cells across a panel of the input mesh in x-direction
  integer(i_def) :: num_cells_y  ! number of cells across a panel of the input mesh in y-direction
  logical(l_def) :: periodic_xy(2) ! Is mesh periodic in the x/y-axes

  ! Get the global mesh object of the input mesh
  global_mesh => global_mesh_collection%get_global_mesh(mesh_name)

  if (geometry == geometry_spherical  .and. &
      topology == topology_fully_periodic ) then  ! Panelled globe mesh such as cubedsphere

    ! Generate the "C" number (as in Cnn) from global mesh data
    ! Assume the panels on the globe are square
    nn = int(sqrt(real(global_mesh%get_ncells()/global_mesh%get_npanels())))

    ! Check that the lowest resolution mesh will line up with the partitions
    ! (xx is the number of cells across a panel of the lowest resolution
    ! multigrid mesh)
    xx = nn / 2**(n_multigrid_levels-1)
    if (mod(xx,xproc) /= 0 .or. mod(xx,yproc) /= 0) then
      write(log_scratch_space, "(a,i0,a,i0,a)") &
          'The level ',n_multigrid_levels, ' multigrid of a C',nn, &
          ' mesh does not align with the partitioning strategy'
      call log_event(log_scratch_space,LOG_LEVEL_ERROR )
    end if
  else  ! Planar or LAM mesh
    ! Planar or LAM mesh mesh might be non-square - so find the dimensions
    ! First find a cell on the west edge of the domain
    ! If periodic, cell id 1 can be used as mesh conectivity loops round
    w_cell = 1

    void_cell = global_mesh%get_void_cell()
    periodic_xy = global_mesh%get_mesh_periodicity()
    if ( .not. periodic_xy(1) ) then
      ! If not periodic in E-W direction then walk West until you reach mesh
      ! edge defined by the void cell.
      call global_mesh%get_cell_next(w_cell,cell_next)
      do while (cell_next(W) /= void_cell)
        w_cell = cell_next(W)
        call global_mesh%get_cell_next(w_cell,cell_next)
      end do
    end if

    ! Work out number of cells in x and y directions
    num_cells_x = 1
    ! Starting at the West edge of the mesh, walk East until you reach either
    ! the cell you started at (periodic) or a void cell (LAM)
    ! - this determines the number of cells in the x direction
    call global_mesh%get_cell_next(w_cell,cell_next)
    cell_next_e = cell_next(E)
    do while (cell_next_e /= w_cell .and. cell_next_e /= void_cell)
      num_cells_x=num_cells_x+1
      call global_mesh%get_cell_next(cell_next_e, cell_next)
      cell_next_e = cell_next(E)
    end do
    ! Infer num_cells_y from the total domin size and num_cells_x
    num_cells_y=global_mesh%get_ncells()/num_cells_x

    xx = num_cells_x / 2**(n_multigrid_levels-1)
    yy = num_cells_y / 2**(n_multigrid_levels-1)
    if (mod(xx,xproc) /= 0 .or. mod(yy,yproc) /= 0) then
      write(log_scratch_space,"(a,i0,a)") &
       'The level ',n_multigrid_levels, &
       ' multigrid of the mesh does not align with the partitioning strategy'
      call log_event(log_scratch_space,LOG_LEVEL_ERROR )
    end if
  end if

end subroutine check_multigrid_partitioning

end module driver_mesh_mod
