!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Set up and destroy partitioned 3D mesh(es)
!> @details Contains routines to:
!>           i) Read global ugrid meshes and set up partitioned 3D mesh(es)
!>          ii) Destroy partitioned mesh(es)
module create_mesh_mod

  use constants_mod,              only: i_def, str_def, l_def, &
                                        r_def, imdi, i_native
  use extrusion_mod,              only: extrusion_type, &
                                        uniform_extrusion_type
  use global_mesh_collection_mod, only: global_mesh_collection, &
                                        global_mesh_collection_type
  use global_mesh_mod,            only: global_mesh_type
  use local_mesh_collection_mod,  only: local_mesh_collection
  use local_mesh_mod,             only: local_mesh_type
  use log_mod,                    only: log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_INFO,    &
                                        LOG_LEVEL_ERROR
  use mesh_collection_mod,        only: mesh_collection
  use mesh_mod,                   only: mesh_type
  use ncdf_quad_mod,              only: ncdf_quad_type
  use partition_mod,              only: partition_type, &
                                        partitioner_interface
  use ugrid_2d_mod,               only: ugrid_2d_type
  use ugrid_mesh_data_mod,        only: ugrid_mesh_data_type
  use ugrid_file_mod,             only: ugrid_file_type

  use base_mesh_config_mod,       only: filename, prime_mesh_name, &
                                        offline_partitioning,      &
                                        key_from_geometry,         &
                                        key_from_topology,         &
                                        geometry,                  &
                                        geometry_spherical,        &
                                        geometry_planar,           &
                                        topology,                  &
                                        topology_fully_periodic,   &
                                        topology_non_periodic

  implicit none

  private
  public  :: init_mesh, final_mesh

  private :: set_partition_parameters,     &
             create_all_base_meshes,       &
             create_base_meshes,           &
             create_mesh_maps,             &
             create_all_3d_meshes,         &
             create_3d_mesh,               &
             add_mesh_maps

contains


!===============================================================================
!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!>
!> @param[in]  local_rank   Number of the MPI rank of this process
!> @param[in]  total_ranks  Total number of MPI ranks in this job
!===============================================================================
!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!> @param[in]  local_rank                    Number of the MPI rank of this process
!> @param[in]  total_ranks                   Total number of MPI ranks in this job
!> @param[in]  stencil_depth                 Stencil depth that local meshes should support
!> @param[out] mesh_id                       Mesh id of partitioned prime mesh
!> @param[out] twod_mesh_id                  Mesh id of the 2D (surface) mesh
!> @param[out] shifted_mesh_id               Mesh id of vertically shifted mesh
!>                                           with an extra level
!> @param[out] double_level_mesh_id          Mesh id of vertically double level mesh
!> @param[out] multigrid_mesh_ids            Gungho multigrid chain mesh ids
!> @param[out] multigrid_2D_mesh_ids         Gungho multigrid chain 2D_mesh ids
!> @param[in]  use_multigrid                 Optional, Configuration switch for multigrid
!> @param[out] multires_coupling_mesh_ids    Multiresolution Coupling miniapp mesh ids
!> @param[out] multires_coupling_2D_mesh_ids Multiresolution Coupling miniapp 2D_mesh ids
!> @param[in]  multires_coupling_mesh_tags   Multiresolution Coupling miniapp mesh names
!> @param[in]  use_multires_coupling         Logical flag to enable multiresolution atmospheric coupling

subroutine init_mesh( local_rank, total_ranks,       &
                      stencil_depth, mesh_id,        &
                      twod_mesh_id,                  &
                      shifted_mesh_id,               &
                      double_level_mesh_id,          &
                      multigrid_mesh_ids,            &
                      multigrid_2D_mesh_ids,         &
                      use_multigrid,                 &
                      multires_coupling_mesh_ids,    &
                      multires_coupling_2D_mesh_ids, &
                      multires_coupling_mesh_tags,   &
                      use_multires_coupling )

  use finite_element_config_mod,  only: cellshape,          &
                                        key_from_cellshape, &
                                        cellshape_triangle, &
                                        cellshape_quadrilateral
  use multigrid_config_mod,       only: chain_mesh_tags

  implicit none

  integer(i_def), intent(in)  :: local_rank
  integer(i_def), intent(in)  :: total_ranks

  ! stencil_depth is the requested depth (of cells outside the cell over
  ! which the stencil is based) of the stencil to be used on fields with
  ! this partition.
  !
  ! A single cell stencil will, therefore, have a stencil_depth=0.
  ! A nine-point square region stencil will have stencil_depth=1
  integer(i_def), intent(in)  :: stencil_depth

  integer(i_def), intent(out) :: mesh_id

  integer(i_def), intent(out), optional :: twod_mesh_id
  integer(i_def), intent(out), optional :: shifted_mesh_id
  integer(i_def), intent(out), optional :: double_level_mesh_id
  integer(i_def), intent(out), optional, allocatable :: multigrid_mesh_ids(:)
  integer(i_def), intent(out), optional, allocatable :: multigrid_2d_mesh_ids(:)
  integer(i_def), intent(out), optional, allocatable :: multires_coupling_mesh_ids(:)
  integer(i_def), intent(out), optional, allocatable :: multires_coupling_2d_mesh_ids(:)

  character(str_def), intent(in), optional :: multires_coupling_mesh_tags(:)

  logical(l_def), intent(in), optional :: use_multigrid
  logical(l_def), intent(in), optional :: use_multires_coupling

  ! Parameters
  integer(i_def), parameter :: max_factor_iters = 10000

  ! Local variables
  integer(i_def) :: xproc  ! Processor ranks in mesh panel x-direction
  integer(i_def) :: yproc  ! Processor ranks in mesh panel y-direction

  procedure(partitioner_interface), pointer :: partitioner_ptr => null()

  logical(l_def) :: create_2d_mesh                     = .false.
  logical(l_def) :: create_shifted_mesh                = .false.
  logical(l_def) :: create_double_level_mesh           = .false.
  logical(l_def) :: create_multigrid_meshes            = .false.
  logical(l_def) :: create_multigrid_2d_meshes         = .false.
  logical(l_def) :: create_multires_coupling_meshes    = .false.
  logical(l_def) :: create_multires_coupling_2d_meshes = .false.
  logical(l_def) :: create_multigrid                   = .false.

  integer(i_def) :: i, j, n_coupling_meshes, n_chain_meshes

  character(str_def) :: mesh_name
  character(str_def) :: mesh_name_A
  character(str_def) :: mesh_name_B

  ! Currently only quad elements are fully functional
  if (cellshape /= cellshape_quadrilateral) then
    call log_event( "Reference_element must be QUAD for now...", &
                    LOG_LEVEL_ERROR )
  end if

  ! 1.0 Use input args to determine which meshes to create
  !=================================================================
  if ( present(twod_mesh_id) )         create_2d_mesh           = .true.
  if ( present(shifted_mesh_id) )      create_shifted_mesh      = .true.
  if ( present(double_level_mesh_id) ) create_double_level_mesh = .true.

  if ( present(use_multigrid) ) then
    create_multigrid = use_multigrid
    n_chain_meshes = size(chain_mesh_tags)
    if ( create_multigrid ) then
      if ( present(multigrid_mesh_ids) ) then
        create_multigrid_meshes    = .true.
        if (allocated(multigrid_mesh_ids)) deallocate(multigrid_mesh_ids)
        allocate( multigrid_mesh_ids(n_chain_meshes) )
      end if
      if ( present(multigrid_2d_mesh_ids) )then
        create_multigrid_2d_meshes = .true.
        if (allocated(multigrid_2d_mesh_ids)) deallocate(multigrid_2d_mesh_ids)
        allocate( multigrid_2d_mesh_ids(n_chain_meshes) )
      end if
    end if
  end if

  if ( present(use_multires_coupling) ) then
    if ( use_multires_coupling               .and. &
         present(multires_coupling_mesh_ids) .and. &
         present(multires_coupling_mesh_tags)) then
      create_multires_coupling_meshes    = .true.
      if (allocated(multires_coupling_mesh_ids)) deallocate(multires_coupling_mesh_ids)
      allocate(multires_coupling_mesh_ids(size(multires_coupling_mesh_tags)))
    end if
    if ( use_multires_coupling                  .and. &
         present(multires_coupling_2d_mesh_ids) .and. &
         present(multires_coupling_mesh_tags) ) then
      create_multires_coupling_2d_meshes = .true.
    if (allocated(multires_coupling_2d_mesh_ids)) deallocate(multires_coupling_2d_mesh_ids)
      allocate(multires_coupling_2d_mesh_ids(size(multires_coupling_mesh_tags)))
    end if
  end if

  if (offline_partitioning) then

    write (log_scratch_space,'(A)') &
        'Use of offline partitioned mesh files is not yet supported'
    call log_event(trim(log_scratch_space), LOG_LEVEL_ERROR )

  else

    call log_event( "Setting up partition mesh(es)", LOG_LEVEL_INFO )

    ! 2.0 Set constants that will control partitioning.
    !=================================================================
    call set_partition_parameters( total_ranks,       &
                                   xproc, yproc,      &
                                   partitioner_ptr )

    ! 3.0 Read in all global meshes and create local meshes from them.
    !=================================================================
    allocate( global_mesh_collection, source = global_mesh_collection_type() )
    call create_all_base_meshes( local_rank, total_ranks,     &
                                 xproc, yproc,                &
                                 stencil_depth,               &
                                 partitioner_ptr,             &
                                 create_multigrid,            &
                                 multires_coupling_mesh_tags )

    ! 4.0 Read in the global intergrid mesh mappings, then create the
    !     associated local mesh maps
    !=================================================================
    call create_mesh_maps()

  end if


  ! 5.0 Extrude all neshes into 3D local meshes
  !=================================================================
  if (present(multires_coupling_mesh_tags)) then
    call create_all_3d_meshes( create_2d_mesh,                  &
                               create_shifted_mesh,             &
                               create_double_level_mesh,        &
                               create_multigrid_meshes ,        &
                               create_multires_coupling_meshes, &
                               multires_coupling_mesh_tags )
  else
    call create_all_3d_meshes( create_2d_mesh,           &
                               create_shifted_mesh,      &
                               create_double_level_mesh, &
                               create_multigrid_meshes , &
                               create_multires_coupling_meshes )
  end if

  ! 6.0 Assign maps to local 3D meshes
  !=================================================================
  if (create_multigrid_meshes) then
    do i=1, n_chain_meshes-1
      mesh_name_A = chain_mesh_tags(i)
      mesh_name_B = chain_mesh_tags(i+1)
      call add_mesh_maps( mesh_name_A, mesh_name_B )
      multigrid_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
    end do
    multigrid_mesh_ids(n_chain_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
  end if

  if (create_multigrid_2d_meshes) then
    do i=1, n_chain_meshes-1
      mesh_name_A = trim(chain_mesh_tags(i))//'_2d'
      mesh_name_B = trim(chain_mesh_tags(i+1))//'_2d'
      call add_mesh_maps( mesh_name_A, mesh_name_B )
      multigrid_2d_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
    end do
    multigrid_2d_mesh_ids(n_chain_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
  end if

  if (create_multires_coupling_meshes) then
    n_coupling_meshes = size(multires_coupling_mesh_tags)
    if ( n_coupling_meshes == 1 ) then
      multires_coupling_mesh_ids(1) = mesh_collection%get_mesh_id(multires_coupling_mesh_tags(1))
    else
      do i=1, n_coupling_meshes-1
        mesh_name_A = multires_coupling_mesh_tags(i)
        do j=i+1, n_coupling_meshes
          mesh_name_B = multires_coupling_mesh_tags(j)
          call add_mesh_maps( mesh_name_A, mesh_name_B )
        end do
        multires_coupling_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
      end do
      multires_coupling_mesh_ids(n_coupling_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
    end if
  end if

  if (create_multires_coupling_2d_meshes) then
    n_coupling_meshes = size(multires_coupling_mesh_tags)
    if ( n_coupling_meshes == 1 ) then
      mesh_name_A = trim(multires_coupling_mesh_tags(1))//'_2d'
      multires_coupling_2d_mesh_ids(1) = mesh_collection%get_mesh_id(mesh_name_A)
    else
      do i=1, n_coupling_meshes-1
        mesh_name_A = trim(multires_coupling_mesh_tags(i))//'_2d'
        do j=i+1, n_coupling_meshes
          mesh_name_B = trim(multires_coupling_mesh_tags(j))//'_2d'
          call add_mesh_maps( mesh_name_A, mesh_name_B )
        end do
        multires_coupling_2d_mesh_ids(i) = mesh_collection%get_mesh_id(mesh_name_A)
      end do
      multires_coupling_2d_mesh_ids(n_coupling_meshes) = mesh_collection%get_mesh_id(mesh_name_B)
    end if
  end if


  ! 7.0 Extract out mesh ids
  !=================================================================
  mesh_name = prime_mesh_name
  mesh_id = mesh_collection%get_mesh_id(mesh_name)

  if ( create_2d_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_2d'
    twod_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  if ( create_shifted_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_shifted'
    shifted_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  if ( create_double_level_mesh ) then
    mesh_name = trim(prime_mesh_name)//'_double'
    double_level_mesh_id = mesh_collection%get_mesh_id(mesh_name)
  end if

  ! 8.0 Discard global mesh collection and all meshes in it.
  !     (They should not be used past this point in the code)
  !=================================================================
  if (allocated(global_mesh_collection)) deallocate(global_mesh_collection)

end subroutine init_mesh

!===============================================================================
!> @brief Sets common partition parameters to be applied to global meshes.
!         (private subroutine)
!>
!> @param[in]   total_ranks  Total number of MPI ranks in this job
!> @param[out]  xproc              Number of ranks in mesh panel x-direction
!> @param[out]  yproc              Number of ranks in mesh panel y-direction
!> @param[out]  partitioner_ptr    Mesh partitioning strategy
!===============================================================================
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

  integer(i_def), intent(in)  :: total_ranks
  integer(i_def), intent(out) :: xproc
  integer(i_def), intent(out) :: yproc

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  ! Locals
  integer(i_def) :: ranks_per_panel
  integer(i_def) :: start_factor
  integer(i_def) :: end_factor
  integer(i_def) :: fact_count
  logical(l_def) :: found_factors

  character(str_def) :: domain_desc

  integer(i_def), parameter :: max_factor_iters = 10000

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


!===============================================================================
!> @brief   Reads in global meshes from ugrid file, partitions them
!>          and creates local meshes.
!           (private subroutine)
!> @param[in]   local_rank         Number of the local MPI rank
!> @param[in]   total_ranks        Total number of MPI ranks in this job
!> @param[in]   xproc              Number of ranks in mesh panel x-direction
!> @param[in]   yproc              Number of ranks in mesh panel y-direction
!> @param[in]   stencil_depth      Depth of cells outside the base cell
!>                                 of stencil.
!> @param[out]  partitioner_ptr    Mesh partitioning strategy
!> @param[in]   multires_coupling_mesh_tags  Multiresolution Coupling miniapp mesh names
!===============================================================================
subroutine create_all_base_meshes( local_rank, total_ranks,     &
                                   xproc, yproc,                &
                                   stencil_depth,               &
                                   partitioner_ptr,             &
                                   create_multigrid,            &
                                   multires_coupling_mesh_tags )

  use multigrid_config_mod, only: chain_mesh_tags

  implicit none

  integer(i_def), intent(in) :: local_rank
  integer(i_def), intent(in) :: total_ranks
  integer(i_def), intent(in) :: xproc
  integer(i_def), intent(in) :: yproc
  integer(i_def), intent(in) :: stencil_depth
  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr
  logical(l_def), intent(in) :: create_multigrid
  character(str_def), intent(in), optional :: multires_coupling_mesh_tags(:)

  integer(i_def) :: n_panels

  if (geometry == geometry_spherical .and. &
      topology == topology_fully_periodic) then
    n_panels = 6
  else
    n_panels = 1
  end if

  write(log_scratch_space, '(A,I0,A)' )        &
      'Creating global meshes comprising of ', &
      n_panels, ' domain(s)'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! 1.0 Read in prime mesh first by default
  !----------------------------------------------
  write(log_scratch_space,'(A)') &
      'Reading prime global mesh: "'//trim(prime_mesh_name)//'"'
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  call create_base_meshes( [prime_mesh_name], n_panels, &
                           local_rank, total_ranks,     &
                           xproc, yproc,                &
                           stencil_depth,               &
                           partitioner_ptr )

  ! 2.0 Read in any other global meshes required
  !     by other additional configuraed schemes
  !----------------------------------------------
  if (create_multigrid) call create_base_meshes( chain_mesh_tags,         &
                                                 n_panels,                &
                                                 local_rank, total_ranks, &
                                                 xproc, yproc,            &
                                                 stencil_depth,           &
                                                 partitioner_ptr )
  if (present(multires_coupling_mesh_tags))                                   &
                        call create_base_meshes( multires_coupling_mesh_tags, &
                                                 n_panels,                    &
                                                 local_rank, total_ranks,     &
                                                 xproc, yproc,                &
                                                 stencil_depth,               &
                                                 partitioner_ptr )
end subroutine create_all_base_meshes

!===============================================================================
!> @brief Loads the given list of global meshes, partitions them
!>        and creates local meshes from them
!         (private subroutine).
!>
!> @param[in]  mesh_names[:]      Array of requested mesh names to load
!>                                from the mesh input file.
!> @param[in]  n_panels           Number of panel domains in global mesh
!> @param[in]  local_rank         Number of the local MPI rank
!> @param[in]  total_ranks        Total number of MPI ranks in this job
!> @param[in]  xproc              Number of ranks in mesh panel x-direction
!> @param[in]  yproc              Number of ranks in mesh panel y-direction
!> @param[in]  stencil_depth      Depth of cells outside the base cell
!>                                of stencil.
!> @param[in]  partitioner_ptr    Mesh partitioning strategy
!===============================================================================
subroutine create_base_meshes( mesh_names, n_panels,    &
                               local_rank, total_ranks, &
                               xproc, yproc,            &
                               stencil_depth,           &
                               partitioner_ptr )

  implicit none

  character(str_def), intent(in) :: mesh_names(:)
  integer(i_def),     intent(in) :: n_panels
  integer(i_def),     intent(in) :: local_rank
  integer(i_def),     intent(in) :: total_ranks
  integer(i_def),     intent(in) :: xproc
  integer(i_def),     intent(in) :: yproc
  integer(i_def),     intent(in) :: stencil_depth

  procedure(partitioner_interface), intent(in), pointer :: partitioner_ptr

  type(ugrid_mesh_data_type)      :: ugrid_mesh_data
  type(global_mesh_type)          :: global_mesh
  type(global_mesh_type), pointer :: global_mesh_ptr
  type(partition_type)            :: partition
  type(local_mesh_type)           :: local_mesh
  integer(i_def)                  :: local_mesh_id
  integer(i_def)                  :: i

  logical(l_def) :: valid_geometry
  logical(l_def) :: valid_topology


  do i=1, size(mesh_names)
    if (.not. global_mesh_collection%check_for(mesh_names(i))) then

      ! Load mesh data into global_mesh
      call ugrid_mesh_data%read_from_file(trim(filename), mesh_names(i))

      global_mesh = global_mesh_type( ugrid_mesh_data, n_panels )
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

!===============================================================================
!> @brief   Reads in and assigns available global intergrid maps from file.
!           (private subroutine)
!> @details Global meshes which have been read into the models global mesh
!>          collection will have a list of target mesh names. These target mesh
!>          names (if any) indicate the valid intergrid maps avaiable in the
!>          mesh file.
!>
!>          This routine will read in in the appropriate intergrid maps and
!>          assign them to the correct global mesh object.
!===============================================================================
subroutine create_mesh_maps()

  implicit none

  type(ncdf_quad_type) :: file_handler

  character(str_def), allocatable :: source_mesh_names(:)
  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: gid_mesh_map(:,:,:)
  integer(i_def),     allocatable :: lid_mesh_map(:,:,:)

  integer(i_def) :: i, j, n, x, y
  integer(i_def) :: n_meshes

  type(global_mesh_type), pointer :: source_global_mesh => null()

  type(local_mesh_type), pointer :: source_local_mesh => null()
  type(local_mesh_type), pointer :: target_local_mesh => null()

  integer(i_def) :: ntarget_per_source_cell_x, ntarget_per_source_cell_y
  integer(i_def) :: ncells
  integer(i_def) :: target_local_mesh_id

  ! Read in the maps for each global mesh
  !=================================================================
  call file_handler%file_open(trim(filename))

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
          ! Convert global cell ids in the global mesh map
          ! into local cell ids in a local mesh map
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

!===============================================================================
!> @brief Generates the 3D-meshes required by the model configuration.
!         (private subroutine)
!>
!> @details The extrusion types are setup for the required configuration
!>          before 3D-meshes are instantiated.
!>
!> @param[in] create_2d_mesh           Create 2D mesh based on prime mesh
!> @param[in] create_shifted_mesh      Create shifted mesh based on prime mesh
!> @param[in] create_double_level_mesh Create double-level mesh based on prime mesh
!> @param[in] create_multigrid_meshes  Create meshes to support multigrid
!> @param[in] create_multires_coupling_meshes Create meshes for multires_coupling miniapp
!> @param[in] multires_coupling_mesh_tags     Multiresolution Coupling miniapp mesh names
!===============================================================================
subroutine create_all_3D_meshes( create_2d_mesh,                  &
                                 create_shifted_mesh,             &
                                 create_double_level_mesh,        &
                                 create_multigrid_meshes,         &
                                 create_multires_coupling_meshes, &
                                 multires_coupling_mesh_tags )

  use extrusion_mod,           only: TWOD, SHIFTED, DOUBLE_LEVEL
  use extrusion_config_mod,    only: domain_top
  use gungho_extrusion_mod,    only: create_extrusion,         &
                                     create_shifted_extrusion, &
                                     create_double_level_extrusion
  use multigrid_config_mod,    only: chain_mesh_tags

  implicit none

  logical(l_def), intent(in) :: create_2d_mesh
  logical(l_def), intent(in) :: create_shifted_mesh
  logical(l_def), intent(in) :: create_double_level_mesh
  logical(l_def), intent(in) :: create_multigrid_meshes
  logical(l_def), intent(in) :: create_multires_coupling_meshes

  character(str_def), intent(in), optional :: multires_coupling_mesh_tags(:)

  class(extrusion_type), allocatable :: extrusion
  class(extrusion_type), allocatable :: extrusion_shifted
  class(extrusion_type), allocatable :: extrusion_double
  type(uniform_extrusion_type)       :: extrusion_2d

  character(str_def) :: mesh_name

  integer(i_native) :: i

  integer(i_def), parameter :: one_layer = 1_i_def
  real(r_def),    parameter :: atmos_bottom = 0.0_r_def

  ! 1.0 Prime Mesh
  !===================================================================
  ! 1.1 Always generate the prime 3D mesh
  allocate( extrusion, &
            source=create_extrusion() )

  call create_3d_mesh( prime_mesh_name, &
                       extrusion )

  ! 2.0 Generate addition 3d mesh partitions based on the prime mesh
  ! NOTE: This includes 2D meshes as they are currently implemented
  !       as a 3D mesh of 1-layer thick.
  if (create_2d_mesh) then
    extrusion_2d = uniform_extrusion_type( atmos_bottom, &
                                           domain_top,   &
                                           one_layer,    &
                                           TWOD )

    mesh_name = trim(prime_mesh_name)//'_2d'

    call create_3d_mesh( prime_mesh_name, &
                         extrusion_2d,    &
                         mesh_name=mesh_name )
  end if

  ! 3.0 Generate addition shifted 3d mesh partition
  !     based on the prime mesh.
  if (create_shifted_mesh) then

    allocate( extrusion_shifted, &
              source=create_shifted_extrusion(extrusion) )

    mesh_name = trim(prime_mesh_name)//'_shifted'

    call create_3d_mesh( prime_mesh_name,   &
                         extrusion_shifted, &
                         mesh_name=mesh_name )
  end if

  ! 4.0 Generate addition double-level 3d mesh partition
  !     based on the prime mesh.
  if (create_double_level_mesh) then

    allocate( extrusion_double, &
              source=create_double_level_extrusion(extrusion) )

    mesh_name = trim(prime_mesh_name)//'_double'

    call create_3d_mesh( prime_mesh_name,  &
                         extrusion_double, &
                         mesh_name=mesh_name )
  end if

  ! 5.0 Generate meshes required by any other schemes
  !     in the model run configuration.
  !===================================================================
  ! 5.1 Dynamics Multigrid
  if (create_multigrid_meshes) then

    do i=1, size(chain_mesh_tags)
      call create_3d_mesh( chain_mesh_tags(i), &
                           extrusion )

      mesh_name = trim(chain_mesh_tags(i))//'_2d'
      call create_3d_mesh( chain_mesh_tags(i), &
                           extrusion_2d,       &
                           mesh_name=mesh_name )
    end do

  end if

  ! 5.2 Multires Coupling Miniapp
  if (create_multires_coupling_meshes) then

    do i=1, size(multires_coupling_mesh_tags)
      call create_3d_mesh( multires_coupling_mesh_tags(i), &
                           extrusion )

      mesh_name = trim(multires_coupling_mesh_tags(i))//'_2d'
      call create_3d_mesh( multires_coupling_mesh_tags(i), &
                           extrusion_2d,                   &
                           mesh_name=mesh_name )
    end do

  end if

  if (allocated(extrusion))         deallocate( extrusion )
  if (allocated(extrusion_shifted)) deallocate( extrusion_shifted )
  if (allocated(extrusion_double))  deallocate( extrusion_double )

  return
end subroutine create_all_3D_meshes

!===============================================================================
!> @brief   Generates a single (partitioned) 3D-mesh.
!           (private subroutine)
!> @details Instantiates a 3d-mesh partition and adds it to the model's
!>          mesh collection. Multiple meshes may be generated in the model
!>          based on the same global mesh but with differing extrusions.
!>
!> @param[in] base_mesh_name  Name of base 2D-mesh
!> @param[in] extrusion       Extrusion type for this 3D-mesh
!> @param[in] mesh_name       Optional, Name of local 3D-mesh,
!>                            defaults to base_mesh_name
!===============================================================================
subroutine create_3d_mesh( base_mesh_name, &
                           extrusion,      &
                           mesh_name )

  use mesh_collection_mod, only: mesh_collection

  implicit none

  character(str_def),    intent(in) :: base_mesh_name
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: mesh_name

  type(local_mesh_type), pointer :: local_mesh_ptr => null()

  type(mesh_type)    :: mesh
  integer(i_def)     :: mesh_id
  character(str_def) :: name

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

!===============================================================================
!> @brief   Assigns intergrid maps to 3D-mesh partitions
!           (private subroutine).
!> @details Adds local ID integrid mappings 3D-mesh parititons.
!>          Local integrid maps are assign to to source (source-> target)
!>          and target (target -> source) meshes. A number of
!>          of assumptions are made.
!>
!>          *  The base global meshes of source/target 3D-mesh
!>             are present in the global_mesh_collection.
!>          *  Integrid maps for the base global meshes where
!>             where read in and assigned.
!>          *  Partitioned 3D-mesh are present in the
!>             mesh_collection.
!>
!> @param[in] source_mesh_name   Name of source 3D mesh partition
!> @param[in] target_mesh_name   Name of target 3D mesh partition
!===============================================================================
subroutine add_mesh_maps( source_mesh_name, &
                          target_mesh_name )

  use mesh_collection_mod, only: mesh_collection

  implicit none

  character(str_def), intent(in) :: source_mesh_name
  character(str_def), intent(in) :: target_mesh_name

  type(mesh_type), pointer :: source_mesh => null()
  type(mesh_type), pointer :: target_mesh => null()


  ! Now add in any mesh maps required by multigrid or
  ! the multires_coupling miniapp
  source_mesh => mesh_collection % get_mesh( source_mesh_name )
  target_mesh => mesh_collection % get_mesh( target_mesh_name )

  if ( associated(source_mesh) .and. &
       associated(target_mesh) ) then

    ! Mesh tag names may be different but "could point to the same mesh
    ! So check the ids are not the same
    if (source_mesh%get_id() == target_mesh%get_id()) then
      write(log_scratch_space,'(A)')                  &
          'Unable to create intergrid map: Source('// &
          trim(source_mesh_name)//' and target('//    &
          trim(target_mesh_name)//') mesh ids are the same'
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

!===============================================================================
!> @brief Finalises the mesh_collection
!===============================================================================
subroutine final_mesh()

  use mesh_collection_mod, only: mesh_collection

  implicit none

  if (allocated(mesh_collection)) then
    call mesh_collection%clear()
    deallocate(mesh_collection)
  end if

  return
end subroutine final_mesh

end module create_mesh_mod
