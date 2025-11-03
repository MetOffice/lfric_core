!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!>
!> @brief   Set up specified mesh(es) from global/local mesh input file(s).
!> @details This routine will create a mesh_object_type(s) from a
!>          specified mesh input file and extrusion.
!>
!>          The algorithm differs depending on whether the input files(s)
!>          are prepartitioned (local mesh files) or not (global mesh files).
!>
!>          The result will be:
!>            * A set of local mesh objects stored in the application local
!>              mesh collection object.
!>            * A set of mesh objects stored in the application mesh collection
!>              object.
!>
!>          Local mesh object names will use the same mesh name as given in the
!>          input file. Extruded meshes are allowed to have an alternative mesh
!>          name as several meshes could exist in memory based on the same local
!>          mesh object.
!>
module driver_mesh_mod

  use add_mesh_map_mod,           only: assign_mesh_maps
  use constants_mod,              only: i_def, l_def, str_def, &
                                        str_max_filename
  use check_global_mesh_mod,      only: check_global_mesh
  use check_local_mesh_mod,       only: check_local_mesh
  use create_mesh_mod,            only: create_extrusion, create_mesh
  use extrusion_mod,              only: extrusion_type
  use global_mesh_mod,            only: global_mesh_type
  use load_global_mesh_mod,       only: load_global_mesh
  use load_local_mesh_mod,        only: load_local_mesh
  use load_local_mesh_maps_mod,   only: load_local_mesh_maps
  use log_mod,                    only: log_event,         &
                                        log_scratch_space, &
                                        log_level_debug,   &
                                        log_level_error
  use namelist_collection_mod,    only: namelist_collection_type
  use namelist_mod,               only: namelist_type
  use panel_decomposition_mod,    only: panel_decomposition_type
  use partition_mod,              only: partitioner_interface

  use runtime_partition_lfric_mod, only: get_partition_parameters
  use runtime_partition_mod,       only: mesh_cubedsphere,       &
                                         mesh_planar,            &
                                         create_local_mesh_maps, &
                                         create_local_mesh

  use global_mesh_collection_mod, only: global_mesh_collection
  use local_mesh_collection_mod,  only: local_mesh_collection

  ! Configuration modules
  use finite_element_config_mod, only: cellshape_quadrilateral
  use base_mesh_config_mod,      only: geometry_spherical, &
                                       topology_fully_periodic

  implicit none

  private
  public :: init_mesh

contains
!=======================================


!===============================================================================
!> @brief  Generates mesh(es) from mesh input file(s) on a given extrusion.
!>
!> @param[in] configuration     Application configuration object.
!>                              This configuration object should contain the
!>                              following defined namelist objects:
!>                                 * base_mesh
!>                                 * finite_element
!>                                 * partititioning
!> @param[in] local_rank        The MPI rank of this process.
!> @param[in] total_ranks       Total number of MPI ranks in this job.
!> @param[in] mesh_names        Mesh names to load from the mesh input file(s).
!> @param[in] extrusion         Extrusion object to be applied to meshes.
!> @param[in] stencil_depth     Required stencil depth for the application.
!> @param[in] check_partitions  Apply check for even partitions with the
!>                              configured partition stratedy.
!>                              (unpartitioned mesh input only)
!> @param[in] alt_names         (Optional), Alternative names for meshes in the
!>                                          application mesh collection object.
!===============================================================================
subroutine init_mesh( configuration,           &
                      local_rank, total_ranks, &
                      mesh_names, extrusion,   &
                      stencil_depth,           &
                      check_partitions,        &
                      alt_names )

  implicit none

  ! Arguments
  type(namelist_collection_type) :: configuration

  integer(i_def),        intent(in) :: local_rank
  integer(i_def),        intent(in) :: total_ranks
  character(str_def),    intent(in) :: mesh_names(:)
  class(extrusion_type), intent(in) :: extrusion

  integer(i_def),    intent(in) :: stencil_depth
  logical(l_def),    intent(in) :: check_partitions

  character(str_def), optional, intent(in) :: alt_names(:)

  ! Parameters
  character(len=9), parameter :: routine_name = 'init_mesh'

  ! Namelist variables
  type(namelist_type), pointer :: base_mesh_nml
  type(namelist_type), pointer :: finite_element_nml
  type(namelist_type), pointer :: partitioning_nml

  character(str_max_filename)  :: file_prefix

  integer(i_def) :: cellshape

  logical :: prepartitioned
  logical :: generate_inner_halos

  integer :: geometry
  integer :: topology
  integer :: mesh_selection

  ! Local variables
  character(str_def), allocatable :: names(:)
  character(str_def), allocatable :: tmp_mesh_names(:)
  character(str_max_filename)     :: input_mesh_file

  procedure(partitioner_interface), pointer :: partitioner_ptr

  class(panel_decomposition_type), allocatable :: decomposition

  integer(i_def)     :: n_digit
  character(str_def) :: fmt_str, number_str

  !============================================================================
  ! 0.0 Extract configuration variables
  !============================================================================
  base_mesh_nml      => configuration%get_namelist('base_mesh')
  finite_element_nml => configuration%get_namelist('finite_element')

  call base_mesh_nml%get_value( 'prepartitioned', prepartitioned )
  call base_mesh_nml%get_value( 'file_prefix',    file_prefix )
  call finite_element_nml%get_value( 'cellshape', cellshape )



  !============================================================================
  ! 0.1 Some basic checks
  !============================================================================
  ! Set up stencil depth
  if (stencil_depth < 0_i_def) then
    write(log_scratch_space,'(A)') &
       'Standard partitioned meshes must support a not -ve stencil_depth'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Currently only quad elements are fully functional
  if (cellshape /= CELLSHAPE_QUADRILATERAL) then
    call log_event( "Reference_element must be QUAD for now...", &
                    LOG_LEVEL_ERROR )
  end if


  !============================================================================
  ! 1.0 Determine which names to apply to resultant meshes.
  !============================================================================
  if (present(alt_names)) then
    if (size(alt_names) == size(mesh_names)) then
      allocate(names, source=alt_names)
    else
      write(log_scratch_space, '(A)')                   &
          'Specified alternative mesh names to does '// &
          'not match number of requested meshes.'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
  else
    allocate(names, source=mesh_names)
  end if


  !===========================================================================
  ! 2.0 Create local mesh objects:
  !     Two code pathes presented, either:
  !     1. The input files have been pre-partitioned.
  !        Meshes and are simply read from file and local mesh objects
  !        are populated.
  !     2. The input files have not been partitioned.
  !        Global meshes are loaded from file and partitioning is applied
  !        at runtime.  NOTE: This option is provided as legacy, and support
  !        is on a best endeavours basis.
  !===========================================================================

  generate_inner_halos = .false.

  if (prepartitioned) then

    !==========================================================================
    ! 2.1 Read in local meshes / partition information / mesh maps
    !     direct from file.
    !==========================================================================
    !
    ! For this local rank, a mesh input file with a common base name
    ! of the following form should exist.
    !
    !   <input_basename>_<local_rank>_<total_ranks>.nc
    !
    ! Where 1 rank is assigned to each mesh partition.
    n_digit = int(log10(real(total_ranks))) + 1
    write(fmt_str, '(A,I0,A,I0,A)') "(I", n_digit, ".", n_digit, ")"
    write(number_str, fmt_str) local_rank
    write(input_mesh_file, '(A, "_", A, "-", I0, ".nc")') &
        trim(file_prefix), trim(number_str), total_ranks

    call log_event( 'Using pre-partitioned mesh file:', log_level_debug )
    call log_event( '   '//trim(input_mesh_file), log_level_debug )
    call log_event( "Loading local mesh(es)", log_level_debug )


    ! 2.1a Read in all local mesh data for this rank and
    !      initialise local mesh objects from them.
    !===========================================================
    ! Each partitioned mesh file will contain meshes of the
    ! same name as all other partitions.
    call load_local_mesh( input_mesh_file, mesh_names )

    ! 2.1b Apply configuration related checks to ensure that these
    !      meshes are suitable for the supplied application
    !      configuration.
    !===========================================================
    call check_local_mesh( configuration, &
                           stencil_depth, &
                           mesh_names )

    ! 2.1c Load and assign mesh maps.
    !===========================================================
    ! Mesh map identifiers are determined by the source/target
    ! mesh IDs they relate to. As a result inter-grid mesh maps
    ! need to be loaded after the relevant local meshes have
    ! been loaded.
    tmp_mesh_names = local_mesh_collection%get_mesh_names()
    call load_local_mesh_maps( input_mesh_file, tmp_mesh_names )
    if (allocated(tmp_mesh_names)) deallocate(tmp_mesh_names)

  else

    !==========================================================================
    ! 2.2 Perform runtime partitioning of global meshes.
    !==========================================================================
    call base_mesh_nml%get_value( 'geometry', geometry )
    call base_mesh_nml%get_value( 'topology', topology )

    partitioning_nml => configuration%get_namelist('partitioning')

    call partitioning_nml%get_value( 'generate_inner_halos', &
                                      generate_inner_halos )

    if ( geometry == geometry_spherical .and. &
         topology == topology_fully_periodic ) then
      mesh_selection = mesh_cubedsphere
      call log_event( "Setting up cubed-sphere partition mesh(es)", &
                      log_level_debug )
    else
      mesh_selection = mesh_planar
      call log_event( "Setting up planar partition mesh(es)", &
                      log_level_debug )
    end if
    write(input_mesh_file,'(A)') trim(file_prefix) // '.nc'

    ! 2.2a Set constants that will control partitioning.
    !===========================================================
    call get_partition_parameters( configuration, mesh_selection, &
                                   total_ranks, decomposition,    &
                                   partitioner_ptr )

    ! 2.2b Read in all global meshes from input file
    !===========================================================
    call load_global_mesh( input_mesh_file, mesh_names )

    ! 2.2c Apply configuration related checks to ensure that these
    !      meshes are suitable for the supplied application
    !      configuration.
    !===========================================================
    call check_global_mesh( configuration, mesh_names )

    ! 2.2e Partition the global meshes
    !===========================================================
    call create_local_mesh( mesh_names,              &
                            local_rank, total_ranks, &
                            decomposition,           &
                            stencil_depth,           &
                            generate_inner_halos,    &
                            partitioner_ptr )


    ! 2.2f Read in the global intergrid mesh mappings,
    !      then create the associated local mesh maps
    !===========================================================
    call create_local_mesh_maps( input_mesh_file )

    ! Clear the global mesh
    call global_mesh_collection%clear()

  end if  ! prepartitioned


  !============================================================================
  ! 3.0 Extrude the specified meshes from local mesh objects into
  !     mesh objects on the given extrusion.
  !============================================================================
  call create_mesh( mesh_names, extrusion, alt_name=names )


  !============================================================================
  ! 4.0 Generate intergrid LiD-LiD maps and assign them to mesh objects.
  !============================================================================
  call assign_mesh_maps(mesh_names)

end subroutine init_mesh

end module driver_mesh_mod
