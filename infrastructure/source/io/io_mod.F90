!-------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-------------------------------------------------------------

!>  @brief Module for IO subroutines
!!
!!  @details Holds all setup routines for the I/O system as
!!           well as the I/O handler routines for reading and
!!           writing LFRic fields
!-------------------------------------------------------------------------------
module io_mod
  use constants_mod,                 only: i_def, i_native, i_halo_index, &
                                           r_def, dp_xios,                &
                                           str_def, str_max_filename,     &
                                           l_def, PI, radians_to_degrees
  use field_mod,                     only: field_type, field_proxy_type
  use field_collection_mod,          only: field_collection_type, &
                                           field_collection_iterator_type
  use finite_element_config_mod,     only: element_order
  use base_mesh_config_mod,          only: geometry, &
                                           geometry_spherical
  use fs_continuity_mod,             only: W0, W1, W2, W3, Wtheta, W2H, &
                                           name_from_functionspace
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection 
  use function_space_mod,            only: function_space_type, BASIS
  use function_space_collection_mod, only: function_space_collection
  use project_output_mod,            only: project_output
  use io_config_mod,                 only: diagnostic_frequency, &
                                           checkpoint_write,     &
                                           checkpoint_read,      &
                                           write_dump
  use initialization_config_mod,     only: init_option,               &
                                           init_option_fd_start_dump, &
                                           ancil_option,              &
                                           ancil_option_aquaplanet

  use files_config_mod,              only: checkpoint_stem_name, &
                                           start_dump_filename, &
                                           start_dump_directory
  use time_config_mod,               only: timestep_start, &
                                           timestep_end
  use runtime_constants_mod,         only: get_coordinates
  use coord_transform_mod,           only: xyz2llr
  use log_mod,                       only: log_event,         &
                                           log_set_level,     &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_DEBUG,   &
                                           LOG_LEVEL_TRACE

  use psykal_lite_mod,               only: invoke_nodal_coordinates_kernel, &
                                           invoke_pointwise_convert_xyz2llr
  use mpi_mod, only: get_comm_size, get_comm_rank, all_gather
  use mpi, only: MPI_SUCCESS
  use xios

  implicit none
  private
  public :: ts_fname,                     &
            checkpoint_write_netcdf,      &
            checkpoint_write_xios,        &
            checkpoint_read_netcdf,       &
            checkpoint_read_xios,         &
            write_checkpoint,             &
            read_checkpoint,              &
            dump_write_xios,              &
            dump_read_xios,               &
            write_state,                  &
            read_state,                   &
            xios_domain_init,             &
            nodal_write_field,            &
            xios_write_field_node,        &
            xios_write_field_face,        &
            xios_write_field_edge,        &
            xios_write_field_single_face, &
            xios_read_field_face,         &
            xios_read_field_single_face

contains

!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS context and domain initialisation
!!
!!  @details  Initialises the XIOS context and calls further routines to
!!            setup various domains for read and write
!!
!!  @param[in]      xios_ctx      XIOS context identifier
!!  @param[in]      mpi_comm      The MPI comm object
!!  @param[in]      dtime         XIOS timestep interval
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      twod_mesh_id  2D Mesh id
!!  @param[in]      chi           Coordinate field
!-------------------------------------------------------------------------------

subroutine xios_domain_init(xios_ctx, mpi_comm, dtime, &
                            mesh_id, twod_mesh_id,  chi)

  use fs_continuity_mod, only : name_from_functionspace

  implicit none

  ! Arguments
  character(len=*),   intent(in)       :: xios_ctx
  integer(i_def),     intent(in)       :: mpi_comm
  integer(i_def),     intent(in)       :: dtime
  integer(i_def),     intent(in)       :: mesh_id
  integer(i_def),     intent(in)       :: twod_mesh_id
  type(field_type),   intent(in)       :: chi(:)


  ! Local variables 
  type(xios_duration)                  :: xios_timestep
  type(xios_duration)                  :: o_freq, cp_freq, dump_freq
  type(xios_duration)                  :: av_freq
  type(xios_context)                   :: xios_ctx_hdl
  type(xios_file)                      :: cpfile_hdl, rsfile_hdl, &
                                          ofile_hdl, dumpfile_hdl
  type(xios_fieldgroup)                :: cpfieldgroup_hdl, &
                                          fdfieldgroup_hdl
  character(len=str_max_filename)      :: checkpoint_write_fname
  character(len=str_max_filename)      :: checkpoint_read_fname
  character(len=str_max_filename)      :: dump_fname
  character(len=str_def)               :: domain_name, domain_fs_name
  integer(i_native), parameter         :: domain_function_spaces(5) &
                                                  = (/W0, W1, W2, W3, Wtheta/)

  integer(i_native) :: fs_index


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Setup context !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call xios_context_initialize(xios_ctx, mpi_comm)
  call xios_get_handle(xios_ctx, xios_ctx_hdl)
  call xios_set_current_context(xios_ctx_hdl)

  !!!!!!!!!!!!!!!!!!!!!!!!!!! Setup diagnostic domains !!!!!!!!!!!!!!!!!!!!!!!!!!

  call xios_diagnostic_domain_init(mesh_id, chi)

  !!!!!!!!!!!!!!!!!!!!!! Setup checkpoint domains !!!!!!!!!!!!!!!!!!!!!

  ! Create all the regular checkpoint domains based on current function spaces
  ! Loop over function spaces we need to create domains for

  do fs_index = lbound(domain_function_spaces, 1), &
                ubound(domain_function_spaces, 1)

    domain_fs_name = name_from_functionspace(domain_function_spaces(fs_index))
    domain_name = "checkpoint_" // trim(domain_fs_name)

    ! Enable use of the XIOS i_index for W3 and Wtheta

    if (domain_function_spaces(fs_index) == W3 .or. &
        domain_function_spaces(fs_index) == Wtheta) then

      call xios_checkpoint_domain_init(domain_function_spaces(fs_index), &
                                           trim(domain_name), mesh_id, chi, .true.)
    else 
    
      call xios_checkpoint_domain_init(domain_function_spaces(fs_index), &
                                       trim(domain_name), mesh_id, chi, .false.)
    end if

  end do


  ! Set up 2D checkpoint domain - only W3 at the moment

  domain_name = "checkpoint_W3_2D"

  call xios_checkpoint_domain_init(W3, domain_name,  twod_mesh_id, chi, .true.)

  !!!!!!!!!!!!! Setup diagnostic output context information !!!!!!!!!!!!!!!!!!

  ! Set diagnostic output (configured in timesteps) frequency in seconds
  o_freq%second =  diagnostic_frequency*dtime

  call xios_get_handle("lfric_diag",ofile_hdl)
  call xios_set_attr(ofile_hdl, output_freq=o_freq)

  ! Set diagnostic output (configured in timesteps) frequency in seconds
  if (xios_is_valid_file("lfric_averages")) then
    av_freq%second = timestep_end*dtime
    
    call xios_get_handle("lfric_averages",ofile_hdl)
    call xios_set_attr(ofile_hdl, output_freq=av_freq)
  end if

 !!!!!!!!!!!!! Setup dump context information !!!!!!!!!!!!!!!!!!

  if ( write_dump ) then

    ! Enable the fd field group

    call xios_get_handle("physics_fd_fields",fdfieldgroup_hdl)
    call xios_set_attr(fdfieldgroup_hdl, enabled=.true.)

    ! Create dump filename from base name and end timestep
    write(dump_fname,'(A,A,I6.6)') &
       trim(start_dump_directory)//'/'//trim(start_dump_filename),"_", timestep_end


    ! Set dump frequency (end timestep) in seconds
    dump_freq%second = timestep_end*dtime

    call xios_get_handle("lfric_fd_dump",dumpfile_hdl)
    call xios_set_attr(dumpfile_hdl,name=dump_fname, enabled=.true.)
    call xios_set_attr(dumpfile_hdl, output_freq=dump_freq)


  end if

  if( init_option == init_option_fd_start_dump .or. &
      ancil_option == ancil_option_aquaplanet ) then    

    ! Enable the fd field group

    call xios_get_handle("physics_fd_fields",fdfieldgroup_hdl)
    call xios_set_attr(fdfieldgroup_hdl, enabled=.true.)

    ! Create dump filename from stem
    write(dump_fname,'(A)') trim(start_dump_directory)//'/'//trim(start_dump_filename)
    ! Set dump filename
    call xios_get_handle("read_lfric_fd_dump",dumpfile_hdl)
    call xios_set_attr(dumpfile_hdl,name=dump_fname, enabled=.true.)


  end if

  !!!!!!!!!!!!! Setup checkpoint context information !!!!!!!!!!!!!!!!!!

  ! Enable the checkpoint field group

  if ( checkpoint_write .or. checkpoint_read ) then

    call xios_get_handle("checkpoint_fields",cpfieldgroup_hdl)
    call xios_set_attr(cpfieldgroup_hdl, enabled=.true.)

  end if

  ! Checkpoint writing
  if ( checkpoint_write ) then

    ! Get the checkpoint file definition handle, set its filename and frequency
    ! and enable/disable as required

    ! Create checkpoint filename from stem and end timestep
    write(checkpoint_write_fname,'(A,A,I6.6)') &
                              trim(checkpoint_stem_name),"_", timestep_end

    ! Set checkpoint frequency (end timestep) in seconds
    cp_freq%second = timestep_end*dtime

    call xios_get_handle("lfric_checkpoint_write",cpfile_hdl)
    call xios_set_attr(cpfile_hdl,name=checkpoint_write_fname, enabled=.true.)
    call xios_set_attr(cpfile_hdl, output_freq=cp_freq)

  end if

  ! Checkpoint reading
  if ( checkpoint_read ) then

    ! Get the checkpoint file definition handle, set its filename
    ! and enable/disable as required

    ! Create checkpoint filename from stem and (start - 1) timestep
    write(checkpoint_read_fname,'(A,A,I6.6)') &
                            trim(checkpoint_stem_name),"_", (timestep_start - 1)

    call xios_get_handle("lfric_checkpoint_read",rsfile_hdl)
    call xios_set_attr(rsfile_hdl, name=checkpoint_read_fname, enabled=.true.)


  end if


  !!!!!!!!!!!!!!!!!!!!!! Setup calendar and finalise context !!!!!!!!!!!!!!!!!!!!

  xios_timestep%second = dtime

  call xios_set_timestep(xios_timestep)

  call xios_close_context_definition()

  return

end subroutine xios_domain_init

!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS diagnostic domain initialisation
!!
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      chi           Coordinate field
!-------------------------------------------------------------------------------

subroutine xios_diagnostic_domain_init(mesh_id, chi)

  implicit none

  ! Arguments

  integer(i_def),   intent(in) :: mesh_id
  type(field_type), intent(in) :: chi(:)

  ! Local variables 

  integer(i_def) :: i


  ! Node domain (W0)
  integer(i_def)             :: ibegin_nodes
  integer(i_def)             :: coord_dim_full
  integer(i_def)             :: coord_dim_owned
  real(r_def),allocatable    :: nodes_lon_full(:)
  real(r_def),allocatable    :: nodes_lat_full(:)
  real(dp_xios),allocatable  :: nodes_lon(:)
  real(dp_xios),allocatable  :: nodes_lat(:)
  real(dp_xios),allocatable  :: bnd_nodes_lon(:,:)
  real(dp_xios),allocatable  :: bnd_nodes_lat(:,:)
  real(dp_xios),pointer      :: fractional_levels_nodes(:) => null()

  ! Face domain (W3)
  integer(i_def)             :: ibegin_faces
  real(dp_xios),allocatable  :: faces_lon(:)
  real(dp_xios),allocatable  :: faces_lat(:)
  real(dp_xios),allocatable  :: bnd_faces_lon(:,:)
  real(dp_xios),allocatable  :: bnd_faces_lat(:,:)
  real(dp_xios),pointer      :: fractional_levels_full_faces(:) => null()
  real(dp_xios),pointer      :: fractional_levels_half_faces(:) => null()

  ! Edge domain on half levels (W2H)
  integer(i_def)             :: ibegin_edges
  real(dp_xios),allocatable  :: edges_lon(:)
  real(dp_xios),allocatable  :: edges_lat(:)
  real(dp_xios),allocatable  :: bnd_edges_lon(:,:)
  real(dp_xios),allocatable  :: bnd_edges_lat(:,:)
  real(dp_xios),pointer      :: fractional_levels_half_edges(:) => null()

  ! Levels variables
  integer(i_def)             :: nfull_levels
  integer(i_def)             :: nhalf_levels

  ! domain index to ensure processor decomposition independent
  ! read for the face domain
  integer(i_def),allocatable :: domain_index(:)

  ! Variables needed to compute output domain coordinates in lat-long

  ! Transformed coords for nodal output
  type( field_type ) :: coord_output(3)
  ! Field proxies (to calculate domain coordinate info)
  type(field_proxy_type), target  :: proxy_coord_output(3)

  ! Variables for mesh information
  type(mesh_type), pointer :: local_mesh => null()
  integer(i_def)           :: num_face_local
  integer(i_def)           :: nodes_per_face
  integer(i_def)           :: nodes_per_edge

  type(function_space_type), pointer :: output_field_fs   => null()
  type(function_space_type), pointer :: w2h_fs   => null()

  ! Variables for the gather to determine global domain sizes
  ! from the local partitioned ones

  integer(i_def)                :: global_undf, size_w2h, levs_w2h
  integer(i_def), allocatable   :: local_undf(:), all_undfs(:)
  integer(i_def)                :: local_annexed_dof

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def)                :: r2d


  if ( geometry == geometry_spherical ) then
   r2d = radians_to_degrees
  else
   r2d = 1.0_r_def
  endif

  ! Set up arrays to hold number of dofs for local and global domains

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))

  
  all_undfs = 0

  ! Set up the 'node' domain.
  ! Here we use information from W0 to calculate the physical coordinates
  ! for the horizontal domain and levels for the vertical domain
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W0 )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do

  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get pointer to local mesh for the coordinate field
  local_mesh => coord_output(1)%get_mesh()

  ! Get mesh information 
  num_face_local = local_mesh%get_last_edge_cell()
  nodes_per_face = local_mesh%get_nverts_per_cell_2d()
  nodes_per_edge = local_mesh%get_nverts_per_edge()

  ! Calculate the local size of a W2H fs in order to determine
  ! how many edge dofs for the current partition
  w2h_fs => function_space_collection%get_fs( mesh_id, element_order, W2H )
  levs_w2h = size(w2h_fs%get_levels())
  size_w2h = w2h_fs%get_last_dof_owned()/levs_w2h

  ! Get the local value for last owned dof

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the local value for last annexed dof

  local_annexed_dof  = proxy_coord_output(1)%vspace%get_last_dof_annexed()

  ! Get the unique fractional levels to set up vertical output domain
  fractional_levels_nodes => proxy_coord_output(1)%vspace%get_levels()

  nfull_levels = size(fractional_levels_nodes)

  ! Allocate coordinate arrays

  ! coord_dim_full is the size of one whole level of the full field
  ! needed to be sure we get all coords for faces on this partition

  coord_dim_full = size(proxy_coord_output(1)%data) / nfull_levels

  ! coord_dim_owned is the size up to last owned dof for a whole level
  ! this is needed to set the node domain for this partition

  coord_dim_owned = local_undf(1) / nfull_levels

  allocate(nodes_lon_full(coord_dim_full))
  allocate(nodes_lat_full(coord_dim_full))

  nodes_lon_full = 0.0_r_def
  nodes_lat_full = 0.0_r_def


  allocate(nodes_lon( coord_dim_owned ))
  allocate(nodes_lat( coord_dim_owned ))

  nodes_lon = 0.0_dp_xios
  nodes_lat = 0.0_dp_xios

  allocate(bnd_nodes_lon(1,size(nodes_lon)))
  allocate(bnd_nodes_lat(1,size(nodes_lat)))

  allocate(bnd_faces_lon(nodes_per_face,num_face_local))
  allocate(bnd_faces_lat(nodes_per_face,num_face_local))  

  allocate(bnd_edges_lon(nodes_per_edge,size_w2h))
  allocate(bnd_edges_lat(nodes_per_edge,size_w2h))  

  ! Calculate the node coords arrays and also the face-node boundary arrays 
  call calc_xios_domain_coords(local_mesh, coord_output, chi,  &
                               nfull_levels, num_face_local,   &
                               nodes_lon_full, nodes_lat_full, &
                               bnd_faces_lon, bnd_faces_lat,   &
                               bnd_edges_lon, bnd_edges_lat)

  ! Get nodal coordinates (owned part of full length arrays)
  nodes_lon =  nodes_lon_full(1:coord_dim_owned)
  nodes_lat =  nodes_lat_full(1:coord_dim_owned)

  ! Construct nodal bounds arrays
  bnd_nodes_lon=(reshape(nodes_lon, (/1, size(nodes_lon)/) ) )
  bnd_nodes_lat=(reshape(nodes_lat, (/1, size(nodes_lat)/) ) )

  !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call all_gather ( local_undf, all_undfs, 1 )

  ! Adjust size of data taking into account how many levels we have (same for each
  ! partition as we only partition horizontally)

  all_undfs = all_undfs/nfull_levels

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for xios node domain

  global_undf = sum(all_undfs)

  ! Calculate ibegin for each rank (as we have the array of undfs in order
  ! we can just sum to get it)

  if (get_comm_rank() == 0) then
    ibegin_nodes = 0
  else
    ibegin_nodes = sum(all_undfs(1:get_comm_rank()))
  end if


 ! Do Node domain setup

  call xios_set_domain_attr("node", ni_glo=global_undf,           &
                            ibegin=ibegin_nodes,                  &
                            ni=local_undf(1)/nfull_levels,        &
                            type='unstructured')
  call xios_set_domain_attr("node", lonvalue_1d=nodes_lon,        &
                            latvalue_1d=nodes_lat)
  call xios_set_domain_attr("node", bounds_lon_1d=bnd_nodes_lon,  &
                            bounds_lat_1d=bnd_nodes_lat)

  call xios_set_axis_attr("vert_axis_full_levels", &
                           n_glo=nfull_levels,     &
                           value=fractional_levels_nodes)


  ! Clean up things not needed or for reuse in face domain setup
  deallocate(local_undf, all_undfs)
  fractional_levels_nodes => null()
  output_field_fs => null()

  ! Set up arrays for AllGather

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))
  
  all_undfs = 0

  ! Set up the face domains
  ! Here we use information from W3 to calculate the physical coordinates
  ! for the horizontal domains and the 'half levels' vertical domain
  ! We use Wtheta to set the 'full levels' vertical domain
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W3 )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do


  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == geometry_spherical ) then
     call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if


  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the unique fractional levels to set up half levels vertical output domain
  fractional_levels_half_faces => proxy_coord_output(1)%vspace%get_levels()

  nhalf_levels = size(fractional_levels_half_faces)

  ! Allocate coordinate arrays for faces

  allocate(faces_lon( num_face_local))
  allocate(faces_lat( num_face_local))

  faces_lon =  proxy_coord_output(1)%data(1: local_undf(1):nhalf_levels) * r2d
  faces_lat =  proxy_coord_output(2)%data(1: local_undf(1):nhalf_levels) * r2d

  !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call all_gather ( local_undf, all_undfs, 1 )

  ! Adjust size of data taking into account how many levels we have (same for each
  ! partition as we only partition horizontally)

  all_undfs = all_undfs/nhalf_levels

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for xios face domain

  global_undf = sum(all_undfs)


  ! Calculate ibegin for each rank as we have the array of undfs in order
  ! we can just sum to get it.

  if (get_comm_rank() == 0) then
    ibegin_faces = 0
  else
    ibegin_faces = sum(all_undfs(1:get_comm_rank()))
  end if


  call xios_set_domain_attr("face_half_levels", ni_glo=global_undf,          &
                            ibegin=ibegin_faces,                             &
                            ni=local_undf(1)/nhalf_levels,                   &
                            type='unstructured')
  call xios_set_domain_attr("face_half_levels", lonvalue_1d=faces_lon,       &
                            latvalue_1d=faces_lat)
  call xios_set_domain_attr("face_half_levels", bounds_lon_1d=bnd_faces_lon, &
                            bounds_lat_1d=bnd_faces_lat)

  ! Set up the domain index (this will be used for face half level and 
  ! face full level domain setup 

  ! Allocate domain_index for faces
  allocate(domain_index(num_face_local))

  ! Populate domain_index for this rank 
  call proxy_coord_output(1)%vspace%get_global_dof_id_2d(domain_index)

  ! Pass local portion of domain_index
  call xios_set_domain_attr("face_half_levels", &
                            i_index=int(domain_index(1:(local_undf(1)/nhalf_levels))))

  call xios_set_axis_attr("vert_axis_half_levels", n_glo=nhalf_levels,       &
                          value=fractional_levels_half_faces)

  ! Clean up things ready to reuse for the vertical full levels domain setup
  deallocate(local_undf, all_undfs)
  output_field_fs => null()

  ! Set up arrays for AllGather

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))

  all_undfs = 0

  ! Calculate the nodal coords for a field on Wtheta
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, Wtheta )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do

  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == geometry_spherical ) then
     call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if

  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the unique fractional levels to set up full levels vertical output domain
  fractional_levels_full_faces => proxy_coord_output(1)%vspace%get_levels()

  nfull_levels = size(fractional_levels_full_faces)

  call xios_set_domain_attr("face_full_levels", ni_glo=global_undf,          &
                            ibegin=ibegin_faces,                             &
                            ni=local_undf(1)/nfull_levels,                   &
                            type='unstructured')
  call xios_set_domain_attr("face_full_levels", lonvalue_1d=faces_lon,       &
                            latvalue_1d=faces_lat)
  call xios_set_domain_attr("face_full_levels", bounds_lon_1d=bnd_faces_lon, &
                            bounds_lat_1d=bnd_faces_lat)

  ! Pass local portion of domain_index
  call xios_set_domain_attr("face_full_levels", &
                            i_index=int(domain_index(1:(local_undf(1)/nfull_levels))))


  ! Clean up things ready to reuse for edge domain setup
  deallocate(local_undf, all_undfs)
  output_field_fs => null()

  ! Set up arrays for AllGather

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))
  
  all_undfs = 0

  ! Set up the edge domain
  ! Here we use information from W2H to calculate the physical coordinates
  ! and the 'half levels' vertical domain

  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W2H )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do


  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == geometry_spherical ) then
     call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if


  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the unique fractional levels to set up half levels vertical output domain
  fractional_levels_half_edges => proxy_coord_output(1)%vspace%get_levels()

  nhalf_levels = size(fractional_levels_half_edges)

  ! Allocate coordinate arrays for edges

  allocate(edges_lon( size_w2h))
  allocate(edges_lat( size_w2h))

  edges_lon =  proxy_coord_output(1)%data(1: local_undf(1):nhalf_levels) * r2d
  edges_lat =  proxy_coord_output(2)%data(1: local_undf(1):nhalf_levels) * r2d


 !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call all_gather ( local_undf, all_undfs, 1 )

  ! Adjust size of data taking into account how many levels we have (same for each
  ! partition as we only partition horizontally)

  all_undfs = all_undfs/nhalf_levels

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for xios face domain

  global_undf = sum(all_undfs)


  ! Calculate ibegin for each rank as we have the array of undfs in order
  ! we can just sum to get it.

  if (get_comm_rank() == 0) then
    ibegin_edges = 0
  else
    ibegin_edges = sum(all_undfs(1:get_comm_rank()))
  end if


  call xios_set_domain_attr("edge_half_levels", ni_glo=global_undf,          &
                            ibegin=ibegin_edges,                             &
                            ni=local_undf(1)/nhalf_levels,                   &
                            type='unstructured')
  call xios_set_domain_attr("edge_half_levels", lonvalue_1d=edges_lon,       &
                            latvalue_1d=edges_lat)
  call xios_set_domain_attr("edge_half_levels", bounds_lon_1d=bnd_edges_lon, &
                            bounds_lat_1d=bnd_edges_lat)


  ! Clean up things that are not needed after domain setup
  deallocate(local_undf, all_undfs)
  deallocate(nodes_lat, nodes_lon, bnd_nodes_lat, bnd_nodes_lon)
  deallocate(faces_lat, faces_lon, bnd_faces_lat, bnd_faces_lon)
  deallocate(bnd_edges_lat, bnd_edges_lon)
  deallocate(nodes_lat_full, nodes_lon_full)
  deallocate(domain_index)
  fractional_levels_half_faces => null()
  fractional_levels_full_faces => null()
  fractional_levels_half_edges => null()
  output_field_fs => null()
  w2h_fs => null()

  return
end subroutine xios_diagnostic_domain_init

!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS checkpoint domain initialisation
!!
!!  @details  Performs checkpoint domain init and returns
!!            Assumes an unstructured 1D domain type
!!
!!  @param[in]      fs_id         Function space id
!!  @param[in]      domain_name   XIOS domain name
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      chi           Coordinate field
!!  @param[in]      use_index     Flag to specify use of domain index
!!                                to preserve order over decomposition

!-------------------------------------------------------------------------------

subroutine xios_checkpoint_domain_init(fs_id, domain_name, mesh_id, chi, use_index)

  implicit none

  ! Arguments

  integer(i_def),   intent(in)         :: fs_id
  character(len=*), intent(in)         :: domain_name
  integer(i_def),   intent(in)         :: mesh_id
  type(field_type), intent(in)         :: chi(3)
  logical(l_def),   intent(in)         :: use_index


  ! Local variables

  integer(i_def)    :: i

  ! Checkpoint domain
  integer(i_def)                      :: ibegin_checkpoint
  real(dp_xios),allocatable           :: checkpoint_lon(:)
  real(dp_xios),allocatable           :: checkpoint_lat(:)
  real(dp_xios),allocatable           :: bnd_checkpoint_lon(:,:)
  real(dp_xios),allocatable           :: bnd_checkpoint_lat(:,:)
  integer(i_halo_index),allocatable   :: domain_index(:)


  ! Variables needed to compute output domain coordinates in lat-long

  ! Transformed coords for nodal output
  type( field_type ) :: coord_output(3)
  ! Field proxies (to calculate domain coordinate info)
  type(field_proxy_type), target  :: proxy_coord_output(3)


  type(function_space_type), pointer :: output_field_fs   => null()

  ! Variables for the gather to determine global domain sizes
  ! from the local partitioned ones

  integer(i_def)                :: global_undf_checkpoint
  integer(i_def), allocatable   :: local_undf(:)
  integer(i_def), allocatable   :: all_undfs_checkpoint_domain(:)

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def) :: r2d


  if ( geometry == geometry_spherical ) then
   r2d = radians_to_degrees
  else
   r2d = 1.0_r_def
  endif


  ! Set up arrays to hold number of dofs for local and global domains


  allocate(local_undf(1))
  allocate(all_undfs_checkpoint_domain(get_comm_size()))

  all_undfs_checkpoint_domain = 0

  ! Create appropriate function space in order to be able to get the
  ! physical coordinates

  output_field_fs => function_space_collection%get_fs( mesh_id,       &
                                                       element_order, &
                                                       fs_id)

  ! Calculate the nodal coords for a field on the function space

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do


  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == geometry_spherical ) then
    call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if


  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()


  !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!

  call all_gather ( local_undf, all_undfs_checkpoint_domain, 1 )

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for checkpoint domain

  global_undf_checkpoint = sum(all_undfs_checkpoint_domain)


  ! Calculate ibegin for each rank as we have the array of undfs in order
  ! we can just sum to get it.

  if (get_comm_rank() == 0) then
    ibegin_checkpoint = 0
  else
    ibegin_checkpoint = sum(all_undfs_checkpoint_domain(1:get_comm_rank()))
  end if

  ! Allocate coordinate arrays to be the size required for checkpoint domain.
  ! Essentially up to last owned dof of the current partition.

  allocate( checkpoint_lon( size( proxy_coord_output(1)%data(1: local_undf(1)))) )
  allocate( checkpoint_lat( size( proxy_coord_output(2)%data(1: local_undf(1)))) )

  ! Populate the arrays with data
  checkpoint_lon =  proxy_coord_output(1)%data(1: local_undf(1)) * r2d
  checkpoint_lat =  proxy_coord_output(2)%data(1: local_undf(1)) * r2d

  allocate(bnd_checkpoint_lon(1,size(checkpoint_lon)))
  allocate(bnd_checkpoint_lat(1,size(checkpoint_lat)))

  ! Construct bounds arrays
  bnd_checkpoint_lon=(reshape(checkpoint_lon, (/1, size(checkpoint_lon)/) ) )
  bnd_checkpoint_lat=(reshape(checkpoint_lat, (/1, size(checkpoint_lat)/) ) )


  call xios_set_domain_attr(trim(domain_name), ni_glo=global_undf_checkpoint,    &
                            ibegin=ibegin_checkpoint, ni=local_undf(1),          &
                            type='unstructured')
  call xios_set_domain_attr(trim(domain_name), lonvalue_1d=checkpoint_lon,       &
                            latvalue_1d=checkpoint_lat)
  call xios_set_domain_attr(trim(domain_name), bounds_lon_1d=bnd_checkpoint_lon, &
                            bounds_lat_1d=bnd_checkpoint_lat)

  ! If we have requested to use domain index then get it and use it

  if (use_index) then

    ! Allocate domain_index - it is of size ndof_glob
    allocate(domain_index(output_field_fs%get_ndof_glob()))

    ! Populate domain_index for this rank 
    call output_field_fs%get_global_dof_id(domain_index)

    ! Pass local portion of domain_index (up to undf)
    call xios_set_domain_attr(domain_name, i_index=int(domain_index(1:local_undf(1))))

  end if 


  if ( allocated(checkpoint_lon) )     deallocate(checkpoint_lon)
  if ( allocated(checkpoint_lat) )     deallocate(checkpoint_lat)
  if ( allocated(domain_index) )    deallocate(domain_index)
  if ( allocated(bnd_checkpoint_lon) ) deallocate(bnd_checkpoint_lon)
  if ( allocated(bnd_checkpoint_lat) ) deallocate(bnd_checkpoint_lat)
  if ( allocated(local_undf) )      deallocate(local_undf)
  if ( allocated(all_undfs_checkpoint_domain) ) deallocate(all_undfs_checkpoint_domain)
  nullify( output_field_fs )

  return
end subroutine xios_checkpoint_domain_init



!> @brief   Compute the node domain coords for this partition
!> @details Samples the chi field at nodal points, calculates cartesian coordinates.
!>          For spherical geometry, converts to lat-lon in degrees for specified layer
!>@param[in] local_mesh the id of the partitioned mesh 
!>@param[in] nodal_coords input field
!>@param[in] chi input coordinate field
!>@param[in] nlayers the number of layers data is output on
!>@param[in] ncells the number of cells on the partition
!>@param[out] lon_coords array of longitude coordinates for the nodes
!>@param[out] lat_coords array of latitude coordinates for the nodes
!>@param[inout] face_bnds_lon_coords array of longitude coords making up the faces 
!>@param[inout] face_bnds_lat_coords array of latitude coords making up the faces 
!>@param[inout] edge_bnds_lon_coords array of coords making up the edges 
!>@param[inout] edge_bnds_lat_coords array of coords making up the edges 

subroutine calc_xios_domain_coords(local_mesh, nodal_coords, chi, &
                                   nlayers, ncells,               &
                                   lon_coords, lat_coords,        &
                                   face_bnds_lon_coords,          &
                                   face_bnds_lat_coords,          &
                                   edge_bnds_lon_coords,          &
                                   edge_bnds_lat_coords)

  implicit none

  type(mesh_type), pointer, intent(in) :: local_mesh
  type(field_type),   intent(in)       :: nodal_coords(3)
  type(field_type),   intent(in)       :: chi(:)
  integer(i_def),     intent(in)       :: nlayers
  integer(i_def),     intent(in)       :: ncells
  real(kind=r_def),   intent(out)      :: lon_coords(:), lat_coords(:)
  real(kind=dp_xios), intent(inout)    :: face_bnds_lon_coords(:,:)
  real(kind=dp_xios), intent(inout)    :: face_bnds_lat_coords(:,:)
  real(kind=dp_xios), intent(inout)    :: edge_bnds_lon_coords(:,:)
  real(kind=dp_xios), intent(inout)    :: edge_bnds_lat_coords(:,:)

  type(field_proxy_type) :: x_p(3), chi_p(3)
   
  integer(i_def)            :: cell, edge_count
  integer(i_def)            :: ndf_chi, ndf_x
  integer(i_def)            :: dim_chi
  integer, pointer          :: map_chi(:)   => null()
  integer, pointer          :: map_x(:)     => null()
  real(kind=r_def), pointer :: nodes_x(:,:) => null()
  real(kind=r_def)          :: xyz(3)
  real(kind=r_def)          :: llr(3)

  real(kind=r_def), allocatable  :: basis_chi(:,:,:)
  integer(i_def)                 :: df_x, df_chi, i
  integer(i_def)                 :: edge1, edge2

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def) :: r2d

  edge_count = 0

  do i = 1,3
    x_p(i)   = nodal_coords(i)%get_proxy()
    chi_p(i) = chi(i)%get_proxy()
  end do

  ndf_x  = x_p(1)%vspace%get_ndf( )
  nodes_x => x_p(1)%vspace%get_nodes()
  ndf_chi  = chi_p(1)%vspace%get_ndf( )

  dim_chi = chi_p(1)%vspace%get_dim_space( )

  allocate(basis_chi(dim_chi, ndf_chi, ndf_x))

  do df_x = 1, ndf_x
    do df_chi = 1, ndf_chi
      basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
    end do
  end do

  if (chi_p(1)%is_dirty(depth=1)) then
    call chi_p(1)%halo_exchange(depth=1)
  end if

  if (chi_p(2)%is_dirty(depth=1)) then
    call chi_p(2)%halo_exchange(depth=1)
  end if

  if (chi_p(3)%is_dirty(depth=1)) then
    call chi_p(3)%halo_exchange(depth=1)
  end if

  ! Loop over cells
  do cell = 1, ncells

    map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
    map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )

    ! Loop over bottom half of the cell dofmap for the given layer
    do df_x = 1,(ndf_x/2)
      xyz(:) = 0.0_r_def
      do df_chi = 1, (ndf_chi/2)
        xyz(1) = xyz(1) + chi_p(1)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi_p(2)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi_p(3)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
      end do

      ! Convert to lat-lon in degrees if required
      if ( geometry == geometry_spherical ) then

        r2d = radians_to_degrees

        call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))

        lon_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(1)*r2d
        lat_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(2)*r2d

        face_bnds_lon_coords(df_x,cell) = llr(1)*r2d
        face_bnds_lat_coords(df_x,cell) = llr(2)*r2d
      else
        r2d = 1.0_r_def

        lon_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(1)*r2d
        lat_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(2)*r2d

        face_bnds_lon_coords(df_x,cell) = xyz(1)*r2d
        face_bnds_lat_coords(df_x,cell) = xyz(2)*r2d

      endif

    end do ! Loop over bottom layer dofs

    ! For this cell compute the edge-bounds coordinates from the face-bounds coordinates

    do df_x = 1,(ndf_x/2)
      
      ! Retrieve the lat / lon coords of the points bounding the edge
      if (df_x == ndf_x/2) then
        edge1 = df_x
        edge2 = 1
      else
        edge1 = df_x
        edge2 = df_x + 1
      endif

      ! Is the edge owned by this cell? 
      if (local_mesh%get_edge_cell_owner(df_x, cell) == cell) then

        edge_count = edge_count + 1

        edge_bnds_lon_coords(1,edge_count) = face_bnds_lon_coords(edge1,cell)
        edge_bnds_lon_coords(2,edge_count) = face_bnds_lon_coords(edge2,cell)
        edge_bnds_lat_coords(1,edge_count) = face_bnds_lat_coords(edge1,cell)
        edge_bnds_lat_coords(2,edge_count) = face_bnds_lat_coords(edge2,cell)


      end if ! Edge is owned by this cell 

    end do ! loop over edges


  end do

  call x_p(1)%set_dirty()
  call x_p(2)%set_dirty()
  call x_p(3)%set_dirty()

  deallocate(basis_chi)

  nullify( map_chi, map_x, nodes_x )

end subroutine calc_xios_domain_coords

!> @brief   Function to determine output filename at a given timestep
!>@param[in] stem_name string file stem
!>@param[in] file_type string used to identify file type (e.g. 'nodal')
!>@param[in] field_name name of the field
!>@param[in] ts time step
!>@param[in] file extension

function ts_fname(stem_name, file_type, field_name, ts, ext)

  implicit none

  character(len=*),    intent(in) :: field_name, stem_name, &
                                     file_type, ext
  integer(i_def),      intent(in) :: ts
  character(len=str_max_filename) :: rank_name
  character(len=str_max_filename) :: ts_fname
  integer(i_def)                  :: total_ranks
  integer(i_def)                  :: local_rank

  total_ranks = get_comm_size()
  local_rank = get_comm_rank()

  if( total_ranks == 1 )then
      rank_name=ext
    else
      write(rank_name,"("".Rank"",I6.6)")local_rank
  end if

  write(ts_fname,'(A,A,A,A,A,I6.6,A)') trim(stem_name),"_", &
         trim(file_type),trim(field_name),"_T",ts,trim(rank_name)//ext

end function ts_fname

!> @brief   Output a field in nodal format to text file
!>@param[in] nodal_coordinates field holding coordinate information
!>@param[in] level field holding level information
!>@param[in] nodal_output field holding diagnostic data
!>@param[in] fspace_dimension dimension of the field's function space
!>@param[in] output_unit file unit to write to
!>@param[in] fname file name to write to
!------------------------------------------------------------------------------- 
subroutine nodal_write_field(nodal_coordinates, level, nodal_output, &
                             fspace_dimension, output_unit, fname)

  implicit none
  type(field_type), intent(in) :: nodal_coordinates(3), level, nodal_output(3)
  integer(i_def),   intent(in) :: fspace_dimension, output_unit
  character(str_max_filename), intent(in) :: fname
  type(field_proxy_type) :: x_p(3), l_p, n_p(3)

  integer(i_def) :: df, undf, i


  do i = 1,3
    x_p(i) = nodal_coordinates(i)%get_proxy()
    n_p(i) = nodal_output(i)%get_proxy()
  end do
  l_p = level%get_proxy()

  undf = n_p(1)%vspace%get_last_dof_owned()

  open(OUTPUT_UNIT, file = trim(fname), status = "replace")   
  write(OUTPUT_UNIT,'(A)') 'x = [' 
  if ( fspace_dimension  == 1 ) then
    do df = 1,undf
      write(OUTPUT_UNIT,'(5e25.15e3)') x_p(1)%data(df), x_p(2)%data(df), &
                          x_p(3)%data(df), l_p%data(df), n_p(1)%data(df)
    end do
  else
    do df = 1,undf
      write(OUTPUT_UNIT,'(7e25.15e3)') x_p(1)%data(df), x_p(2)%data(df), &
                                       x_p(3)%data(df), l_p%data(df),    &
                                       n_p(1)%data(df), n_p(2)%data(df), &
                                       n_p(3)%data(df)
    end do
  end if
  write(OUTPUT_UNIT,'(A)') '];'
  close(OUTPUT_UNIT)

end subroutine nodal_write_field


!> @brief   I/O handler for reading a netcdf checkpoint
!> @details Legacy method for reading checkpoints
!           Note this routine accepts a field name but
!           doesn't use it - this is to keep the interface
!           the same for all methods 
!>@param[in] field_name Name of the field to read
!>@param[in] file_name Name of the file to read from
!>@param[in,out] field_proxy the proxy of the field to read data into
subroutine checkpoint_read_netcdf(field_name, file_name, field_proxy)
  use field_io_ncdf_mod,    only : field_io_ncdf_type

  implicit none

  character(len=*), intent(in) :: field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(inout) :: field_proxy

  type(field_io_ncdf_type), allocatable :: ncdf_file

  allocate(ncdf_file)

  call ncdf_file%file_open( file_name )

  call ncdf_file%read_field_data ( field_proxy%data(:) )

  call ncdf_file%file_close()

  deallocate(ncdf_file)


end subroutine checkpoint_read_netcdf

!> @brief   I/O handler for writing a netcdf checkpoint
!> @details Legacy method for writing checkpoints
!           Note this routine accepts a field name but
!           doesn't use it - this is to keep the interface
!           the same for all methods 
!>@param[in] field_name Name of the field to write
!>@param[in] file_name Name of the file to write to
!>@param[in,out] field_proxy the proxy of the field to write
subroutine checkpoint_write_netcdf(field_name, file_name, field_proxy)
  use field_io_ncdf_mod,    only : field_io_ncdf_type

  implicit none

  character(len=*), intent(in) :: field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(in) :: field_proxy

  type(field_io_ncdf_type), allocatable :: ncdf_file

  allocate(ncdf_file)

  call ncdf_file%file_new( file_name )

  call ncdf_file%write_field_data ( field_proxy%data(:) )
  call ncdf_file%file_close()

  deallocate(ncdf_file)


end subroutine checkpoint_write_netcdf


!> @brief   I/O handler for reading an XIOS netcdf checkpoint
!> @details Note this routine accepts a filename but doesn't
!           use it - this is to keep the interface the same
!           for all methods
!>@param[in] xios_field_name XIOS unique id for the field
!>@param[in] file_name Name of the file to read
!>@param[in,out] field_proxy the proxy of the field to read into
subroutine checkpoint_read_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: undf


  ! We only read in up to undf for the partition
  undf = field_proxy%vspace%get_last_dof_owned()

  call xios_recv_field(xios_field_name, field_proxy%data(1:undf))


end subroutine checkpoint_read_xios


!> @brief   I/O handler for writing an XIOS netcdf checkpoint
!> @details Note this routine accepts a filename but doesn't
!           use it - this is to keep the interface the same
!           for all methods
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] file_name Name of the file to write
!>@param[in,out] field_proxy the proxy of the field to write
subroutine checkpoint_write_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: undf

  undf = field_proxy%vspace%get_last_dof_owned()

  call xios_send_field(xios_field_name, field_proxy%data(1:undf))


end subroutine checkpoint_write_xios


!> @brief   Write a checkpoint from a collection of fields
!> @details Iterate over a field collection and checkpoint each field
!>          if it is enabled for checkpointing
!>@param[in] state - a collection of fields to checkpoint
!>@param[in] timestep the current timestep
subroutine write_checkpoint(state, timestep)

    implicit none

    type( field_collection_type ), intent(inout) :: state
    integer(i_def), intent(in) :: timestep

    type( field_collection_iterator_type) :: iter
    type( field_type ), pointer :: fld => null()

    iter=state%get_iterator()
    do
      if(.not.iter%has_next())exit
      fld=>iter%next()
      if (fld%can_checkpoint()) then
        write(log_scratch_space,'(3A,I6)') &
            "Checkpointing ", trim(adjustl(fld%get_name())), " at timestep ", &
            timestep
        call log_event(log_scratch_space,LOG_LEVEL_INFO)
        call fld%write_checkpoint("checkpoint_"//trim(adjustl(fld%get_name())), &
                                  trim(ts_fname(checkpoint_stem_name,&
                                  "", trim(adjustl(fld%get_name())),timestep,"")))
      else

        call log_event( 'Checkpointing for  '// trim(adjustl(fld%get_name())) // &
                        ' not set up', LOG_LEVEL_INFO )

      end if

    end do

    nullify(fld)

end subroutine write_checkpoint

!> @brief   Read from a checkpoint into a collection of fields
!> @details Iterate over a field collection and read each field
!>          into a collection, if it is enabled for checkpointing
!>@param[in] state -  the collection of fields to populate
!>@param[in] timestep the current timestep
subroutine read_checkpoint(state, timestep)

    implicit none

    type( field_collection_type ), intent(inout) :: state
    integer(i_def), intent(in) :: timestep

    type( field_collection_iterator_type) :: iter
    type( field_type ), pointer :: fld => null()

    iter=state%get_iterator()
    do
      if(.not.iter%has_next())exit
      fld=>iter%next()
      if (fld%can_checkpoint()) then
        call log_event( &
          'Reading checkpoint file to restart '//trim(adjustl(fld%get_name())), &
          LOG_LEVEL_INFO)
        call fld%read_checkpoint("restart_"//trim(adjustl(fld%get_name())), &
                              trim(ts_fname(checkpoint_stem_name, &
                              "", trim(adjustl(fld%get_name())),timestep,"")))
      else
        call log_event( 'Checkpointing for  '// trim(adjustl(fld%get_name())) // &
                        ' not set up', LOG_LEVEL_INFO )
      end if
    end do

    nullify(fld)

end subroutine read_checkpoint

! The dump read and write subroutines are wrappers that call 
! XIOS read and write I/O subroutines

!> @brief   Write a field to a dump via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
subroutine dump_write_xios(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  call xios_write_field_face(xios_field_name, field_proxy)

end subroutine dump_write_xios

!> @brief   Read a field from a dump via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in,out] field_proxy a field proxy containing the data to output
subroutine dump_read_xios(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(inout) :: field_proxy

  call xios_read_field_face("read_"//xios_field_name, field_proxy)

end subroutine dump_read_xios

!> @brief   Write a collection of fields
!> @details Iterate over a field collection and write each field
!>          if it is enabled for write
!>@param[in] state - a collection of fields
subroutine write_state(state)

    implicit none

    type( field_collection_type ), intent(inout) :: state

    type( field_collection_iterator_type) :: iter
    type( field_type ), pointer :: fld => null()

    iter=state%get_iterator()
    do
      if(.not.iter%has_next())exit
      fld=>iter%next()
      if (fld%can_write()) then
        write(log_scratch_space,'(3A,I6)') &
            "Writing ", trim(adjustl(fld%get_name()))
        call log_event(log_scratch_space,LOG_LEVEL_INFO)
        call fld%write_field(trim(adjustl(fld%get_name())))
      else

        call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                        ' not set up', LOG_LEVEL_INFO )

      end if

    end do

    nullify(fld)

end subroutine write_state

!> @brief   Read into a collection of fields
!> @details Iterate over a field collection and read each field
!>          into a collection, if it is enabled for read
!>@param[in,out] state -  the collection of fields to populate
subroutine read_state(state)

    implicit none

    type( field_collection_type ), intent(inout) :: state

    type( field_collection_iterator_type) :: iter
    type( field_type ), pointer :: fld => null()

    iter=state%get_iterator()
    do
      if(.not.iter%has_next())exit
      fld=>iter%next()
      if (fld%can_read()) then
        call log_event( &
          'Reading '//trim(adjustl(fld%get_name())), &
          LOG_LEVEL_INFO)
        call fld%read_field(trim(adjustl(fld%get_name())))
      else
        call log_event( 'Read method for  '// trim(adjustl(fld%get_name())) // &
                        ' not set up', LOG_LEVEL_INFO )
      end if
    end do

    nullify(fld)

end subroutine read_state


!> @brief   Output a field in UGRID format on the node domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!------------------------------------------------------------------------------- 
subroutine xios_write_field_node(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)


  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal domain size for the rank
  call xios_get_domain_attr('node', ni=domain_size)
  ! Get the expected vertical axis size
  call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)


  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))


  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change 

  ! First get the data on the same level ordered in chunks  
  do i=0, axis_size-1
    send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
               field_proxy%data(i+1:undf:axis_size)
  end do

  ! Reshape into 2D horizontal + vertical levels for output
  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine xios_write_field_node

!> @brief   Output a single level field in UGRID format on the face domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!-------------------------------------------------------------------------------
subroutine xios_write_field_single_face(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: undf
  integer(i_def) :: domain_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal size
  ! all 2D fields are nominally in W3, hence half levels
  call xios_get_domain_attr('face_half_levels', ni=domain_size)

  ! Size the arrays to be what is expected
  allocate(send_field(domain_size))

  ! All data are scalar fields

  send_field(1:domain_size) = field_proxy%data(1:undf)

  call xios_send_field(xios_field_name, send_field)

  deallocate(send_field)

end subroutine xios_write_field_single_face

!> @brief   Output a field in UGRID format on the face domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!------------------------------------------------------------------------------- 
subroutine xios_write_field_face(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! Get the expected horizontal and vertical axis size

  if (fs_id == W3) then
    call xios_get_domain_attr('face_half_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)
  else
    call xios_get_domain_attr('face_full_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)
  end if


  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))


  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change 

  ! First get the data on the same level ordered in chunks  
  do i=0, axis_size-1
    send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
               field_proxy%data(i+1:undf:axis_size)
  end do

  ! Reshape into 2D horizontal + vertical levels for output

  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine xios_write_field_face

!> @brief   Read a field in UGRID format on the face domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy to read data into
subroutine xios_read_field_face(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: recv_field(:)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! get the horizontal / vertical domain sizes

  if (fs_id == W3) then
    call xios_get_domain_attr('face_half_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)
  else
    call xios_get_domain_attr('face_full_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)
  end if


  ! Size the array to be what is expected
  allocate(recv_field(domain_size*axis_size))

  ! Read the data into a temporary array - this should be in the correct order
  ! as long as we set up the horizontal domain using the global index
  call xios_recv_field(xios_field_name, recv_field)

  ! We need to reshape the incoming data to get the correct data layout for the LFRic
  ! field - It is the reverse of what we did for writing

  do i=0, axis_size-1

    field_proxy%data(i+1:undf:axis_size) = &
                                   recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size)   
  end do

  deallocate(recv_field)

end subroutine xios_read_field_face

!> @brief   Read a single level field in UGRID format on the face domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy to read data into
subroutine xios_read_field_single_face(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size
  real(dp_xios), allocatable :: recv_field(:)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()


  ! Get the expected horizontal size
  ! all 2D fields are nominally in W3, hence half levels
  call xios_get_domain_attr('face_half_levels', ni=domain_size)


  ! Size the array to be what is expected
  allocate(recv_field(domain_size))

  ! Read the data into a temporary array - this should be in the correct order
  ! as long as we set up the horizontal domain using the global index
  call xios_recv_field(xios_field_name, recv_field)

  field_proxy%data(1:undf) = recv_field(1:domain_size)

  deallocate(recv_field)


end subroutine xios_read_field_single_face



!> @brief   Output a field in UGRID format on the edge domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!------------------------------------------------------------------------------- 
subroutine xios_write_field_edge(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()
  
  ! Get the expected horizontal and vertical axis size

  call xios_get_domain_attr('edge_half_levels', ni=domain_size)
  call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)

  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))

  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change 

  ! First get the data on the same level ordered in chunks  
  do i=0, axis_size-1
    send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
               field_proxy%data(i+1:undf:axis_size)
  end do

  ! Reshape into 2D horizontal + vertical levels for output

  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)
  
end subroutine xios_write_field_edge

end module io_mod

