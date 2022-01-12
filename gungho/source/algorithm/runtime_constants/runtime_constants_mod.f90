!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief A module that controls set-up of various run time constants.
!>
!> @details This module controls the set-up of various objects that are
!>          created at setup and are not changed thereafter but are needed
!>          throughout the algorithm layers.
module runtime_constants_mod

  use boundaries_config_mod,             only: limited_area
  use constants_mod,                     only: i_def, r_def, str_def
  use field_mod,                         only: field_type
  use formulation_config_mod,            only: moisture_conservation
  use io_config_mod,                     only: subroutine_timers
  use log_mod,                           only: log_event, LOG_LEVEL_INFO
  use runtime_tools_mod,                 only: primary_mesh_label,      &
                                               shifted_mesh_label,      &
                                               double_level_mesh_label, &
                                               twod_mesh_label,         &
                                               multigrid_mesh_label,    &
                                               extra_mesh_label
  use timer_mod,                         only: timer

  implicit none

  private

  integer(i_def), target :: global_shifted_mesh_id
  integer(i_def), target :: global_double_level_mesh_id

  ! Public functions to create and access the module contents
  public :: create_runtime_constants
  public :: final_runtime_constants

contains
  !>@brief Subroutine to create the runtime constants
  !> @param[in] mesh_id              Mesh_id
  !> @param[in] twod_mesh_id         Mesh_id for 2D domain
  !> @param[in] chi                  Coordinate field on primary mesh
  !> @param[in] panel_id             panel id
  !> @param[in] dt                   The model timestep length
  !> @param[in] shifted_mesh_id      Mesh_id for vertically shifted field
  !> @param[in] shifted_chi          Coordinate field for vertically shifted field
  !> @param[in] double_level_mesh_id Mesh_id for double level field
  !> @param[in] double_level_chi     Coordinate field for double level field
  !> @param[in] mg_mesh_ids          A list of mesh IDs for the multigrid meshes
  !> @param[in] mg_2D_mesh_ids       A list of mesh IDs for the 2D MG meshes
  !> @param[in] chi_mg               The coordinate fields for the MG meshes
  !> @param[in] panel_id_mg          The panel_id fields for the MG meshes
  !> @param[in] extra_mesh_ids       A list of mesh IDs for any extra meshes
  !> @param[in] extra_2D_mesh_ids    A list of mesh IDs for extra 2D meshes
  !> @param[in] chi_extra            The coordinate fields for any extra meshes
  !> @param[in] panel_id_extra       The panel_id fields for any extra MG meshes
  subroutine create_runtime_constants(mesh_id, twod_mesh_id, &
                                      chi,                   &
                                      panel_id, dt,          &
                                      shifted_mesh_id,       &
                                      shifted_chi,           &
                                      double_level_mesh_id,  &
                                      double_level_chi,      &
                                      mg_mesh_ids,           &
                                      mg_2D_mesh_ids,        &
                                      chi_mg,                &
                                      panel_id_mg,           &
                                      extra_mesh_ids,        &
                                      extra_2D_mesh_ids,     &
                                      chi_extra,             &
                                      panel_id_extra         )

    ! Other runtime_constants modules
    use advective_update_alg_mod,    only: advective_update_set_num_meshes
    use fem_constants_mod,           only: create_fem_constants
    use flux_alg_mod,                only: flux_alg_set_num_meshes
    use geometric_constants_mod,     only: create_geometric_constants
    use intermesh_constants_mod,     only: create_intermesh_constants
    use limited_area_constants_mod,  only: create_limited_area_constants
    use physical_op_constants_mod,   only: create_physical_op_constants
    use rk_transport_rho_mod,        only: rk_transport_rho_set_num_meshes
    use rk_transport_theta_mod,      only: rk_transport_theta_set_num_meshes
    use runge_kutta_init_mod,        only: runge_kutta_init
    use runtime_tools_mod,           only: init_mesh_id_list, &
                                           init_hierarchical_mesh_id_list

    implicit none

    integer(kind=i_def),             intent(in) :: mesh_id, twod_mesh_id
    type(field_type),      target,   intent(in) :: chi(:)
    type(field_type),      target,   intent(in) :: panel_id
    real(r_def),                     intent(in) :: dt
    integer(kind=i_def),   optional, intent(in) :: shifted_mesh_id
    type(field_type),      optional, intent(in) :: shifted_chi(:)
    integer(kind=i_def),   optional, intent(in) :: double_level_mesh_id
    type(field_type),      optional, intent(in) :: double_level_chi(:)
    integer(kind=i_def),   optional, intent(in) :: mg_mesh_ids(:)
    integer(kind=i_def),   optional, intent(in) :: mg_2D_mesh_ids(:)
    type(field_type),      optional, intent(in) :: chi_mg(:,:)
    type(field_type),      optional, intent(in) :: panel_id_mg(:)
    integer(kind=i_def),   optional, intent(in) :: extra_mesh_ids(:)
    integer(kind=i_def),   optional, intent(in) :: extra_2D_mesh_ids(:)
    type(field_type),      optional, intent(in) :: chi_extra(:,:)
    type(field_type),      optional, intent(in) :: panel_id_extra(:)

    ! Internal variables
    integer(kind=i_def)                         :: num_meshes, mesh_counter, i, j
    integer(kind=i_def),            allocatable :: mesh_id_list(:)
    integer(kind=i_def),            allocatable :: label_list(:)
    type(field_type),               allocatable :: chi_list(:,:)
    type(field_type),               allocatable :: panel_id_list(:)

    if ( subroutine_timers ) call timer('runtime_constants_alg')
    call log_event( "Gungho: creating runtime_constants", LOG_LEVEL_INFO )

    !==========================================================================!
    ! Turn all the mesh IDs and coordinate fields into lists
    !==========================================================================!

    ! Count the number of meshes that we have
    num_meshes = 2_i_def ! We should always have primary mesh_id and twod_mesh_id
    if ( present(shifted_mesh_id) .and. present(shifted_chi) ) num_meshes = num_meshes + 1_i_def
    if ( present(double_level_mesh_id) .and. present(double_level_chi) ) num_meshes = num_meshes + 1_i_def
    if ( present(mg_mesh_ids) .and. present(chi_mg) ) num_meshes = num_meshes + size(mg_mesh_ids)
    if ( present(mg_2D_mesh_ids) .and. present(chi_mg) ) num_meshes = num_meshes + size(mg_2D_mesh_ids)
    if ( present(extra_mesh_ids) .and. present(chi_extra) ) num_meshes = num_meshes + size(extra_mesh_ids)
    if ( present(extra_2D_mesh_ids) .and. present(chi_extra) ) num_meshes = num_meshes + size(extra_2D_mesh_ids)

    allocate(mesh_id_list(num_meshes))
    allocate(chi_list(3,num_meshes))
    allocate(panel_id_list(num_meshes))
    allocate(label_list(num_meshes))

    ! Populate these lists
    mesh_counter = 1_i_def
    label_list(mesh_counter) = primary_mesh_label
    mesh_id_list(mesh_counter) = mesh_id
    call panel_id%copy_field(panel_id_list(mesh_counter))
    do j = 1, 3
      call chi(j)%copy_field(chi_list(j, mesh_counter))
    end do

    ! Primary 2D mesh
    mesh_counter = mesh_counter + 1_i_def
    label_list(mesh_counter) = twod_mesh_label
    mesh_id_list(mesh_counter) = twod_mesh_id
    call panel_id%copy_field(panel_id_list(mesh_counter))
    do j = 1, 3
      call chi(j)%copy_field(chi_list(j, mesh_counter))
    end do

    if ( present(shifted_mesh_id) .and. present(shifted_chi) ) then
      global_shifted_mesh_id = shifted_mesh_id
      mesh_counter = mesh_counter + 1_i_def
      label_list(mesh_counter) = shifted_mesh_label
      mesh_id_list(mesh_counter) = shifted_mesh_id
      call panel_id%copy_field(panel_id_list(mesh_counter)) ! Same as for primary mesh
      do j = 1, 3
        call shifted_chi(j)%copy_field(chi_list(j, mesh_counter))
      end do
    end if

    if ( present(double_level_mesh_id) .and. present(double_level_chi) ) then
      global_double_level_mesh_id = double_level_mesh_id
      mesh_counter = mesh_counter + 1_i_def
      label_list(mesh_counter) = double_level_mesh_label
      mesh_id_list(mesh_counter) = double_level_mesh_id
      call panel_id%copy_field(panel_id_list(mesh_counter)) ! Same as for primary mesh
      do j = 1, 3
        call double_level_chi(j)%copy_field(chi_list(j, mesh_counter))
      end do
    end if

    if ( present(mg_mesh_ids) .and. present(chi_mg) ) then
      do i = 1, size(mg_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = multigrid_mesh_label
        mesh_id_list(mesh_counter) = mg_mesh_ids(i)
        call panel_id_mg(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          call chi_mg(j,i)%copy_field(chi_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(mg_2D_mesh_ids) .and. present(chi_mg) ) then
      do i = 1, size(mg_2D_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = twod_mesh_label
        mesh_id_list(mesh_counter) = mg_2D_mesh_ids(i)
        call panel_id_mg(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          call chi_mg(j,i)%copy_field(chi_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(extra_mesh_ids) .and. present(chi_extra) ) then
      do i = 1, size(extra_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = extra_mesh_label
        mesh_id_list(mesh_counter) = extra_mesh_ids(i)
        call panel_id_extra(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          call chi_extra(j,i)%copy_field(chi_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(extra_2D_mesh_ids) .and. present(chi_extra) ) then
      do i = 1, size(extra_2D_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = twod_mesh_label
        mesh_id_list(mesh_counter) = extra_2D_mesh_ids(i)
        call panel_id_extra(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          call chi_extra(j,i)%copy_field(chi_list(j, mesh_counter))
        end do
      end do
    end if

    !==========================================================================!
    ! Set up runtime_constants for each category
    !==========================================================================!

    call init_mesh_id_list(mesh_id_list)

    if ( present(mg_mesh_ids) ) then
      ! mg_mesh_ids contains all mesh ids used in the multigrid chain
      ! including the primary mesh
      call init_hierarchical_mesh_id_list(mg_mesh_ids)
    else
      ! Just create a list with the primary mesh id in it
      call init_hierarchical_mesh_id_list( (/ mesh_id /) )
    end if

    call create_geometric_constants(mesh_id_list,      &
                                    chi_list,          &
                                    panel_id_list,     &
                                    label_list         )

    ! Finite element constants should be created after geometric constants
    ! The chi fields set up in geometric constants are used here
    call create_fem_constants(mesh_id_list,      &
                              chi_list,          &
                              panel_id_list,     &
                              label_list, dt     )

    call create_physical_op_constants(mesh_id_list,  &
                                      chi_list,      &
                                      panel_id_list, &
                                      label_list, dt )

    if ( limited_area ) then
      call create_limited_area_constants(mesh_id_list, &
                                         chi_list,     &
                                         label_list    )
    end if

    if ( moisture_conservation ) then
      call create_intermesh_constants(mesh_id,               &
                                      chi,                   &
                                      panel_id,              &
                                      shifted_mesh_id,       &
                                      shifted_chi,           &
                                      double_level_mesh_id,  &
                                      double_level_chi)
    end if

    ! Set-up arrays for transport coefficients
    call runge_kutta_init()
    call flux_alg_set_num_meshes( num_meshes )
    call advective_update_set_num_meshes( num_meshes )
    call rk_transport_rho_set_num_meshes( num_meshes )
    call rk_transport_theta_set_num_meshes( num_meshes )

    deallocate(mesh_id_list)
    deallocate(label_list)

    call log_event( "Gungho: created runtime_constants", LOG_LEVEL_INFO )
    if ( subroutine_timers ) call timer('runtime_constants_alg')

  end subroutine create_runtime_constants


  !> @brief Explicitly reclaim memory from module scope variables
  !
  subroutine final_runtime_constants()

    ! Other runtime_constants modules
    use advective_update_alg_mod,    only: advective_update_alg_final
    use fem_constants_mod,           only: final_fem_constants
    use flux_alg_mod,                only: flux_alg_final
    use geometric_constants_mod,     only: final_geometric_constants
    use intermesh_constants_mod,     only: final_intermesh_constants
    use limited_area_constants_mod,  only: final_limited_area_constants
    use physical_op_constants_mod,   only: final_physical_op_constants
    use rk_transport_rho_mod,        only: rk_transport_rho_final
    use rk_transport_theta_mod,      only: rk_transport_theta_final
    use runge_kutta_init_mod,        only: runge_kutta_final
    use runtime_tools_mod,           only: final_mesh_id_list, &
                                           final_hierarchical_mesh_id_list

    implicit none

    call final_geometric_constants()
    call final_fem_constants()
    call final_physical_op_constants()
    if ( limited_area ) call final_limited_area_constants()
    if ( moisture_conservation ) call final_intermesh_constants()
    call rk_transport_theta_final()
    call rk_transport_rho_final()
    call flux_alg_final()
    call advective_update_alg_final()
    call runge_kutta_final()
    call final_hierarchical_mesh_id_list()
    call final_mesh_id_list()


  end subroutine final_runtime_constants

end module runtime_constants_mod
