!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @todo Remove create_io_demo_constants if/when #2389 moves dx_at_w2 related
!>       kernels outside of gungho.

!> @brief Initialisation functionality for the io_demo miniapp
!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators
module init_io_demo_mod

  use sci_assign_field_random_range_alg_mod,  only: assign_field_random_range
  use constants_mod,                          only : i_def, r_def, l_def
  use driver_modeldb_mod,                     only : modeldb_type
  use field_collection_mod,                   only : field_collection_type
  use field_mod,                              only : field_type
  use field_parent_mod,                       only : write_interface
  use function_space_collection_mod,          only : function_space_collection
  use function_space_mod,                     only : function_space_type
  use fs_continuity_mod,                      only : Wtheta
  use key_value_mod,                          only : abstract_value_type
  use log_mod,                                only : log_event,       &
                                                     LOG_LEVEL_TRACE, &
                                                     LOG_LEVEL_ERROR
  use mesh_mod,                               only : mesh_type
  use lfric_xios_write_mod,                   only : write_field_generic
  use io_demo_constants_mod,                  only : create_io_demo_constants
  use random_number_generator_mod,            only : random_number_generator_type

  implicit none

  contains

  !> @details Initialises everything needed to run the io_demo miniapp
  !> @param[in,out] modeldb The structure that holds model state
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in,out] chi The co-ordinate field
  !> @param[in,out] panel_id 2d field giving the id for cubed sphere panels
  subroutine init_io_demo(modeldb, mesh, chi, panel_id)

    implicit none

    type(modeldb_type), intent(inout)       :: modeldb
    type(mesh_type),    intent(in), pointer :: mesh

    ! Coordinate field
    type(field_type), intent(inout) :: chi(:)
    type(field_type), intent(inout) :: panel_id

    class(abstract_value_type), pointer         :: abstract_value
    type(random_number_generator_type), pointer :: rng
    type(field_type)                            :: diffusion_field
    type(field_collection_type), pointer        :: depository
    procedure(write_interface), pointer         :: tmp_ptr
    real(kind=r_def), parameter                 :: min_val = 280.0_r_def
    real(kind=r_def), parameter                 :: max_val = 330.0_r_def

    type(function_space_type), pointer :: fs

    integer(i_def) :: order_h, order_v
    logical(l_def) :: write_diag
    logical(l_def) :: use_xios_io

    call log_event( 'io_demo: Initialising miniapp ...', LOG_LEVEL_TRACE )

    order_h = modeldb%config%finite_element%element_order_h()
    order_v = modeldb%config%finite_element%element_order_v()

    write_diag  = modeldb%config%io%write_diag()
    use_xios_io = modeldb%config%io%use_xios_io()

    ! seed the random number generator
    call modeldb%values%get_value("rng", abstract_value)
    select type(abstract_value)
      type is (random_number_generator_type)
        rng => abstract_value
      class default
        call log_event( &
          "Error: the value called 'rng' is not a random_number_generator", &
          LOG_LEVEL_ERROR &
        )
    end select
    call rng%check_seed()

    ! Create prognostic fields
    ! Creates a field in the Wtheta function space
    fs => function_space_collection%get_fs(mesh, order_h, order_v, Wtheta)
    call diffusion_field%initialise(fs, name="diffusion_field")

    ! Set up field with an IO behaviour (XIOS only at present)
    if (write_diag .and. use_xios_io) then
       tmp_ptr => write_field_generic
       call diffusion_field%set_write_behaviour(tmp_ptr)
    end if

    ! Add field to modeldb
    depository => modeldb%fields%get_field_collection("depository")

    ! Initialising field
    call assign_field_random_range( diffusion_field, min_val, max_val )

    call depository%add_field(diffusion_field)

    ! Create io_demo runtime constants. This creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_io_demo_constants(modeldb, mesh, chi, panel_id)

    call log_event( 'io_demo: Miniapp initialised', LOG_LEVEL_TRACE )

  end subroutine init_io_demo

end module init_io_demo_mod
