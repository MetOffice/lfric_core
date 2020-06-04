!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief init functionality for the io_dev miniapp

!> @details Handles creation and initialisation of IO test fields
!>
module io_dev_init_mod

  ! Infrastructure
  use constants_mod,                  only : i_def, str_def
  use field_mod,                      only : field_type, field_proxy_type
  use field_parent_mod,               only : field_parent_type, read_interface, write_interface
  use field_collection_mod,           only : field_collection_type, field_collection_iterator_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W2, Wtheta, W3
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  ! Configuration
  use finite_element_config_mod,      only : element_order
  ! I/O methods
  use read_methods_mod,               only : read_field_face, &
                                             read_field_single_face, &
                                             read_state
  use write_methods_mod,              only : write_field_node, &
                                             write_field_face, &
                                             write_field_single_face

  implicit none

  private
  public :: create_io_dev_fields, &
            init_io_dev_fields

  contains

  !> @details Creates fields used for inputting and outputting IO_Dev data
  !> @param[in]  mesh_id       The identifier given to the current 3d mesh
  !> @param[in]  twod_mesh_id  The identifier given to the current 2d mesh
  !> @param[out] depository    The depository field collection
  !> @param[out] input_fields  Collection of fields to be read
  !> @param[out] output_fields Collection of fields to be written
  subroutine create_io_dev_fields( mesh_id,      &
                                   twod_mesh_id, &
                                   depository,   &
                                   input_fields, &
                                   output_fields )

    implicit none

    ! Arguments
    integer(i_def),              intent(in)  :: mesh_id
    integer(i_def),              intent(in)  :: twod_mesh_id
    type(field_collection_type), intent(out) :: depository
    type(field_collection_type), intent(out) :: input_fields
    type(field_collection_type), intent(out) :: output_fields

    ! Local variables
    type(field_type)       :: input_face_field, input_single_face_field
    type(field_type)       :: node_field, face_field, single_face_field
    type(field_proxy_type) :: tmp_proxy
    integer(i_def)         :: u

    ! Pointers
    class(pure_abstract_field_type), pointer :: tmp_field_ptr => null()
    procedure(read_interface),       pointer :: tmp_read_ptr => null()
    procedure(write_interface),      pointer :: tmp_write_ptr => null()

    call log_event( 'IO_Dev: creating model data', LOG_LEVEL_INFO )

    !----------------------------------------------------------------------------
    ! Create core fields to send/recieve data from file and set I/O behaviours
    !----------------------------------------------------------------------------
    ! Create the depository, input and output field collections.
    depository    = field_collection_type(name='depository')
    input_fields  = field_collection_type(name="input_fields")
    output_fields = field_collection_type(name="output_fields")

    ! W0 (node) field
    call node_field%initialise( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, W0), &
                   name = 'node_field' )
    tmp_write_ptr => write_field_node
    call node_field%set_write_behaviour( tmp_write_ptr )
    call depository%add_field( node_field )

    ! W3 (face) fields
    call face_field%initialise( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, W3), &
                   name = 'face_field' )
    tmp_write_ptr => write_field_face
    call face_field%set_write_behaviour( tmp_write_ptr )
    call depository%add_field( face_field )

    call input_face_field%initialise( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, W3), &
                   name = 'input_face_field' )
    tmp_read_ptr => read_field_face
    call input_face_field%set_read_behaviour( tmp_read_ptr )
    call depository%add_field( input_face_field )

    ! W3_2D (single_face) fields
    call single_face_field%initialise( vector_space = &
                   function_space_collection%get_fs(twod_mesh_id, element_order, W3), &
                   name = 'single_face_field' )
    tmp_write_ptr => write_field_single_face
    call single_face_field%set_write_behaviour( tmp_write_ptr )
    call depository%add_field( single_face_field )

    call input_single_face_field%initialise( vector_space = &
                   function_space_collection%get_fs(twod_mesh_id, element_order, W3), &
                   name = 'input_single_face_field' )
    tmp_read_ptr => read_field_single_face
    call input_single_face_field%set_read_behaviour( tmp_read_ptr )
    call depository%add_field( input_single_face_field )

    !----------------------------------------------------------------------------
    ! Output fields
    !----------------------------------------------------------------------------
    ! Add fields to output_fields collection

    tmp_field_ptr => depository%get_field( 'node_field' )
    call output_fields%add_reference_to_field( tmp_field_ptr )

    tmp_field_ptr => depository%get_field( 'face_field' )
    call output_fields%add_reference_to_field( tmp_field_ptr )

    tmp_field_ptr => depository%get_field( 'single_face_field' )
    call output_fields%add_reference_to_field( tmp_field_ptr )

    nullify( tmp_read_ptr )
    nullify( tmp_write_ptr )
    nullify( tmp_field_ptr )

    call log_event( 'IO_Dev: fields created', LOG_LEVEL_INFO )

  end subroutine create_io_dev_fields

  !> @details Initialises model fields by reading data from file or
  !>          initialising to default value
  !> @param[in,out] input_fields  Collection of fields to be read
  !> @param[in,out] output_fields Collection of fields to be written
  subroutine init_io_dev_fields( input_fields, &
                                 output_fields )

    implicit none

    ! Arguments
    type(field_collection_type), intent(inout) :: input_fields
    type(field_collection_type), intent(inout) :: output_fields

    ! Local variables
    integer(i_def) :: dof_index
    character(str_def) :: input_field_name
    type(field_collection_iterator_type) :: iter
    type(field_proxy_type) :: input_proxy, core_proxy

    ! Pointers
    class(field_parent_type), pointer :: fld => null()
    class(field_type),        pointer :: input_fld_ptr => null()


    ! Read input fields from file
    ! This line will be enabled with the second round of rose-stem tests, which
    ! will include the reading of data generated in the first round.
    !call read_state( input_fields )

    ! Populate core fields from input fields
    iter = output_fields%get_iterator()
    do
      if ( .not.iter%has_next() ) exit
      fld => iter%next()

      select type(fld)
        type is (field_type)
          if ( input_fields%field_exists( 'input_'//trim(adjustl(fld%get_name())) ) ) then

            ! Get field from input collection
            input_fld_ptr => input_fields%get_field( 'input_'//trim(adjustl(fld%get_name())) )
            ! Get field proxies for input and core field
            input_proxy = input_fld_ptr%get_proxy()
            core_proxy = fld%get_proxy()
            ! Copy input field data to field in core collection
            core_proxy%data = input_proxy%data

          else
            ! Initialise field to default value (number of dofs)
            core_proxy = fld%get_proxy()
            do dof_index = 1, core_proxy%vspace%get_last_dof_owned()
              core_proxy%data(dof_index) = core_proxy%vspace%get_last_dof_owned()
            end do

          end if

      end select

    end do

    nullify(fld)
    nullify(input_fld_ptr)

  end subroutine init_io_dev_fields

end module io_dev_init_mod
