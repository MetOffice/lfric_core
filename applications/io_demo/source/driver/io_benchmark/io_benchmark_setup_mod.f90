!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Setup infrastructure used for I/O benchmark
!> @details Handles the setup of all the fields that will be passed used to
!!          benchmark the speed of XIOS reading and writing
module io_benchmark_setup_mod

  use constants_mod,                 only: i_def, str_def
  use driver_modeldb_mod,            only: modeldb_type
  use field_collection_mod,          only: field_collection_type
  use field_mod,                     only: field_type
  use field_parent_mod,              only: read_interface, write_interface
  use file_mod,                      only: FILE_MODE_WRITE
  use fs_continuity_mod,             only: Wtheta
  use function_space_mod,            only: function_space_type
  use function_space_collection_mod, only: function_space_collection
  use lfric_xios_file_mod,           only: lfric_xios_file_type, OPERATION_TIMESERIES
  use lfric_xios_read_mod,           only: read_field_generic
  use lfric_xios_write_mod,          only: write_field_generic
  use linked_list_mod,               only: linked_list_type
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection

  implicit none

  public create_io_benchmark_fields, setup_io_benchmark_files

contains

  !> @details Creates the fields needed for the IO benchmark
  !> @param[in,out] modeldb The model database in which to store model data.
  subroutine create_io_benchmark_fields(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(mesh_type), pointer :: mesh
    type(field_collection_type), pointer :: io_benchmark_fields
    type(field_type)                     :: tmp_io_field
    procedure(read_interface),  pointer  :: tmp_read_ptr
    procedure(write_interface), pointer  :: tmp_write_ptr
    type(function_space_type),  pointer  :: wtheta_fs

    character(str_def) :: prime_mesh_name, tmp_field_name
    integer(i_def) :: element_order_h
    integer(i_def) :: element_order_v
    integer(i_def) :: i
    integer(i_def) :: n_benchmark_fields
    integer(i_def) :: diagnostic_frequency

    prime_mesh_name  = modeldb%config%base_mesh%prime_mesh_name()
    element_order_h  = modeldb%config%finite_element%element_order_h()
    element_order_v  = modeldb%config%finite_element%element_order_v()
    n_benchmark_fields = modeldb%config%io_demo%n_benchmark_fields()
    diagnostic_frequency = modeldb%config%io%diagnostic_frequency()

    mesh => mesh_collection%get_mesh(prime_mesh_name)

    call modeldb%fields%add_empty_field_collection("io_benchmark_fields")
    io_benchmark_fields => modeldb%fields%get_field_collection("io_benchmark_fields")
    wtheta_fs => function_space_collection%get_fs( mesh, element_order_h, &
                                                   element_order_v, Wtheta )

    do i = 1, n_benchmark_fields
      write(tmp_field_name, "(A19, I3.3)") 'io_benchmark_field_', i
      call tmp_io_field%initialise( vector_space = wtheta_fs, &
                                    name=tmp_field_name )
      tmp_read_ptr => read_field_generic
      tmp_write_ptr => write_field_generic
      call tmp_io_field%set_read_behaviour(tmp_read_ptr)
      call tmp_io_field%set_write_behaviour(tmp_write_ptr)
      call io_benchmark_fields%add_field(tmp_io_field)
    end do

    nullify( mesh, io_benchmark_fields, wtheta_fs )

  end subroutine create_io_benchmark_fields

  subroutine setup_io_benchmark_files(file_list, modeldb)

    implicit none

    type(linked_list_type),       intent(out)   :: file_list
    type(modeldb_type), optional, intent(inout) :: modeldb

    integer(i_def)                       :: diagnostic_frequency
    type(field_collection_type), pointer :: io_benchmark_fields

    diagnostic_frequency = modeldb%config%io%diagnostic_frequency()

    io_benchmark_fields => modeldb%fields%get_field_collection("io_benchmark_fields")

    file_list = linked_list_type()
    call file_list%insert_item( lfric_xios_file_type( "lfric_xios_write_benchmark",         &
                                                      xios_id="lfric_xios_write_benchmark", &
                                                      io_mode=FILE_MODE_WRITE,              &
                                                      operation=OPERATION_TIMESERIES,       &
                                                      freq=diagnostic_frequency,            &
                                                      fields_in_file=io_benchmark_fields ) )

  nullify(io_benchmark_fields)

  end subroutine setup_io_benchmark_files

end module io_benchmark_setup_mod