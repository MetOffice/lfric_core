!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_ancils_mod

  use constants_mod,                  only : i_def, l_def
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO
  use init_tstar_analytic_alg_mod,    only : init_tstar_analytic_alg

  use field_mod,                      only : field_type,     &
                                             read_interface, &
                                             write_interface
  use io_config_mod,                  only : use_xios_io
  use io_mod,                         only : xios_read_field_face, &
                                             xios_read_field_single_face, &
                                             xios_write_field_face, &
                                             xios_write_field_single_face
  use field_collection_mod,           only : field_collection_type, &
                                             field_collection_iterator_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3

  implicit none

  private  :: setup_ancil_field
  public   :: init_analytic_ancils,     &
              init_aquaplanet_ancils,   &
              create_fd_ancils,         &
              test_ancil_output

contains


  !> @details Initialises ancillary fields analytically
  !> @param[in,out] surface_fields the 2D field collection
  subroutine init_analytic_ancils(surface_fields)

    implicit none

    type( field_collection_type ), intent( inout ) :: surface_fields

    type( field_type ), pointer ::    tstar_ptr  => null()

    tstar_ptr => surface_fields%get_field('tstar')

    call init_tstar_analytic_alg(tstar_ptr)

    nullify(tstar_ptr)

  end subroutine init_analytic_ancils

  !> @details Initialises ancillary fields from
  !>          an aquaplanet dump
  !> @param[in,out] surface_fields the 2D field collection
  subroutine init_aquaplanet_ancils(surface_fields)

    implicit none

    type( field_collection_type ), intent( inout ) :: surface_fields

    ! local variables
    ! Pointer to the 2D tstar in the surface field collection
    type( field_type ), pointer ::    tstar_ptr  => null()

    procedure(read_interface), pointer  :: tmp_read_ptr => null()

    call log_event("Reading tstar from dump", LOG_LEVEL_INFO)

    tstar_ptr => surface_fields%get_field('tstar')


    ! Need to set the I/O handler for read. Any ancils here
    ! are currently read from a UM2LFRic dump

    tmp_read_ptr => xios_read_field_single_face
    call tstar_ptr%set_read_behaviour(tmp_read_ptr)
    call tstar_ptr%read_field("read_tstar")

    nullify( tstar_ptr, tmp_read_ptr )

  end subroutine init_aquaplanet_ancils

  !> @details Organises fields to be read from ancils into ancil_fields
  !           collection then reads them.
  !> @param[in,out] depository The depository field collection
  !> @param[out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh_id The identifier given to the current 3d mesh
  !> @param[in] twod_mesh_id The identifier given to the current 2d mesh
  subroutine create_fd_ancils(depository, ancil_fields, mesh_id, twod_mesh_id)

    implicit none

    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( out )   :: ancil_fields
    integer(i_def), intent(in) :: mesh_id
    integer(i_def), intent(in) :: twod_mesh_id

    ! Set up ancil_fields collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Setting up ancil field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    ancil_fields = field_collection_type(name='ancil_fields')

    ! For fields that are not already created, we create them here then call the subroutine
    ! to set up the read behaviour and field collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Creating test - land area fraction & soil surface fields."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    !=====  SURFACE ANCILS  =====
    call setup_ancil_field("land_area_fraction", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)

    !=====  SOIL ANCILS  =====
    call setup_ancil_field("soil_albedo", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("soil_carbon_content", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("soil_thermal_cond", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)

    !=====  OROGRAPHY ANCILS  =====
    call setup_ancil_field("sd_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("grad_xx_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("grad_xy_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("grad_yy_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("peak_to_trough_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)
    call setup_ancil_field("silhouette_area_orog", depository, ancil_fields, mesh_id, twod_mesh_id, twod=.true.)

    ! Now the field collection is set up, we will read all the fields in the collection
    ! in gungho_model_data_mod (call to init_fd_ancils)

  end subroutine create_fd_ancils

  !> @details Creates fields to be read into from ancillary files -
  !           only used if field is not already created elsewhere
  !> @param[in] name The field name
  !> @param[in,out] depository The depository field collection
  !> @param[in,out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh_id The identifier given to the current 3d mesh
  !> @param[in] twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in] twod Flag to determine if field is on 2D mesh or regular
  subroutine setup_ancil_field(name, depository, ancil_fields, mesh_id, twod_mesh_id, twod)

    implicit none

    character(*), intent(in)                       :: name
    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( inout ) :: ancil_fields
    integer(i_def), intent(in)                     :: mesh_id
    integer(i_def), intent(in)                     :: twod_mesh_id
    logical(l_def), optional, intent(in)           :: twod

    !Local variables
    type(field_type)           :: new_field
    integer(i_def), parameter  :: fs_order = 0

    ! Pointers
    type(function_space_type), pointer  :: w3_space => null()
    type(function_space_type), pointer  :: twod_space => null()
    procedure(read_interface), pointer  :: tmp_read_ptr => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()
    type(field_type), pointer           :: tgt_ptr => null()

    w3_space => function_space_collection%get_fs(mesh_id, fs_order, W3)
    twod_space => function_space_collection%get_fs(twod_mesh_id, fs_order, W3)

    ! If field does not yet exist, then create it
    if (.not. depository%field_exists( name ) ) then
      if (present(twod)) then
        call new_field%initialise( twod_space, name=trim(name) )
      else
        call new_field%initialise( w3_space, name=trim(name) )
      end if
      ! Add the new field to the field depository
      call depository%add_field(new_field)
    end if

    ! Get a field pointer from the depository
    tgt_ptr => depository%get_field(name)

    !Set up field read behaviour for 2D and 3D fields
    if (present(twod))then
      tmp_read_ptr => xios_read_field_single_face
      tmp_write_ptr => xios_write_field_single_face
    else
      tmp_read_ptr => xios_read_field_face
      tmp_write_ptr => xios_write_field_face
    end if

    ! Set field read behaviour for target field
    call tgt_ptr%set_read_behaviour(tmp_read_ptr)
    call tgt_ptr%set_write_behaviour(tmp_write_ptr)

    ! Add the field pointer to the ancil_fields collection
    call ancil_fields%add_reference_to_field( tgt_ptr )

    ! Nullify pointers
    nullify(w3_space)
    nullify(twod_space)
    nullify(tmp_read_ptr)
    nullify(tmp_write_ptr)
    nullify(tgt_ptr)

  end subroutine setup_ancil_field

  !> @details Writes ancils to file so that they can be checked
  !> @param[in] fld_collection The field collection used for ancils
  subroutine test_ancil_output(fld_collection)

    implicit none

    type( field_collection_type ), intent(inout) :: fld_collection

    type( field_collection_iterator_type) :: iter
    type( field_type ), pointer :: fld => null()

    iter=fld_collection%get_iterator()
    do
      if(.not.iter%has_next())exit
      fld=>iter%next()
      if (fld%can_write()) then
        write(log_scratch_space,'(3A,I6)') &
            "Writing ", trim(adjustl(fld%get_name()))//'_out'
        call log_event(log_scratch_space,LOG_LEVEL_INFO)
        call fld%write_field( trim(adjustl(fld%get_name()))//'_out' )
      else
        call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                        ' not set up', LOG_LEVEL_INFO )
      end if

    end do

    nullify(fld)

  end subroutine test_ancil_output

end module init_ancils_mod
