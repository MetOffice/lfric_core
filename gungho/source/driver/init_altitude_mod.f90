!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Create and initialise surface altitude fields
!> @details Creates the surface altitude fields and populates them either with
!>          data read from an ancillary file or data calculated analytically
module init_altitude_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type,  &
                                             read_interface,  &
                                             write_interface
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use orography_config_mod,           only : orog_init_option,  &
                                             orog_init_option_analytic, &
                                             orog_init_option_ancil
  use orography_control_mod,          only : set_orography_option
  use io_config_mod,                  only : use_xios_io
  use io_mod,                         only : xios_read_field_single_face, &
                                             xios_write_field_single_face
  use log_mod,                        only : log_event, &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,  &
                                             LOG_LEVEL_ERROR
  use xios,                           only : xios_is_valid_file

  implicit none

  private

  public :: init_altitude

  contains
  !> @brief Routine to initialise surface altitude field used for orography
  !>        - either from ancil file or by calculating it analytically.
  !> @param[in]   twod_mesh_id The identifier given to the current 2d mesh
  !> @param[out]  surface_altitude Surface altitude field
  subroutine init_altitude( twod_mesh_id, surface_altitude )

    implicit none

    integer(i_def), intent(in)      :: twod_mesh_id
    type(field_type), intent(out)   :: surface_altitude

    ! Pointers to vector space and read interface
    type(function_space_type), pointer  :: vector_space => null()
    procedure(read_interface), pointer  :: read_ptr => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Set the function space order to be 0 - ancil files currently store data
    ! on 2D face centres
    integer(i_def), parameter :: element_order = 0

    write(log_scratch_space,'(A,A)') "Initialise surface altitude: "// &
          "Create surface altitude."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Here we create the 2D surface altitude field
    vector_space => function_space_collection%get_fs( twod_mesh_id, &
                                                      element_order, W3)

    call surface_altitude%initialise(vector_space, name=trim('surface_altitude'))

    select case ( orog_init_option )

      ! Analytic orography
      case ( orog_init_option_analytic )

        ! Set up analytic orography parameters
        call set_orography_option()
        ! In the future we will initialise the surface altitude field based on the
        ! orography configuration switch.

      ! Orography read from ancil file
      case ( orog_init_option_ancil )

        ! Set up the read behaviour for our field and read it
        if (use_xios_io) then

          ! Altitude field should read in on W3_2D function space
          read_ptr => xios_read_field_single_face
          call surface_altitude%set_read_behaviour(read_ptr)

          write(log_scratch_space,'(A,A)') "Initialise surface altitude: "// &
                "Read surface altitude from ancil."
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
          call surface_altitude%read_field(trim(adjustl(surface_altitude%get_name())))

          ! Write behaviour also set up for diagnostic output of altitude if required
          tmp_write_ptr => xios_write_field_single_face
          call surface_altitude%set_write_behaviour(tmp_write_ptr)

        ! Call error if XIOS is not being used - ancils cannot be read without it
        else
          write(log_scratch_space,'(A,A)') "Initialise surface altitude: "// &
                "XIOS must be enabled to read surface altitude from ancil."
          call log_event(log_scratch_space, LOG_LEVEL_INFO)

        end if

      ! Default - no orography
      case default

        write(log_scratch_space,'(A,A)') "Initialise surface altitude: "// &
              "No initialisation option set for orography."
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

    end select

    nullify(vector_space)
    nullify(read_ptr)
    nullify(tmp_write_ptr)

  end subroutine init_altitude

end module init_altitude_mod
