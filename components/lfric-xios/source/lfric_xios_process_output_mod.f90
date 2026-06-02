!-------------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief  Module to only hold routines for re-processing files output by XIOS
!>         where the configuration is unable to be delivered by XIOS.
!>         The only facet in scope for this module is post processing projected
!>         coordinates for UGRID specification files, pending new feature
!>         development in XIOS.  This is technical debt, and this
!>         post processing approach should not be used for other file
!>         manipulations.
!>
module lfric_xios_process_output_mod

  use constants_mod,              only: i_def, str_max_filename
  use file_mod,                   only: FILE_MODE_WRITE,     &
                                        FILE_OP_OPEN
  use io_config_mod,              only: file_convention,       &
                                        file_convention_ugrid
  use lfric_mpi_mod,              only: global_mpi
  use lfric_ncdf_field_mod,       only: lfric_ncdf_field_type
  use lfric_ncdf_file_mod,        only: lfric_ncdf_file_type
  use lfric_xios_constants_mod,   only: dp_xios
  use lfric_xios_file_mod,        only: lfric_xios_file_type
  use log_mod,                    only: log_event, log_level_debug

  implicit none

  public :: process_output_file
  private

  ! Public scaling factor for planar mesh coordinates to circumvent XIOS issue
  real(kind=dp_xios), public, parameter :: xyz_scaling_factor = 1.0e-4_dp_xios

contains

!> @brief Processes a NetCDF file produced by XIOS to work around
!>        limmitations in UGRID projected coordinates only.
!>
!> @param[in] file  The lfric_xios file object for the NetCDF file to be edited
subroutine process_output_file(file)

  implicit none

  character(len=str_max_filename) :: file_path
  class(lfric_xios_file_type), intent(in) :: file
  type(lfric_ncdf_file_type) :: file_ncdf
  logical                    :: file_exists

  ! Output processing must be done in serial
  if (global_mpi%get_comm_rank() /= 0) return

  file_path = file%get_filepath()
  call log_event("Processing output file to format_mesh: "//trim(file_path), &
                  log_level_debug)

  ! If file has not been written out, then don't attempt to process it
  inquire(file=trim(file_path), exist=file_exists)
  if (.not. file_exists) return

  ! Open output file
  file_ncdf = lfric_ncdf_file_type( trim(file_path),           &
                                    open_mode=FILE_OP_OPEN, &
                                    io_mode=FILE_MODE_WRITE )

  call format_mesh(file_ncdf)

  call file_ncdf%close_file()

end subroutine process_output_file

!> @brief Formats the mesh object in the output file
!>
!> @param[in] file_ncdf  The netcdf file to be edited
subroutine format_mesh(file_ncdf)

  implicit none

  type(lfric_ncdf_file_type), intent(inout) :: file_ncdf

  type(lfric_ncdf_field_type) :: mesh_var

  if (.not. file_ncdf%contains_var("Mesh2d")) return

  mesh_var = lfric_ncdf_field_type("Mesh2d", file_ncdf)

  ! This post process code should only be called on a file where
  ! the model has Planar geometry and is UGRID encoded.
  ! This is controlled by lfric_xios_context_mod,
  ! lfric_xios_file_mod and lfric_xios_setup_mod
  call fix_planar_coordinates(file_ncdf)

end subroutine format_mesh


!> @brief Fixes issues with planar coordinates in output file
!>
!> @param[in] file_ncdf  The netcdf file to be edited
subroutine fix_planar_coordinates(file_ncdf)

  implicit none

  type(lfric_ncdf_file_type), intent(inout) :: file_ncdf

  type(lfric_ncdf_field_type) :: coord_field
  character(len=1)            :: dim
  character(len=13)           :: coord_field_names(6) = [ "Mesh2d_node_x", &
                                                          "Mesh2d_node_y", &
                                                          "Mesh2d_edge_x", &
                                                          "Mesh2d_edge_y", &
                                                          "Mesh2d_face_x", &
                                                          "Mesh2d_face_y"  ]
  integer(i_def)              :: i

  ! Modify coordinate fields to correctly represent output from model
  do i = 1, size(coord_field_names)
    ! Trim the last character from the dim name to get the dimension
    dim = coord_field_names(i)(len(coord_field_names(i)):len(coord_field_names(i)))

    coord_field = lfric_ncdf_field_type(trim(coord_field_names(i)), file_ncdf)
    call coord_field%set_char_attribute("standard_name", "projection_"//trim(dim)//"_coordinate")
    call coord_field%set_char_attribute("long_name", trim(dim)//" coordinate of projection")
    call coord_field%set_char_attribute("units","m")
    call coord_field%set_real_attribute("scale_factor", xyz_scaling_factor**(-1.0_dp_xios))
  end do

end subroutine fix_planar_coordinates

end module lfric_xios_process_output_mod
