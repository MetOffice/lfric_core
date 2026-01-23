!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Subroutines to for interfacing runtime partitioning using LFRic
!>        infrastructure objects.

module runtime_partition_lfric_mod

  use constants_mod,           only: i_def, l_def
  use config_mod,              only: config_type
  use partitioning_nml_mod,    only: partitioning_nml_type
  use partition_mod,           only: partitioner_interface
  use runtime_partition_mod,   only: get_partition_strategy
  use panel_decomposition_mod, only: panel_decomposition_type,           &
                                     auto_decomposition_type,            &
                                     row_decomposition_type,             &
                                     column_decomposition_type,          &
                                     custom_decomposition_type,          &
                                     auto_nonuniform_decomposition_type, &
                                     guided_nonuniform_decomposition_type
  use partitioning_config_mod, only: panel_decomposition_auto,            &
                                     panel_decomposition_row,             &
                                     panel_decomposition_column,          &
                                     panel_decomposition_custom,          &
                                     panel_decomposition_auto_nonuniform, &
                                     panel_decomposition_guided_nonuniform

  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public :: get_partition_parameters


contains

subroutine get_partition_parameters( partitioning_nml, &
                                     mesh_selection,   &
                                     total_ranks,      &
                                     decomposition,    &
                                     partitioner_ptr )

  implicit none


  type(partitioning_nml_type), intent(in) :: partitioning_nml

  integer,        intent(in) :: mesh_selection
  integer(i_def), intent(in) :: total_ranks

  class(panel_decomposition_type), intent(inout), allocatable :: decomposition

  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  integer(i_def) :: panel_xproc, panel_yproc

  integer(i_def) :: panel_decomposition

  panel_decomposition = partitioning_nml%panel_decomposition()
  panel_xproc = partitioning_nml%panel_xproc()
  panel_yproc = partitioning_nml%panel_yproc()

  select case (panel_decomposition)

  case ( panel_decomposition_auto )
    decomposition = auto_decomposition_type()

  case ( panel_decomposition_row )
    decomposition = row_decomposition_type()

  case ( panel_decomposition_column )
    decomposition = column_decomposition_type()

  case ( panel_decomposition_custom )
    decomposition = custom_decomposition_type( panel_xproc, panel_yproc )

  case ( panel_decomposition_auto_nonuniform )
    decomposition = auto_nonuniform_decomposition_type()

  case ( panel_decomposition_guided_nonuniform )
    decomposition = guided_nonuniform_decomposition_type( panel_xproc )

  case default
    ! Not clear it's possible to still error at this point but no harm in checking
    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select

  call get_partition_strategy(mesh_selection, total_ranks, partitioner_ptr)

end subroutine get_partition_parameters


end module runtime_partition_lfric_mod
