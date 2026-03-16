!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Subroutines to for interfacing runtime partitioning using LFRic
!>        infrastructure objects.

module runtime_partition_lfric_mod

  use constants_mod,           only: i_def, l_def
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
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

  interface get_partition_parameters
    procedure get_partition_parameters_cfg
    procedure get_partition_parameters_nml
  end interface get_partition_parameters

contains

subroutine get_partition_parameters_cfg( configuration,  &
                                         mesh_selection, &
                                         total_ranks,    &
                                         decomposition,  &
                                         partitioner_ptr )

  implicit none

  type(namelist_collection_type), intent(in) :: configuration

  integer,        intent(in) :: mesh_selection
  integer(i_def), intent(in) :: total_ranks

  class(panel_decomposition_type), intent(inout), allocatable :: decomposition

  type(namelist_type), pointer :: partitioning

  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  partitioning => configuration%get_namelist('partitioning')

  call get_partition_parameters_nml( partitioning,     &
                                     mesh_selection,   &
                                     total_ranks,      &
                                     decomposition,    &
                                     partitioner_ptr )


end subroutine get_partition_parameters_cfg

subroutine get_partition_parameters_nml( partitioning,   &
                                         mesh_selection, &
                                         total_ranks,    &
                                         decomposition,  &
                                         partitioner_ptr )

  implicit none


  type(namelist_type), intent(in), pointer :: partitioning

  integer,        intent(in) :: mesh_selection
  integer(i_def), intent(in) :: total_ranks

  class(panel_decomposition_type), intent(inout), allocatable :: decomposition

  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  integer(i_def) :: panel_xproc, panel_yproc

  integer :: panel_decomposition

  call partitioning%get_value( 'panel_decomposition', panel_decomposition )

  select case (panel_decomposition)

  case ( panel_decomposition_auto )
    decomposition = auto_decomposition_type()

  case ( panel_decomposition_row )
    decomposition = row_decomposition_type()

  case ( panel_decomposition_column )
    decomposition = column_decomposition_type()

  case ( panel_decomposition_custom )
    call partitioning%get_value( 'panel_xproc', panel_xproc )
    call partitioning%get_value( 'panel_yproc', panel_yproc )
    decomposition = custom_decomposition_type( panel_xproc, panel_yproc )

  case ( panel_decomposition_auto_nonuniform )
    decomposition = auto_nonuniform_decomposition_type()

  case ( panel_decomposition_guided_nonuniform )
    call partitioning%get_value( 'panel_xproc', panel_xproc )
    decomposition = guided_nonuniform_decomposition_type( panel_xproc )

  case default
    ! Not clear it's possible to still error at this point but no harm in checking
    call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select

  call get_partition_strategy(mesh_selection, total_ranks, partitioner_ptr)

end subroutine get_partition_parameters_nml


end module runtime_partition_lfric_mod
