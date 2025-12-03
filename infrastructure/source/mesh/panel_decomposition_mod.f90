!-----------------------------------------------------------------------------
! Copyright (c) 2025,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Provides strategies for decomposing rectangular panels into partitions
!>
module panel_decomposition_mod

  use global_mesh_mod, only: global_mesh_type
  use global_mesh_collection_mod, only: global_mesh_collection_type
  use constants_mod, only: i_def, l_def, r_def
  use log_mod, only: log_event, log_scratch_space, &
                     LOG_LEVEL_ERROR, LOG_LEVEL_INFO, LOG_LEVEL_DEBUG

  implicit none

  !> @brief Parent for panel decomposition types
  type, public, abstract :: panel_decomposition_type
  contains
    procedure(get_partition_interface), deferred :: get_partition
  end type panel_decomposition_type

  !> @brief Decomposition that accepts user specified number of xprocs and yprocs
  type, extends(panel_decomposition_type), public :: custom_decomposition_type
    integer(i_def) :: num_xprocs, num_yprocs
  contains
    procedure, public :: get_partition => get_custom_partition
  end type custom_decomposition_type
  ! Constructor
  interface custom_decomposition_type
    module procedure custom_decomposition_constructor
  end interface

  !> @brief Decomposition that automatically determines number of xprocs and yprocs
  type, extends(panel_decomposition_type), public :: auto_decomposition_type
  contains
    procedure, public :: get_partition => get_auto_partition
  end type auto_decomposition_type

  !> @brief Decomposition only in x direction
  type, extends(panel_decomposition_type), public :: row_decomposition_type
  contains
    procedure, public :: get_partition => get_row_partition
  end type row_decomposition_type

  !> @brief Decomposition only in y direction
  type, extends(panel_decomposition_type), public :: column_decomposition_type
  contains
    procedure, public :: get_partition => get_column_partition
  end type column_decomposition_type

  !> @brief Decomposition that automatically generates a nonuniform decomposition
  type, extends(panel_decomposition_type), public :: auto_nonuniform_decomposition_type
  contains
    procedure, public :: get_partition => get_auto_nonuniform_partition
  end type auto_nonuniform_decomposition_type

  ! @brief Decomposition that accepts user specified number of xprocs to
  !        generate a nonuniform partition
  type, extends(panel_decomposition_type), public :: guided_nonuniform_decomposition_type
    integer(i_def) :: num_xprocs
  contains
    procedure, public :: get_partition => get_guided_nonuniform_partition
  end type guided_nonuniform_decomposition_type
  ! Constructor
  interface guided_nonuniform_decomposition_type
    module procedure guided_nonuniform_decomposition_constructor
  end interface


  ! Interface for routines that generate partition shape and location
  abstract interface

    subroutine get_partition_interface( self,             &
                                        relative_rank,    &
                                        panel_ranks,      &
                                        mapping_factor,   &
                                        num_cells_x,      &
                                        num_cells_y,      &
                                        any_maps,         &
                                        partition_width,  &
                                        partition_height, &
                                        partition_x_pos,  &
                                        partition_y_pos )
      use constants_mod, only: i_def
      import :: panel_decomposition_type

      class(panel_decomposition_type), intent(in) :: self

      integer(i_def), intent(in)    :: relative_rank,    &
                                       panel_ranks,      &
                                       mapping_factor,   &
                                       num_cells_x,      &
                                       num_cells_y
      logical,        intent(in)    :: any_maps

      integer(i_def), intent(inout) :: partition_width,  &
                                       partition_height, &
                                       partition_x_pos,  &
                                       partition_y_pos

    end subroutine get_partition_interface

  end interface

contains

  !> @brief Partition the panel into a given number of x and y processes
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_custom_partition( self,             &
                                   relative_rank,    &
                                   panel_ranks,      &
                                   mapping_factor,   &
                                   num_cells_x,      &
                                   num_cells_y,      &
                                   any_maps,         &
                                   partition_width,  &
                                   partition_height, &
                                   partition_x_pos,  &
                                   partition_y_pos )
    implicit none

    class(custom_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    call log_event("Using custom decomposition", LOG_LEVEL_INFO)

    num_xprocs = self%num_xprocs
    num_yprocs = self%num_yprocs

    if ( panel_ranks /= num_xprocs * num_yprocs ) then
      write( log_scratch_space, "(a,i0,a,i0,a,i0)" ) "Total ranks per panel ", panel_ranks, " must be the product of xprocs ", self%num_xprocs, " and yprocs ", self%num_yprocs
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if


    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_custom_partition


  !> @brief Constructor for custom_decomposition_type
  !> @param[in] xprocs The requested number of partitions in the x direction
  !> @param[in] yprocs The requested number of partitions in the y direction
  function custom_decomposition_constructor(xprocs, yprocs) result(self)
    implicit none

    type(custom_decomposition_type), target :: self
    integer(i_def), intent(in) :: xprocs, yprocs

    self%num_xprocs = xprocs
    self%num_yprocs = yprocs

  end function custom_decomposition_constructor


  !> @brief Partition the panel into an automatically determined number of x and
  !         y processes
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_auto_partition( self,             &
                                 relative_rank,    &
                                 panel_ranks,      &
                                 mapping_factor,   &
                                 num_cells_x,      &
                                 num_cells_y,      &
                                 any_maps,         &
                                 partition_width,  &
                                 partition_height, &
                                 partition_x_pos,  &
                                 partition_y_pos )
    implicit none

    class(auto_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs
    integer(i_def) :: mp_num_cells_x, mp_num_cells_y
    integer(i_def) :: start_xprocs, start_width, i
    logical :: found_partition

    call log_event("Using auto decomposition", LOG_LEVEL_INFO)

    ! For automatic partitioning, try to partition into the squarest
    ! possible partitions.

    mp_num_cells_x = num_cells_x / mapping_factor
    mp_num_cells_y = num_cells_y / mapping_factor

    ! Find width of squarest possible partitions and corresponding xprocs
    start_width = nint( sqrt( &
                    real(mp_num_cells_x * mp_num_cells_y / panel_ranks, kind=r_def) &
                  ), kind=i_def)
    start_xprocs = mp_num_cells_x / start_width

    found_partition = .false.

    do i = 0, start_xprocs - 1
      ! Check higher values
      num_xprocs = start_xprocs + i

      ! num_xprocs must divide panel_ranks
      if ( mod(panel_ranks, num_xprocs) == 0 ) then
        num_yprocs = panel_ranks / num_xprocs

        ! If we have any maps then x and y procs must divide the coarsest panel
        if ( (.not. any_maps) .or. &
             ( mod( mp_num_cells_x, num_xprocs ) == 0 .and. &
               mod( mp_num_cells_y, num_yprocs ) == 0 ) &
        ) then
          found_partition = .true.
          exit
        end if

      end if

      ! Check lower values
      num_xprocs = start_xprocs - i

      ! num_xprocs must divide panel_ranks
      if ( mod(panel_ranks, num_xprocs) == 0 ) then
        num_yprocs = panel_ranks / num_xprocs

        ! If we have any maps then x and y procs must divide the coarsest panel
        if ( (.not. any_maps) .or. &
             ( mod( mp_num_cells_x, num_xprocs ) == 0 .and. &
               mod( mp_num_cells_y, num_yprocs ) == 0 ) &
        ) then
          found_partition = .true.
          exit
        end if

      end if

    end do

    if (.not. found_partition) call log_event( &
      "Could not automatically partition domain.", LOG_LEVEL_ERROR )

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_auto_partition


  !> @brief Partition the panel only in the x direction
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_row_partition( self,             &
                                relative_rank,    &
                                panel_ranks,      &
                                mapping_factor,   &
                                num_cells_x,      &
                                num_cells_y,      &
                                any_maps,         &
                                partition_width,  &
                                partition_height, &
                                partition_x_pos,  &
                                partition_y_pos )
    implicit none

    class(row_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    call log_event("Using row decomposition", LOG_LEVEL_INFO)

    num_xprocs = panel_ranks
    num_yprocs = 1_i_def

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_row_partition


  !> @brief Partition the panel only in the y direction
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_column_partition( self,             &
                                   relative_rank,    &
                                   panel_ranks,      &
                                   mapping_factor,   &
                                   num_cells_x,      &
                                   num_cells_y,      &
                                   any_maps,         &
                                   partition_width,  &
                                   partition_height, &
                                   partition_x_pos,  &
                                   partition_y_pos )
    implicit none

    class(column_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    call log_event("Using column decomposiiton", LOG_LEVEL_INFO)

    num_xprocs = 1_i_def
    num_yprocs = panel_ranks

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_column_partition

  !> @brief Partition the panel into an automatically determined number of
  !         columns of partitions of variable size.
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_auto_nonuniform_partition( self,             &
                                            relative_rank,    &
                                            panel_ranks,      &
                                            mapping_factor,   &
                                            num_cells_x,      &
                                            num_cells_y,      &
                                            any_maps,         &
                                            partition_width,  &
                                            partition_height, &
                                            partition_x_pos,  &
                                            partition_y_pos )
    implicit none

    class(auto_nonuniform_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, mp_num_cells_x, mp_num_cells_y
    integer(i_def) :: start_xprocs, start_width, i
    logical ::found_factors

    call log_event("Using auto_nonuniform decomposition", LOG_LEVEL_INFO)

    mp_num_cells_x = num_cells_x / mapping_factor
    mp_num_cells_y = num_cells_y / mapping_factor

    ! For automatic partitioning, try to partition into the squarest
    ! possible partitions

    ! Find width of squarest possible partitions and corresponding xprocs
    start_width = nint( sqrt( &
                    real(mp_num_cells_x * mp_num_cells_y / panel_ranks, kind=r_def) &
                  ), kind=i_def)
    start_xprocs = mp_num_cells_x / start_width
    found_factors = .false.

    do i = 0, start_xprocs - 1
      num_xprocs = start_xprocs - i

      ! If there are any intermesh maps then xprocs must also divide the domain.
      if ( ( mod(mp_num_cells_x, num_xprocs) == 0 ) .or. .not. any_maps ) then
        found_factors = .true.
        exit
      end if

      num_xprocs = start_xprocs + i

      ! If there are any intermesh maps then xprocs must also divide the domain.
      if ( ( mod(mp_num_cells_x, num_xprocs) == 0 ) .or. .not. any_maps ) then
        found_factors = .true.
        exit
      end if

    end do

    if (.not. found_factors) call log_event( &
      "Could not automatically partition domain.", LOG_LEVEL_ERROR )

    call nonuniform_decomposition(relative_rank,    &
                                  panel_ranks,      &
                                  num_cells_x,      &
                                  num_cells_y,      &
                                  num_xprocs,       &
                                  mapping_factor,   &
                                  partition_width,  &
                                  partition_height, &
                                  partition_x_pos,  &
                                  partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_auto_nonuniform_partition

  !> @brief Partition the panel into a given number of columns of partitions of
  !         variable size.
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    mapping_factor   The ratio between this and coarsest mesh
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_guided_nonuniform_partition( self,             &
                                              relative_rank,    &
                                              panel_ranks,      &
                                              mapping_factor,   &
                                              num_cells_x,      &
                                              num_cells_y,      &
                                              any_maps,         &
                                              partition_width,  &
                                              partition_height, &
                                              partition_x_pos,  &
                                              partition_y_pos )
    implicit none

    class(guided_nonuniform_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     mapping_factor,   &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs

    call log_event("Using guided_nonuniform decomposition", LOG_LEVEL_INFO)

    num_xprocs = self%num_xprocs

    ! Defensive checks
    if ( num_xprocs <= 0 ) then
      call log_event("Number of x processes must be strictly positive.", LOG_LEVEL_ERROR)
    end if

    if ( num_cells_x < num_xprocs ) then
      write(log_scratch_space, "(a,i0,a,i0)") "Must have more cells than partitions in x direction."
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    if ( any_maps .and. ( mod(num_cells_x, num_xprocs) /= 0 ) ) then
      write(log_scratch_space, "(a,i0,a,i0)") "Requested number of ranks in x direction ", num_xprocs, &
        " must divide panel x dimension ", num_cells_x
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    call nonuniform_decomposition(relative_rank,    &
                                  panel_ranks,      &
                                  num_cells_x,      &
                                  num_cells_y,      &
                                  num_xprocs,       &
                                  mapping_factor,   &
                                  partition_width,  &
                                  partition_height, &
                                  partition_x_pos,  &
                                  partition_y_pos )

    write(log_scratch_space, "(a,i0,a,i0,a,i0,a,i0)") &
      " partition_x_pos ",  partition_x_pos, &
      " partition_y_pos ",  partition_y_pos, &
      " partition_width ",  partition_width, &
      " partition_height ", partition_height
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  end subroutine get_guided_nonuniform_partition

  !> @brief Constructor for guided_nonuniform_decomposition_type
  !> @param[in] xprocs The requested number of partitions in the x direction
  function guided_nonuniform_decomposition_constructor( xprocs ) result(self)
    implicit none

    type(guided_nonuniform_decomposition_type), target :: self
    integer(i_def), intent(in) :: xprocs

    self%num_xprocs = xprocs

  end function guided_nonuniform_decomposition_constructor


  !> @brief Helper function for generating identical partitions arranged in a
  !         rectangular grid
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    num_xprocs       The number of partitions in the x direction
  !> @param[in]    num_yprocs       The number of partitions in the y direction
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine xy_decomposition( relative_rank,    &
                               num_cells_x,      &
                               num_cells_y,      &
                               num_xprocs,       &
                               num_yprocs,       &
                               partition_width,  &
                               partition_height, &
                               partition_x_pos,  &
                               partition_y_pos )
    implicit none

    integer(i_def), intent(in)    :: relative_rank,    &
                                     num_cells_x,      &
                                     num_cells_y,      &
                                     num_xprocs,       &
                                     num_yprocs
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: xproc, yproc

    xproc = mod(relative_rank - 1, num_xprocs)
    yproc = (relative_rank - 1) / num_xprocs

    partition_x_pos  = ( (xproc * num_cells_x) / num_xprocs ) + 1
    partition_width  = ( ((xproc+1)*num_cells_x) / num_xprocs ) - partition_x_pos + 1
    partition_y_pos  = ( (yproc * num_cells_y) / num_yprocs ) + 1
    partition_height = ( ((yproc+1)*num_cells_y) / num_yprocs ) - partition_y_pos + 1

  end subroutine xy_decomposition


  !> @brief Defensive checks common for all decomposition strategies into a
  !         rectangular grid
  !> @param[in] num_cells_x The panel's size in the x direction
  !> @param[in] num_cells_y The panel's size in the y direction
  !> @param[in] num_xprocs  The number of partitions in the x direction
  !> @param[in] num_yprocs  The number of partitions in the y direction
  !> @param[in] panel_ranks The number of ranks the panel is to be split into
  !> @param[in] any_maps    Whether any mesh maps exist for this mesh
  subroutine xy_defensive_checks( num_cells_x, &
                                  num_cells_y, &
                                  num_xprocs,  &
                                  num_yprocs,  &
                                  panel_ranks, &
                                  any_maps )
    implicit none

    integer(i_def), intent(in) :: num_cells_x, &
                                  num_cells_y, &
                                  num_xprocs,  &
                                  num_yprocs,  &
                                  panel_ranks
    logical       , intent(in) :: any_maps

    if ( num_xprocs <=0 .or. num_yprocs <= 0 ) then
      write(log_scratch_space, "(a,i0,a,i0,a)") &
        "Number of x processes ", num_xprocs, " and y processes ", num_yprocs, &
        " must both be strictly positive."
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    if ( num_cells_x < num_xprocs ) then
      write(log_scratch_space, "(a,i0,a,i0)") &
        "Number of cells in x direction ", num_cells_x, &
        " must be greater than number of processes in x direction ", num_xprocs
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    if ( num_cells_y < num_yprocs ) then
      write(log_scratch_space, "(a,i0,a,i0)") &
        "Number of cells in y direction ", num_cells_y, &
        " must be greater than number of processes in y direction ", num_yprocs
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Equal divisions are only required if there are maps between meshes
    if (any_maps) then
      if ( mod(num_cells_x, num_xprocs) /= 0 ) then
        write(log_scratch_space, "(a,i0,a,i0)") "Requested number of ranks in x direction ", num_xprocs, &
          " must divide panel x dimension ", num_cells_x
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if

      if ( mod(num_cells_y, num_yprocs) /= 0 ) then
        write(log_scratch_space, "(a,i0,a,i0)") "Requested number of ranks in y direction ", num_yprocs, &
          " must divide panel y dimension ", num_cells_y
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end if

    if ( num_xprocs * num_yprocs /= panel_ranks ) then
      write(log_scratch_space, "(a,i0,a,i0)") "Requested number of partitions ", num_xprocs * num_yprocs, &
        " must equal available number of ranks per panel ", panel_ranks
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine xy_defensive_checks

  !> @brief Helper function for generating nonuniform partitions arranged in
  !         columns of unequal numbers of partitions
  !> @details Heights of partitions are always divisible by the mapping_factor
  !>          and differ by no more than the mapping_factor. This ensures all
  !>          cells with intermesh maps between them are on the same partition.
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    num_xprocs       The number of partitions in the x direction
  !> @param[in]    mapping_factor   The ratio between this and coarsest mapped mesh
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine nonuniform_decomposition( relative_rank,    &
                                       panel_ranks,      &
                                       num_cells_x,      &
                                       num_cells_y,      &
                                       num_xprocs,       &
                                       mapping_factor,   &
                                       part_width,       &
                                       part_height,      &
                                       part_x_pos,       &
                                       part_y_pos )
    implicit none

    integer(i_def), intent(in)    :: relative_rank,   &
                                     panel_ranks,     &
                                     num_cells_x,     &
                                     num_cells_y,     &
                                     num_xprocs,      &
                                     mapping_factor
    integer(i_def), intent(inout) :: part_width,      &
                                     part_height,     &
                                     part_x_pos,      &
                                     part_y_pos

    integer(i_def) :: mp_num_cells_y, column_idx, row_idx, rows
    real(r_def)    :: avg_rows

    ! Reduce cells by mapping factor
    mp_num_cells_y = num_cells_y / mapping_factor

    ! Column and row idx are 0-ordered
    ! Decide which column this part is in, rounding favours earlier columns
    column_idx = ( ( relative_rank - 1 ) * num_xprocs ) / panel_ranks
    ! Decide which row this part is in within its column
    row_idx = int( &
                mod( real((relative_rank - 1), r_def), real(panel_ranks, r_def) &
                / real(num_xprocs, r_def) )                                     &
              , i_def)

    ! Determine number of rows in this column
    ! Ceiling here counters use of floor (implied in int division) in column_idx
    avg_rows = real(panel_ranks, r_def) / real(num_xprocs, r_def)
    rows = ceiling( avg_rows * ( column_idx + 1 ) ) &
         - ceiling( avg_rows * column_idx )

    part_width = num_cells_x / num_xprocs
    part_x_pos = 1 + part_width * column_idx
    part_y_pos = 1 + ( ( row_idx * mp_num_cells_y ) / rows ) * mapping_factor
    ! Height is difference in start positions of adjacent partitions,
    ! or implied adjacent in case of top partition
    part_height = ( ( (row_idx + 1) * mp_num_cells_y ) / rows ) * mapping_factor &
                - part_y_pos + 1

  end subroutine nonuniform_decomposition


  !> @brief Calculate the ratio in resolution between this and the coarsest mesh
  !>        to align partitions for mapped grids
  !> @param[in] global-mesh_collection The global mesh collection
  !> @param[in] global_mesh            The global mesh to calculate the factor for
  function calc_mapping_factor( global_mesh_collection, global_mesh ) result(mp)
    implicit none

    type(global_mesh_collection_type), intent(in) :: global_mesh_collection
    type(global_mesh_type), intent(in), pointer :: global_mesh

    integer(i_def) :: mp

    type(global_mesh_type), pointer :: comparison_global_mesh
    integer(i_def) :: this_panel_width, shortest_panel_width, n_meshes, i

    this_panel_width = calc_panel_width(global_mesh)

    n_meshes = global_mesh_collection%n_meshes()

    shortest_panel_width = huge(0_i_def)

    do i = 1, n_meshes
      comparison_global_mesh => global_mesh_collection%get_mesh_by_id(i)
      if ( associated(comparison_global_mesh) ) then
        shortest_panel_width = min( shortest_panel_width, &
                                    calc_panel_width(comparison_global_mesh))
      end if
    end do

    ! If no meshes were found, return 1. This is relevant for JEDI where mesh
    ! initialisation is run twice.
    if ( shortest_panel_width < huge(0_i_def) ) then
      mp = this_panel_width / shortest_panel_width
    else
      mp = 1
    end if

  end function calc_mapping_factor


  !> @brief Calculate the width of the mesh panel. On a spherical mesh this is
  !>        the C number.
  !> @param[in] gloabl_mesh The mesh to calculate the panel width of
  function calc_panel_width( global_mesh ) result(panel_edge_ncells_x)
    use reference_element_mod, only : W, E

    implicit none

    type(global_mesh_type), intent(in), pointer :: global_mesh

    integer(i_def) :: void_cell    ! Cell id that marks the cell as a cell outside of the partition.
    integer(i_def) :: w_cell       ! The id of a cell on the western edge of the domain
    integer(i_def) :: cell_next(4) ! The cells around the cell being queried
    integer(i_def) :: cell_next_e  ! The cell to the east of the cell being queried
    logical :: periodic_xy(2)      ! Is mesh periodic in the x/y-axes
    logical :: valid_for_global_model

    integer(i_def) :: panel_edge_ncells_x  ! number of cells across a panel of the input mesh in x-direction
    integer(i_def) :: npanels

    valid_for_global_model = ( global_mesh%is_topology_periodic() .and. &
                               global_mesh%is_coord_sys_ll()      .and. &
                               global_mesh%is_geometry_spherical() )

    if ( valid_for_global_model ) then
      ! Treat as panelled sphere mesh

      ! Number of panels is sometimes unset (0) on cubesphere meshes, particularly
      ! some test cases. Just assume 6 for these.
      npanels = global_mesh%get_npanels()
      if ( npanels < 1 ) npanels = 6

      ! Generate the "C" number (as in Cnn) from global mesh data
      ! Assume the panels on the globe are square
      panel_edge_ncells_x = nint(sqrt(real( global_mesh%get_ncells() / npanels )))

    else
      ! Treat as a single panel mesh.

      ! First find a cell on the west edge of the domain
      ! If periodic, cell id 1 can be used as mesh conectivity loops round
      w_cell = 1

      void_cell   = global_mesh%get_void_cell()
      periodic_xy = global_mesh%get_mesh_periodicity()
      if ( .not. periodic_xy(1) ) then
        ! If not periodic in E-W direction then walk West until you reach mesh
        ! edge defined by the void cell.
        call global_mesh%get_cell_next(w_cell,cell_next)
        do while (cell_next(W) /= void_cell)
          w_cell = cell_next(W)
          call global_mesh%get_cell_next(w_cell,cell_next)
        end do
      end if

      ! Work out number of cells in x direction
      panel_edge_ncells_x = 1
      ! Starting at the West edge of the mesh, walk East until you reach either
      ! the cell you started at (periodic) or a void cell (LAM)
      ! - this determines the number of cells in the x direction
      call global_mesh%get_cell_next(w_cell,cell_next)
      cell_next_e = cell_next(E)
      do while (cell_next_e /= w_cell .and. cell_next_e /= void_cell)
        panel_edge_ncells_x = panel_edge_ncells_x + 1
        call global_mesh%get_cell_next(cell_next_e, cell_next)
        cell_next_e = cell_next(E)
      end do

    end if

  end function calc_panel_width

end module panel_decomposition_mod