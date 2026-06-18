!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Abstract base type for linear operator which can be coarsened
!> @details Implements an abstract class for a linear operator which
!> defines an interface for the linear operator application y = A.x and which
!> can be coarsened to the next level of the multigrid hierarchy.

module sci_hierarchical_linear_operator_mod

  use sci_linear_operator_mod, only : abstract_linear_operator_type
  use config_mod, only: config_type

  implicit none

  private

  !>@brief Abstract linear operator type which can be coarsened
  type, public, abstract, extends(abstract_linear_operator_type) :: &
                                  abstract_hierarchical_linear_operator_type
   contains
     procedure (coarsen_interface), deferred :: coarsen
  end type abstract_hierarchical_linear_operator_type

  abstract interface
     !> Abstract interface defined for the coarsening linear operator to the
     !> next level of the multigrid hierarchy.
     !>
     !>@param[in] self a hierarchical linear operator
     !>@param[in] fs_coarse coarse level function space
     !>@param[inout] coarse_operator coarsened version on the next multigrid level
     subroutine coarsen_interface(self, config, coarse_operator)
       import :: abstract_hierarchical_linear_operator_type
       import :: config_type
       class(abstract_hierarchical_linear_operator_type), intent(in) :: self
       type(config_type), intent(in) :: config
       class(abstract_hierarchical_linear_operator_type), allocatable, intent(inout) :: coarse_operator
     end subroutine coarsen_interface
  end interface

end module sci_hierarchical_linear_operator_mod
