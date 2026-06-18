!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Abstract base type for preconditioner which can be coarsned
!>
!> @details Implements an abstract class for a preconditioner which
!>          defines an interface for the preconditioner application
!>          y = P^{-1}.x and which can be coarsenend to the next level of a
!>          multigrid hierarchy

module sci_hierarchical_preconditioner_mod

  use sci_preconditioner_mod, only: abstract_preconditioner_type
  use config_mod, only: config_type

  implicit none
  private

  !>@brief Abstract preconditioner type which can be coarsened
  type, public, abstract, extends(abstract_preconditioner_type) :: &
                                  abstract_hierarchical_preconditioner_type
   contains
     procedure (coarsen_interface), deferred :: coarsen
  end type abstract_hierarchical_preconditioner_type

  abstract interface
     !> Abstract interface defined for the coarsening preconditioner to the
     !> next level of the multigrid hierarchy.
     !>
     !>@param[inout] self a hierarchical preconditioner
     !>@param[inout] other coarsened version on the next multigrid level
     subroutine coarsen_interface(self, config,other)
       import :: abstract_hierarchical_preconditioner_type
       import :: config_type
       class(abstract_hierarchical_preconditioner_type), intent(inout) :: self
type(config_type), intent(in) :: config
       class(abstract_hierarchical_preconditioner_type), allocatable, intent(inout) :: other
     end subroutine coarsen_interface
  end interface

end module sci_hierarchical_preconditioner_mod
