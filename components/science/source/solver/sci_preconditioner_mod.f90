!-------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------

!> @brief Abstract base type for preconditioner
!>
!> @details Implements an abstract class for a preconditioner which
!>          defines an interface for the preconditioner application y = P^{-1}.x

module sci_preconditioner_mod
  use vector_mod,    only : abstract_vector_type

  implicit none
  private

  !>@brief Abstract preconditioner type for the solver API
  type, public, abstract :: abstract_preconditioner_type
     private
   contains
     procedure (apply_interface), deferred :: apply
  end type abstract_preconditioner_type

  abstract interface
     !> abstract interface defined for the apply procedure of a preconditioner
     !! y = P.x
     !> @param[in] self a preconditioner
     !> @param[in]    x a vector that the preconditioner is applied to.
     !> @param[inout] y a vector, the result.
     subroutine apply_interface(self, x, y)
       import :: abstract_vector_type
       import :: abstract_preconditioner_type
       class(abstract_preconditioner_type), intent(inout) :: self
       class(abstract_vector_type),         intent(in)    :: x
       class(abstract_vector_type),         intent(inout) :: y
     end subroutine apply_interface
  end interface

end module sci_preconditioner_mod
