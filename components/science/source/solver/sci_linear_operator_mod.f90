!-------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------

!> @brief Abstract base type for linear operator
!> @details Implements an abstract class for a linear operator which
!! defines an interface for the linear operator application y = A.x

module sci_linear_operator_mod

  use config_mod,         only: config_type
  use function_space_mod, only: function_space_type
  use vector_mod,         only: abstract_vector_type

  implicit none
  private

  !>@brief Abstract linear operator type to define interfaces for the solver API
  type, public, abstract :: abstract_linear_operator_type
     private
   contains
     procedure (apply_interface), deferred :: apply
  end type abstract_linear_operator_type

  abstract interface
     !> @brief Abstract interface defined for the apply method.
     !> @param[in] self a linear operator
     !> @param[in] config Application configuration object
     !> @param[in] x a vector the linear operator is applied to
     !> @param[inout] y a vector, the result.
     subroutine apply_interface(self, config, x, y)
       import :: abstract_linear_operator_type
       import :: abstract_vector_type
       import :: config_type
       class(abstract_linear_operator_type), intent(inout) :: self
       type(config_type),                    intent(in)    :: config
       class(abstract_vector_type),          intent(in)    :: x
       class(abstract_vector_type),          intent(inout) :: y
     end subroutine apply_interface
  end interface

end module sci_linear_operator_mod
