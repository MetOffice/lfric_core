!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to sample a flux at nodal points: F = u*q

!> @details Samples a flux field F at nodal points where the flux
!> is defined as the product of a velocity u and a scalar q
!> The quadrature rule is overloaded to give the nodal points 
!> of the flux space
module sample_flux_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    W0, W2, ANY_SPACE_1,                     &
                                    GH_BASIS, GH_DIFF_BASIS,                 &
                                    CELLS, GH_EVALUATOR, EVALUATOR
use constants_mod,           only : r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_flux_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_1)                      &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: evaluator_shape = EVALUATOR
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass ::sample_flux_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface sample_flux_kernel_type
   module procedure sample_flux_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public sample_flux_code
contains

type(sample_flux_kernel_type) function sample_flux_kernel_constructor() result(self)
  return
end function sample_flux_kernel_constructor

!> @brief Kernel to sample a flux at nodal points: F = u*q
!! @param[in] nlayers Number of layers
!! @param[in] ndf_f Number of degrees of freedom per cell for w2
!! @param[in] undf_f Number of unique degrees of freedom for w2
!! @param[in] map_f Dofmap for the cell at the base of the column for w2
!! @param[inout] flux Field to contain the right hand side to be computed
!! @param[in] multiplicity How many times the dof has been visited in total
!! @param[in] u Advecting wind
!! @param[in] ndf_q Number of degrees of freedom per cell for the field to be advected
!! @param[in] undf_q  Number of unique degrees of freedom for the advected field
!! @param[in] map_q Dofmap for the cell at the base of the column for the field to be advected
!! @param[in] basis_q Basis functions evaluated at gaussian quadrature points 
!! @param[in] q Advected field
subroutine sample_flux_code(nlayers,                                           &
                            flux, multiplicity, u, q,                          &
                            ndf_f, undf_f, map_f,                              &
                            ndf_q, undf_q, map_q, basis_q                      &
                            )
                           
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_f, ndf_q, undf_f, undf_q
  integer, dimension(ndf_f), intent(in) :: map_f
  integer, dimension(ndf_q), intent(in) :: map_q
  real(kind=r_def), dimension(1,ndf_q,ndf_f), intent(in)    :: basis_q
  real(kind=r_def), dimension(undf_f),        intent(inout) :: flux
  real(kind=r_def), dimension(undf_f),        intent(in)    :: u, multiplicity
  real(kind=r_def), dimension(undf_q),        intent(in)    :: q

  !Internal variables
  integer               :: df, df_q, k, loc
  
  real(kind=r_def), dimension(ndf_q) :: q_cell
  real(kind=r_def)                   :: q_at_node

  do k = 0, nlayers-1
    do df_q = 1, ndf_q
      q_cell(df_q) = q( map_q(df_q) + k )
    end do    
    do df = 1, ndf_f
      q_at_node = 0.0_r_def
      do df_q = 1,ndf_q
        q_at_node = q_at_node + q_cell(df_q)*basis_q(1,df_q,df)
      end do
      loc = map_f(df) + k
      flux( loc ) = flux( loc ) + u( loc )*q_at_node/multiplicity( loc )
    end do
  end do

end subroutine sample_flux_code

end module sample_flux_kernel_mod
