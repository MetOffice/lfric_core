!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various mesh details.

!> @details information about the mesh is held in  here.
module mesh_mod
use constants_mod,  only : r_def
implicit none

integer :: num_cells
integer :: num_cells_x, num_cells_y
integer :: num_layers
integer :: element_order
logical :: l_spherical

integer :: w_unique_dofs(4,2)
integer :: w_dof_entity(4,0:3)    

real(kind=r_def)  :: dx, dy, dz

end module mesh_mod

