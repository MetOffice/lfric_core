!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A dummy set up routine for the mesh and elements
!> @details A dummy routine to set up the reference element, domain size  and compute the grid
!> This is to be replaced at some point by the correct preprocesser/mesh reading routines
!> @param[in] filename the filename of the grid file to read in
!> @param[out] num_cells the number of cells in a horizontal layer
!> @param[out] num_layers the number of vertical layers
!> @param[out] element_order the polynomial order of the elmenets
!> @param[out] w_unique_dofs the number of unique dofs for the function spaces
!> @param[out] w_dof_entity the number of dofs for each function space on each grid entity

module dummy_read_setup_mod

use constants_mod, only : r_def

implicit none

contains

subroutine dummy_read_setup(filename ,num_cells, num_layers, element_order, &
     w_unique_dofs,w_dof_entity)

  use mesh_generator_mod,    only : mesh_generator_init,        &
                                    mesh_generator_cubedsphere, &
                                    mesh_generator_biperiodic,  &
                                    mesh_connectivity
  use reference_element_mod, only : reference_cube
  use num_dof_mod,           only : num_dof_init
  ! dummy read the setup file, or at least the header
  implicit none
  character(*), intent(in)    :: filename 
  integer , intent(out)       :: num_cells,  num_layers, element_order
  integer, intent(out)        :: w_unique_dofs(4,2)
  integer, intent(out)        :: w_dof_entity(4,0:3)
  logical                     :: l_spherical
  real(kind=r_def), parameter :: delta = 1.0_r_def

  ! no file or file format, just coded for now 
  ! Bi-linear plane, 3x3x3 
  num_cells     = 9
  num_layers    = 3
  element_order = 0
  l_spherical   = .false.
  
! Setup reference cube  
  call reference_cube()
! Initialise mesh
  call mesh_generator_init(num_cells,num_layers)
! Genereate mesh  
  if ( l_spherical ) then
    call mesh_generator_cubedsphere(filename,num_cells,num_layers,delta)
  else
    call mesh_generator_biperiodic(num_cells,3,3,num_layers,delta,delta,delta)
  end if
! Extend connectivity ( cells->faces, cells->edges )  
  call mesh_connectivity(num_cells)    
! initialise numbers of dofs    
  call num_dof_init(num_cells,num_layers,element_order,w_unique_dofs,w_dof_entity)
    
  return
  
end subroutine dummy_read_setup

end module dummy_read_setup_mod
