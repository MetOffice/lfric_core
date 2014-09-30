!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Sets up data required for the the function spaces

!> @details This code generates a mesh and determines the basis functions and
!> dofmaps. This will be replaced with code that reads this in from a mesh
!> generation and paritioning pre-processor stage.  

! There are no tests for this code as this will be replaced.

module set_up_mod

  use constants_mod,              only : r_def, str_def
  use function_space_mod,         only : function_space_type
  use reference_element_mod,      only : reference_cube

  use mesh_generator_mod,         only : mesh_generator_init,        &
                                         mesh_generator_cubedsphere, &
                                         mesh_generator_biperiodic,  &
                                         mesh_connectivity
  use num_dof_mod,                only : num_dof_init
  use basis_function_mod,         only : get_basis, &
              w0_basis, w1_basis, w2_basis, w3_basis, &
              w0_diff_basis, w1_diff_basis, w2_diff_basis, w3_diff_basis, &
              w0_nodal_coords, w1_nodal_coords, w2_nodal_coords, w3_nodal_coords

  use dofmap_mod,                 only : get_dofmap, get_orientation, &
              w0_dofmap, w1_dofmap, w2_dofmap, w3_dofmap
  use gaussian_quadrature_mod,    only : ngp_h, ngp_v
  use mass_matrices_mod,          only : mass_matrix_init
  implicit none
  
contains 

!> Generates a mesh and determines the basis functions and dofmaps (this will 
!> be replaced with code that reads the information in)
  subroutine set_up( )

    use log_mod,  only : log_event, LOG_LEVEL_INFO
    use mesh_mod, only : num_cells, num_layers, element_order, l_spherical, &
                         w_unique_dofs, w_dof_entity, dx, dy, dz,           &
                         num_cells_x, num_cells_y

    implicit none

    real(kind=r_def), parameter              :: delta = 1.0_r_def
    character(len = str_def)                 :: filename

    ! hard-coded these numbers are
    num_cells_x = 50
    num_cells_y = 4
    num_layers = 5
    element_order = 0
    l_spherical = .false.
    dx = 6000.0_r_def
    dy = 1000.0_r_def
    dz = 2000.0_r_def
    filename = 'ugrid_quads_2d.nc' 
    call log_event( "set_up: generating/reading the mesh", LOG_LEVEL_INFO )
    
    ! total number of horizontal cells ( num_cells is currently cells along one edge )   
    if ( l_spherical ) then 
      num_cells_y = num_cells_x
      num_cells = 6*num_cells_x**2
    else
      num_cells = num_cells_x*num_cells_y
    end if

!  ----------------------------------------------------------
!  Mesh generation, really a preprocessor step for reading
! -----------------------------------------------------------

    ! Setup reference cube  
    call reference_cube()
    ! Initialise mesh
    call mesh_generator_init(num_cells,num_layers)
    ! Genereate mesh  
    if ( l_spherical ) then
       call mesh_generator_cubedsphere(filename,num_cells,num_layers,dz)
    else
       call mesh_generator_biperiodic(num_cells,num_cells_x,num_cells_y,num_layers,dx,dy,dz)
    end if
    ! Extend connectivity ( cells->faces, cells->edges )  
    call mesh_connectivity(num_cells)    

! -----------------------------------------------------------
! Initialise FE elements on the mesh constructed above
! really another pre-processor step
! ----------------------------------------------------------

    ! initialise numbers of dofs    
    call num_dof_init(num_cells,num_layers,element_order,w_unique_dofs,w_dof_entity)
         
    call mass_matrix_init(w_unique_dofs(1,2),w_unique_dofs(2,2),w_unique_dofs(3,2),num_cells*num_layers)

    call log_event( "set_up: computing basis functions", LOG_LEVEL_INFO )

    ! read the values of the basis functions. 
    call get_basis( k=element_order, &
                    w_unique_dofs=w_unique_dofs,w_dof_entity=w_dof_entity )  

    call log_event( "set_up: computing the dof_map", LOG_LEVEL_INFO )
    ! compute the dof maps for each function space
    call get_dofmap(nlayers=num_layers,w_dof_entity=w_dof_entity, &
                    ncell=num_cells,w_unique_dofs=w_unique_dofs)
    
    ! compute cell local orientations for vector spaces
    call get_orientation(num_cells, w_unique_dofs)

    return

  end subroutine set_up

end module set_up_mod
