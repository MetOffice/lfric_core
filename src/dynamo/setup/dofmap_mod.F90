!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A module that holds dofmaps for the four element spaces 
!>
!> @detail The dofmaps for the four element spaces are stored in this module. The 
!>         module also contains the code to calculate the dofmaps. This will eventually
!>         be replaced with code that reads them in from a file.

!-------------------------------------------------------------------------------
! Computes the dofmaps for the 4 element spaces given grid connectivity information
! requires: list of cell next to current cell
!           list of vertices on this cell
!-------------------------------------------------------------------------------
module dofmap_mod

use num_dof_mod
use reference_element_mod
use mesh_generator_mod, only: nedge_h_g, nvert_h_g, face_on_cell, edge_on_cell, vert_on_cell

implicit none

!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W0 function space over the bottom level of the domain.
integer, allocatable :: w0_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W1 function space over the bottom level of the domain.
integer, allocatable :: w1_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W2 function space over the bottom level of the domain.
integer, allocatable :: w2_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W3 function space over the bottom level of the domain.
integer, allocatable :: w3_dofmap(:,:)

!> A two dim integer array which holds the orientation data for the
!> W0 function space
integer, allocatable :: w0_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W1 function space
integer, allocatable :: w1_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W2 function space
integer, allocatable :: w2_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W3 function space
integer, allocatable :: w3_orientation(:,:)

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains 

!> Subroutine to get the dofmap and copy is into the function space
!> @param[in] nlayers the number of vertical layers
!> @param[in] function_space the function space
!> @param[in] ndf_entity the number of dofs on each grid entity
subroutine get_dofmap(nlayers, w_dof_entity, &
                      ncell, w_unique_dofs )
  
  implicit none
  
  integer, intent(in)       :: nlayers
  integer, intent(in)       :: w_dof_entity(4,0:3)
  integer, intent(in)       :: ncell
  integer, intent(in)       :: w_unique_dofs(4,2) 

  allocate( w0_dofmap(0:ncell,w_unique_dofs(1,2)) )
  allocate( w1_dofmap(0:ncell,w_unique_dofs(2,2)) )
  allocate( w2_dofmap(0:ncell,w_unique_dofs(3,2)) )
  allocate( w3_dofmap(0:ncell,w_unique_dofs(4,2)) )
  
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(1,2), w_dof_entity(1,:), w0_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(2,2), w_dof_entity(2,:), w1_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(3,2), w_dof_entity(3,:), w2_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(4,2), w_dof_entity(4,:), w3_dofmap)

end subroutine get_dofmap

!> Subroutine to compute the dofmap based upon grid connectivities for a function space
!> @param[in] ncells the number of horizontal cells
!> @param[in] nlayers the number of vertical layers
!> @param[in] ndof_sum the total number of dofs associated with a single cell
!> @param[in] ndf_entity the number of dofs on each grid entity
!> @param[out] dofmap
subroutine dofmap_populate(ncells,nlayers,ndof_sum,ndof_entity,dofmap)

  integer, intent(in) :: ncells, nlayers

! number of dofs per entity for this space
  integer, intent(in) :: ndof_entity(0:3)
! total number of dofs associated with each cell  
  integer, intent(in) :: ndof_sum

! output dofmap for this space
  integer, intent(out) :: dofmap(0:ncells,ndof_sum)

! loop counters
  integer :: i, j, k

  integer :: id, id0, jd, jdp, dof_idx
! Number of entities for a single layer  
  integer :: nvert_layer, nedge_layer, nface_layer

! entity dofmaps
  integer, allocatable :: dofmap_d0(:,:), dofmap_d1(:,:), dofmap_d2(:,:), dofmap_d3(:,:)

! dofmaps for a 3D horizontal layer
  nvert_layer = 2*nvert_h_g 
  nedge_layer = 2*nedge_h_g + nvert_h_g
  nface_layer = nedge_h_g + 2*ncells
  
  if ( ndof_entity(0) > 0 ) then
    allocate( dofmap_d0(nvert_layer,ndof_entity(0)) )
  else
    allocate( dofmap_d0(nvert_layer,1) )
  end if
  if ( ndof_entity(1) > 0 ) then  
    allocate( dofmap_d1(nedge_layer,ndof_entity(1)) )
  else
    allocate( dofmap_d1(nedge_layer,1) )
  end if  
  if ( ndof_entity(2) > 0 ) then  
    allocate( dofmap_d2(nface_layer,ndof_entity(2)) )
  else
    allocate( dofmap_d2(nface_layer,1) )
  end if
  if ( ndof_entity(3) > 0 ) then    
    allocate( dofmap_d3(ncells,     ndof_entity(3)) )
  else
    allocate( dofmap_d3(ncells,1) )
  end if 

! initialise entity dofmaps
  dofmap_d0(:,:) = 0
  dofmap_d1(:,:) = 0
  dofmap_d2(:,:) = 0
  dofmap_d3(:,:) = 0

! assume we have all possible global connectivity information
! in practice this requires connectivity
! (3,2) -> faces on cells
! (3,1) -> edges on cells
! (3,0) -> vertices on cells

  id = 1
! loop over 3 entities (cells)
  do i=1,ncells
    dof_idx = 1
! assign dofs for connectivity (3,3) (dofs in cell)
    do j=1,ndof_entity(3)
      dofmap_d3(i,j) = id
      dofmap(i,dof_idx) = dofmap_d3(i,j)
      id = id + nlayers
      dof_idx = dof_idx + 1
    end do
  
! assign dofs for connectivity (3,2) (dofs on faces)
    do j=1,nfaces_h
      jd = face_on_cell(i,j) 
      if ( dofmap_d2(jd,1) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(jd,k) = id        
          id = id + nlayers
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(i,dof_idx) = dofmap_d2(jd,k)
        dof_idx = dof_idx + 1
      end do
    end do
    id0 = id
    do j=nfaces_h+1,nfaces
      jd = face_on_cell(i,j) 
      if ( dofmap_d2(jd,1) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(jd,k) = id        
          id = id + nlayers + 1
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(i,dof_idx) = dofmap_d2(jd,k)
        dof_idx = dof_idx + 1
      end do
      if (j==nfaces_h+1) then
        id = id0 + 1
      else
        id = id - 1
      end if
    end do
! assign dofs for connectivity (3,1) (dofs on edges)  
    do j=1,nedges_h
      jd  = edge_on_cell(i,j)   
      jdp = edge_on_cell(i,j+nedges-nedges_h)  
      if ( dofmap_d1(jd,1) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(jd,k)  = id
          dofmap_d1(jdp,k) = id+1
          id = id + nlayers + 1
        end do
      end if
    end do
    do j=5,8
      jd  = edge_on_cell(i,j) 
      if ( dofmap_d1(jd,1) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(jd,k)  = id
          id = id + nlayers 
        end do
      end if
    end do
    do j=1,nedges
      jd  = edge_on_cell(i,j) 
      do k=1,ndof_entity(1)
        dofmap(i,dof_idx) = dofmap_d1(jd,k)
        dof_idx = dof_idx + 1  
      end do
    end do 
! assign dofs for connectivity (3,0) (dofs on verts)    
    do j=1,nverts_h
      jd  = vert_on_cell(i,j)
      jdp = vert_on_cell(i,j+nverts_h)
      if ( dofmap_d0(jd,1) == 0 ) then
        do k=1,ndof_entity(0)
          dofmap_d0(jd, k)  = id
          dofmap_d0(jdp,k)  = id + 1
          id = id + nlayers + 1  
        end do
      end if
    end do
    do j=1,nverts
      jd  = vert_on_cell(i,j)
      do k=1,ndof_entity(0)
        dofmap(i,dof_idx) = dofmap_d0(jd,k) 
        dof_idx = dof_idx + 1  
      end do
    end do
  end do
  
  dofmap(0,:) = 0

  if (allocated(dofmap_d0) ) deallocate( dofmap_d0 )
  if (allocated(dofmap_d1) ) deallocate( dofmap_d1 )
  if (allocated(dofmap_d2) ) deallocate( dofmap_d2 )
  if (allocated(dofmap_d3) ) deallocate( dofmap_d3 )

end subroutine dofmap_populate

!> Subroutine to compute the orientation of vectors
!> @param[in] ncell the number of horizontal cells
!> @param[in] w_unique_dofs The number of dofs in each function space
subroutine get_orientation(ncell,w_unique_dofs)
!-----------------------------------------------------------------------------
! Subroutine to read orientation
!-----------------------------------------------------------------------------
  
  implicit none

  integer, intent(in) :: ncell
  integer, intent(in) :: w_unique_dofs(4,2)

  allocate( w0_orientation(0:ncell,w_unique_dofs(1,2)) )
  allocate( w1_orientation(0:ncell,w_unique_dofs(2,2)) )
  allocate( w2_orientation(0:ncell,w_unique_dofs(3,2)) )
  allocate( w3_orientation(0:ncell,w_unique_dofs(4,2)) )

  call orientation_populate(ncell, w_unique_dofs(1,2), w0_orientation)
  call orientation_populate(ncell, w_unique_dofs(2,2), w1_orientation)
  call orientation_populate(ncell, w_unique_dofs(3,2), w2_orientation)
  call orientation_populate(ncell, w_unique_dofs(4,2), w3_orientation)

end subroutine get_orientation

!> Subroutine to compute the orientation of vectors
!> @param[in] ncells the number of horizontal cells
!> @param[in] ndof_sum the total number of dofs associated with a single cell
!> @param[out] orientation The output orientation
subroutine orientation_populate(ncells,ndof_sum,orientation)

! total number of cells
  integer, intent(in) :: ncells
! total number of dofs associated with each cell  
  integer, intent(in) :: ndof_sum
! output orientation 
  integer, intent(out) :: orientation(0:ncells,ndof_sum)

  integer :: cell, df

  do cell = 0, ncells
     do df = 1, ndof_sum
       orientation(cell,df) = 1
     end do
  end do

end subroutine orientation_populate

end module dofmap_mod
