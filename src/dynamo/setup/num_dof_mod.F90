!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Computes the global and local number of dofs for the 4 element spaces
!>        W0..W3
!>
module num_dof_mod

  use mesh_generator_mod, only : nface_g,nedge_g,nvert_g

contains 

  !> Compute the local and global number of dofs.
  !>
  !> @param ncells Number of cells.
  !> @param nlayers Number of vertical cells.
  !> @param k Order of RT space ( = 0 for lowest order )
  !>
  subroutine num_dof_init( ncells, nlayers, k, w_unique_dofs, w_dof_entity )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    integer, intent( in ) :: ncells
    integer, intent( in ) :: nlayers
    integer, intent( in ) :: k

    integer, intent( out ) :: w_unique_dofs(4,2) ! there are 4 vspaces
    integer, intent( out ) :: w_dof_entity(4,0:3)

    ! numbers of dofs in each space for an element 
    integer :: nw0,nw0_cell,nw0_face,nw0_edge,nw0_vert,            &
               nw1,nw1_cell,nw1_face,nw1_edge,                     &
               nw2,nw2_cell,nw2_face,                              &
               nw3,nw3_cell

    ! global numbers of unique dofs
    integer :: nw0_g, nw1_g, nw2_g, nw3_g

    integer :: ndof_entity_w0(0:3), ndof_entity_w1(0:3),           &
               ndof_entity_w2(0:3), ndof_entity_w3(0:3)

    ! local values
    nw0 = (k+2)*(k+2)*(k+2)
    nw0_cell = k*k*k
    nw0_face = k*k
    nw0_edge = k
    nw0_vert = 1

    nw1 = 3*(k+2)*(k+2)*(k+1)
    nw1_cell = 3*k*k*(k+1)
    nw1_face = 2  *k*(k+1)
    nw1_edge =       (k+1)

    nw2  = 3*(k+2)*(k+1)*(k+1)
    nw2_cell = 3*k*(k+1)*(k+1)
    nw2_face =     (k+1)*(k+1)

    nw3 = (k+1)*(k+1)*(k+1)
    nw3_cell = nw3

    ! global numbers of dofs per function space
    nw3_g = ncells*nlayers*nw3_cell
    nw2_g = ncells*nlayers*nw2_cell + nface_g*nw2_face
    nw1_g = ncells*nlayers*nw1_cell + nface_g*nw1_face + nedge_g*nw1_edge
    nw0_g = ncells*nlayers*nw0_cell + nface_g*nw0_face + nedge_g*nw0_edge + nvert_g*nw0_vert

    ! populate the returned arrays
    w_unique_dofs(1,1) = nw0_g
    w_unique_dofs(2,1) = nw1_g
    w_unique_dofs(3,1) = nw2_g
    w_unique_dofs(4,1) = nw3_g

    w_unique_dofs(1,2) = nw0
    w_unique_dofs(2,2) = nw1
    w_unique_dofs(3,2) = nw2
    w_unique_dofs(4,2) = nw3

    ! Number of dofs per mesh entity for each space
    ndof_entity_w0(:) = (/ nw0_vert, nw0_edge, nw0_face, nw0_cell /)
    ndof_entity_w1(:) = (/ 0       , nw1_edge, nw1_face, nw1_cell /)
    ndof_entity_w2(:) = (/ 0       , 0       , nw2_face, nw2_cell /)
    ndof_entity_w3(:) = (/ 0       , 0       , 0       , nw3_cell /)

    !populate the returned arrays
    w_dof_entity(1,:) = ndof_entity_w0(:)
    w_dof_entity(2,:) = ndof_entity_w1(:)
    w_dof_entity(3,:) = ndof_entity_w2(:)
    w_dof_entity(4,:) = ndof_entity_w3(:)

    ! diagnostic output
    write( log_scratch_space, '(A, I0, A, I0)' ) &
        'ncells = ', ncells, ', nlayers = ', nlayers
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( '   space     |   W0   |   W1   |   W2   |   W3   |', &
                    LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'global dof    ', nw0_g, '   ', nw1_g, '   ', nw2_g, '   ', nw3_g
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'local dof     ', nw0, '   ', nw1, '   ', nw2, '   ', nw3
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof in volume ', nw0_cell, '   ', nw1_cell, '   ', nw2_cell, &
        '   ', nw3
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on face   ', nw0_face, '   ', nw1_face, '   ', nw2_face, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on edge   ', nw0_edge, '   ', nw1_edge, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on vert   ', nw0_vert, '   ', 0, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine num_dof_init

end module num_dof_mod
