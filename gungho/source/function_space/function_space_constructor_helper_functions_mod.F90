!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!>
!> @brief Holds support routines for instantiating a function space.
!>
!> @details 
!
module function_space_constructor_helper_functions_mod

use constants_mod,         only: i_def, i_halo_index, r_def, dp_xios, IMDI
use mesh_mod,              only: mesh_type
use fs_continuity_mod,     only: W0, W1, W2, W3, Wtheta, W2V, W2H, Wchi
use reference_element_mod, only: nfaces_h, nedges_h, nverts_h                  &
                               , nfaces,   nedges,   nverts,   x_vert          &
                               , edge_on_face, tangent_to_edge, normal_to_face &
                               , select_entity_type, select_entity_all         &
                               , select_entity_theta, select_entity_w2h        &
                               , select_entity_w2v

use reference_element_data_mod, only :  select_data_entity_all, &
                                        select_data_entity_theta, &
                                        select_data_entity_w2h, &
                                        select_data_entity_w2v, &
                                        vert_coords, edge_coords, &
                                        face_coords, vol_coords

implicit none

private
public :: ndof_setup, basis_setup, dofmap_setup, levels_setup

contains

!-----------------------------------------------------------------------------
! Compute the number of dofs for various entities
!-----------------------------------------------------------------------------
!> @brief   Private routine (function_space_type construction)
!> @details Compute number of dofs for this function space. This subroutine
!>          computes the number of dofs in this function space
!>          and the number of dofs associated with various entities, primitive
!>          and composite. This is a private routine used in the construction
!>          of function spaces.
!> @param[in] mesh Mesh to define the function space on
!> @param[in] element_order Polynomial order of the function space
!> @param[in] dynamo_fs Integer enumeration of the function space
!> @param[out] ndof_vert Number dofs on each vertex
!> @param[out] ndof_edge Number of dofs on each edge
!> @param[out] ndof_face Number of dofs on each face
!> @param[out] ndof_vol Number of dofs in each volume
!> @param[out] ndof_cell Total number of dofs associated with a cell
!> @param[out] ndof_glob Total number of global dofs
!> @param[out] ndof_interior Number of dofs with no  vertical connectivity
!> @param[out] ndof_exterior Number of dofs with vertical connectivity
 subroutine ndof_setup( mesh, element_order, dynamo_fs,                         &
                       ndof_vert, ndof_edge, ndof_face, ndof_vol, ndof_cell,   &
                       ndof_glob, ndof_interior, ndof_exterior )

  ! NOTE: ndofs will be used as short hand for Number of Degrees Of Freedom
  implicit none

  ! Input
  type(mesh_type), intent(in), pointer :: mesh
  integer(i_def),  intent(in)          :: element_order
  integer(i_def),  intent(in)          :: dynamo_fs
 
  ! Output
  ! Number of dofs per ...
  integer(i_def), intent(out) :: ndof_vert     ! vertex entity
  integer(i_def), intent(out) :: ndof_edge     ! edge entity
  integer(i_def), intent(out) :: ndof_face     ! face entity
  integer(i_def), intent(out) :: ndof_vol      ! volume entity

  integer(i_def), intent(out) :: ndof_cell     ! 3D-cell entity
  integer(i_def), intent(out) :: ndof_interior ! interior entity (in vertical)
  integer(i_def), intent(out) :: ndof_exterior ! exterior entity (in vertical)
  integer(i_def), intent(out) :: ndof_glob     ! 3D-mesh (on a rank)


  ! Local variables

  ! Variables for properties of the local 3D-Mesh
  integer(i_def) :: ncells           ! No. of 2D-cells in 3D-mesh partition
  integer(i_def) :: nlayers          ! No. of layers of 3D-cells
  integer(i_def) :: nface_g          ! No. of faces
  integer(i_def) :: nedge_g          ! No. of edges
  integer(i_def) :: nvert_g          ! No. of vertices
  integer(i_def) :: nedges_per_level ! No. of edges per level

  ! Variables for Exterior-Interior topology (vertical direction)
  integer(i_def) :: nverts_exterior  ! No. of vertices per exterior entity
  integer(i_def) :: nedges_exterior  ! No. of edges    per exterior entity
  integer(i_def) :: nfaces_exterior  ! No. of faces    per exterior entity
  integer(i_def) :: nedges_interior  ! No. of edges    per interior entity
  integer(i_def) :: nfaces_interior  ! No. of faces    per interior entity

  integer(i_def) :: k

  ! Adding ndof for exterior and interior composite entities
  !
  !   ndof_exterior = ndof_edge*nedges_exterior 
  !                 + ndof_face*nfaces_exterior 
  !                 + ndof_vert*nverts_exterior
  !
  !   ndof_interior = ndof_edge*nedges_interior 
  !                 + ndof_face*nfaces_interior
  !                 + ndof_vol
  !
  ! Elements on interior/exterior cell decomposition in vertical,
  ! the horizontal faces and associated edges/vertices
  ! (i.e. top OR bottom ) are classed as exterior entities.
  ! The vertical faces/edges are classed as an interior entities.

  nverts_exterior = nverts_h
  nedges_exterior = nedges_h
  nfaces_exterior = 1

  nedges_interior = nverts_h
  nfaces_interior = nedges_h


  ! Local values
  nlayers  = mesh % get_nlayers()
  ncells   = mesh % get_ncells_2d_with_ghost()
  nface_g  = mesh % get_nfaces()
  nedge_g  = mesh % get_nedges()
  nvert_g  = mesh % get_nverts()
  nedges_per_level = mesh % get_nedges_2d()

  ndof_vert = 0
  ndof_edge = 0
  ndof_face = 0
  ndof_vol  = 0
  ndof_cell = 0
  ndof_glob = 0

  ndof_interior  = 0
  ndof_exterior  = 0

  k = element_order

  ! Possible modifications to number of dofs
  ! on edges depending on presets
  select case (dynamo_fs)

  case (W0)
    ! H1 locates dofs on the element vertices for a element order = 0,
    ! though the order for the H1 function space is k+1, i.e.
    ! linear across the element on each axis
    ndof_vert = 1
    ndof_edge = k
    ndof_face = k*k
    ndof_vol  = k*k*k
    ndof_cell = (k+2)*(k+2)*(k+2)


  case (W1)
    ! Dofs located on edges, as vectors
    ! in direction of edge.

    ! For order 0, the vector is constant along the
    ! edge, but can vary linearly normal to it.
    ndof_edge =   (k+1)
    ndof_face = 2*(k+1)*k
    ndof_vol  = 3*(k+1)*k*k
    ndof_cell = 3*(k+1)*(k+2)*(k+2)


  case (W2)
    ! Dofs are located on faces for vector fields
    ! and direction is normal to the face.
    !
    ! For order 0 the value of the vector normal to the
    ! face is constant across the face(tangential) but can
    ! varying linearly passing through the face(normal) to
    ! the next cell.
    !
    ! So linear   in normal: 1-dim, ndof = 2
    ! So constant in tangengial: 2-dim, each ndof = 1
    ! So 3 dimensions each with ndof (k+2)(k+1)(k+1)
    ndof_face =   (k+1)*(k+1)
    ndof_vol  = 3*(k+1)*(k+1)*k
    ndof_cell = 3*(k+1)*(k+1)*(k+2)


  case (W3)
    ! Order of this function space is same as base order
    ! This function space is discontinuous so all dofs are
    ! located on the cell volume, not the edges or vertices

    ! Dofs located on cell volume entities/discontinuous
    ! between cells.

    ! Number of dofs on each dimension is lowest order + 1
    ndof_vol  = (k+1)*(k+1)*(k+1)
    ndof_cell = ndof_vol


  case (WTHETA)
    nfaces_interior = 0
    ndof_face  =   (k+1)*(k+1)
    ndof_vol   = k*(k+1)*(k+1)
    ndof_cell  = (k+2)*(k+1)*(k+1)


  case (W2H)
    nfaces_exterior = 0
    ndof_face  =     (k+1)*(k+1)
    ndof_vol   = 2*k*(k+1)*(k+1)
    ndof_cell  = 2*(k+2)*(k+1)*(k+1)


  case (W2V)
    nfaces_interior = 0
    ndof_face  =     (k+1)*(k+1)
    ndof_vol   = 1*k*(k+1)*(k+1)
    ndof_cell  = 1*(k+2)*(k+1)*(k+1)

  case (WCHI)
    ndof_vol  = (k+1)*(k+1)*(k+1)
    ndof_cell = ndof_vol

  end select

  ndof_exterior = ndof_vert*nverts_exterior &
                + ndof_edge*nedges_exterior &
                + ndof_face*nfaces_exterior

  ndof_interior = ndof_edge*nedges_interior &
                + ndof_face*nfaces_interior &
                + ndof_vol  

  ! Calculated the global number of dofs on the function space
  select case (dynamo_fs)
  case (W0,W1,W2,W3,WCHI)
    ndof_glob = ncells*nlayers*ndof_vol + nface_g*ndof_face                    &
                + nedge_g*ndof_edge     + nvert_g*ndof_vert

  case (WTHETA,W2V)
    ndof_glob = ncells*nlayers*ndof_vol + ncells*(nlayers+1)*ndof_face

  case (W2H)
    ndof_glob = ncells*nlayers*ndof_vol + nedges_per_level*nlayers*ndof_face   &
              + nedge_g*ndof_edge       + nvert_g*ndof_vert
  end select

  return
end subroutine ndof_setup



!-----------------------------------------------------------------------------
! Setup basis functions for the function space
!-----------------------------------------------------------------------------
!> @brief   Private routine (function_space_type construction)
!> @details Setup arrays to for basis function generation. This subroutine
!>          computes the arrays required for "on-the-fly" basis function
!>          generation. This is a private routine used in the construction
!>          of function spaces. This routine is only valid for cube elements
!> @param[in] element_order Polynomial order of the function space
!> @param[in] dynamo_fs Enumeration of the function space
!> @param[in] ndof_vert Number dofs on each vertex
!> @param[in] ndof_cell Total number of dofs associated with a cell 
!> @param[out] basis_index Array containing index of polynomial function
!> @param[out] basis_order Polynomial order of basis function
!> @param[out] basis_vector Direction of basis for vector functions
!> @param[out] basis_x Array of nodal points of the basis functions
!> @param[out] nodal_coords 3D coordinates of zeros of the basis functions
!> @param[out] dof_on_vert_boundary Array indication if a dof is on the top or
!>            bottom boundary of a cell
subroutine basis_setup( element_order, dynamo_fs, ndof_vert,  ndof_cell,       &
                        basis_index, basis_order, basis_vector, basis_x,       &
                        nodal_coords, dof_on_vert_boundary )

  implicit none

  ! Input
  integer(i_def), intent(in) :: element_order
  integer(i_def), intent(in) :: dynamo_fs
 
  ! Number of dofs per entity
  integer(i_def), intent(in) :: ndof_vert ! ndofs per vertex
  integer(i_def), intent(in) :: ndof_cell ! ndofs per 3D-cell

  ! Output
  integer(i_def), intent(out) :: basis_index  (:,:)
  integer(i_def), intent(out) :: basis_order  (:,:)
  real(r_def),    intent(out) :: basis_vector (:,:)
  real(r_def),    intent(out) :: basis_x      (:,:,:)
  real(r_def),    intent(out) :: nodal_coords (:,:)
  integer(i_def), intent(out) :: dof_on_vert_boundary (:,:)

  integer(i_def) :: k


  integer(i_def) :: i, jx, jy, jz, poly_order, idx, j1, j2
  integer(i_def) :: j(3), j2l_edge(12,3), j2l_face(6,3), face_idx(6), edge_idx(12,2)
  integer(i_def), allocatable :: lx(:), ly(:), lz(:)
  real(kind=r_def), allocatable :: unit_vec(:,:)

  real(r_def) :: x1(element_order+2)
  real(r_def) :: x2(element_order+2)

  ! To uniquely specify a 3D tensor product basis function the folling is needed:
  ! basis_order(3): The polynomial order in the x,y,z directions
  ! basis_x(3,basis_order+1): The nodal points of the polynomials in each direction
  ! basis_index(3): The index of the nodal points array at which the basis function is unity
  ! basis_vector(3): Additionally if the function space is a vectro then a unit vector is needed.
  
  ! Although not strictly needed the nodal coordinates at which each basis
  ! function equals 1 is stored as nodal_coords
  ! A flag is also set to 0 if a basis function is associated with an entity on
  ! the top or bottom of the cell, i.e has nodal_coord(3) = 0 or 1

  k = element_order

  ! Allocate to be larger than should be needed
  allocate ( lx(3*(k+2)**3) )
  allocate ( ly(3*(k+2)**3) )
  allocate ( lz(3*(k+2)**3) )

  lx(:) = 0
  ly(:) = 0
  lz(:) = 0

  ! Positional arrays - need two, i.e quadratic and linear for RT1
  do i=1,k+2
    x1(i) = real(i-1,r_def)/real(k+1,r_def)
  end do

  if ( k == 0 ) then
    x2(1) = 0.5_r_def
  else
    do i=1,k+1
      x2(i) = real(i-1,r_def)/real(k,r_def)
    end do
  end if

  if ( k == 0 ) x2(1) = 0.5_r_def
  ! This value isn't needed and is always multipled by 0
  x2(k+2) = 0.0_r_def

  ! Some look arrays based upon reference cube topology
  ! index of nodal points for dofs located on faces.
  ! Faces are defined as having one coodinate fixed, 
  ! i.e. for face 1 x = 0 for all points on the face
  ! and for face 4 y = 1 for all points on the face
  ! This array give the index for the fixed coordinate for each face.
  ! If a face has fixed coordinate = 0 the index is 1
  ! If a face has fixed coordinate = 1 the index is k+2
  face_idx = (/ 1, 1, k+2, k+2, 1, k+2 /)

  ! index of nodal points for dofs located on edges
  ! edges are defined as having two coodinates fixed, 
  ! i.e. for edge 1 x = 0 & z = 0 for all points on the edge
  ! and for edge 6 x = 1 y = 0 for all points on the edge
  ! These arrays give the index for the two fixed coordinates for each edge.
  ! If an edge has fixed coordinate = 0 the index is 1
  ! If an edge has fixed coordinate = 1 the index is k+2
  edge_idx(:,1) = (/ 1, 1, k+2, k+2, 1, k+2, k+2, 1,   1,   1,   k+2, k+2 /)
  edge_idx(:,2) = (/ 1, 1, 1,   1,   1, 1,   k+2, k+2, k+2, k+2, k+2, k+2 /)

  ! Each dof living on a face or edge will have its index defined by three
  ! integers (j1, j2, j3) where:
  !  for faces one j will be the face index and the other two can vary 
  !  for edges two j's will be the edge indices and the final one can vary
  ! These j's need to be converted to the indices lx ,ly, lz
  ! For faces the first value of j2l is the l that corresponds to the 
  ! constant coordinate, so for face 1 lx = j3, ly = j2 and lz = j1/
  ! for edge 1: lx = j2, ly = j1, and lz = j3  
  j2l_face(1,:) = (/ 3, 2, 1 /)
  j2l_face(2,:) = (/ 2, 3, 1 /)
  j2l_face(3,:) = (/ 3, 2, 1 /)
  j2l_face(4,:) = (/ 2, 3, 1 /)
  j2l_face(5,:) = (/ 1, 2, 3 /)
  j2l_face(6,:) = (/ 1, 2, 3 /)

  j2l_edge(1 ,:) = (/ 2, 1, 3 /)
  j2l_edge(2 ,:) = (/ 1, 2, 3 /)
  j2l_edge(3 ,:) = (/ 2, 1, 3 /)
  j2l_edge(4 ,:) = (/ 1, 2, 3 /)
  j2l_edge(5 ,:) = (/ 2, 3, 1 /)
  j2l_edge(6 ,:) = (/ 2, 3, 1 /)
  j2l_edge(7 ,:) = (/ 2, 3, 1 /)
  j2l_edge(8 ,:) = (/ 2, 3, 1 /)
  j2l_edge(9 ,:) = (/ 2, 1, 3 /)
  j2l_edge(10,:) = (/ 1, 2, 3 /)
  j2l_edge(11,:) = (/ 2, 1, 3 /)
  j2l_edge(12,:) = (/ 1, 2, 3 /)

  ! Array to flag vertices on the top or bottom boundaries
  ! If dof j is on the bottom boundary then  dof_on_vert_boundary(j,1) = 0
  ! If dof j is on the top boundary then  dof_on_vert_boundary(j,2) = 0
  dof_on_vert_boundary(:,:) = 1

  ! Allocate arrays to allow on the fly evaluation of basis functions
  select case (dynamo_fs)
  case (W1, W2, W2H, W2V)
    allocate( unit_vec(ndof_cell,3) )
  end select


  select case (dynamo_fs)
  case (W0)

    !---------------------------------------------------------------------------
    ! Section for test/trial functions of CG spaces
    !---------------------------------------------------------------------------
    poly_order = k+1

    ! Compute indices of functions
    idx = 1

    ! ===============================
    ! dofs in volume
    ! ===============================
    do jz=2, k+1
      do jy=2, k+1
        do jx=2, k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          idx = idx + 1
        end do
      end do
    end do

    ! ===============================
    ! dofs on faces
    ! ===============================
    do i=1, nfaces
      do j1=2, k+1
        do j2=2, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          idx = idx + 1
        end do
      end do
    end do

    ! ===============================
    ! dofs on edges
    ! ===============================
    do i=1, nedges
      do j1=2, k+1
        j(1)    = j1
        j(2)    = edge_idx(i,1)
        j(3)    = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        idx     = idx + 1
      end do 
    end do

    ! ===============================
    ! dofs on vertices
    ! ===============================    
    do i=1, nverts
      do j1=1, ndof_vert
        lx(idx) = 1+(k+1)*int(x_vert(i,1))
        ly(idx) = 1+(k+1)*int(x_vert(i,2))
        lz(idx) = 1+(k+1)*int(x_vert(i,3))
        idx     = idx + 1
      end do
    end do

    do i=1, ndof_cell

      ! Explicitly for quads, as ngp_h = ngp_v * ngp_v
      nodal_coords(1,i) = x1(lx(i))
      nodal_coords(2,i) = x1(ly(i))
      nodal_coords(3,i) = x1(lz(i))

      basis_order(:,i)  = poly_order
      basis_x(:,1,i)    = x1
      basis_x(:,2,i)    = x1
      basis_x(:,3,i)    = x1

    end do

    basis_index(1,:)   = lx(1:ndof_cell)
    basis_index(2,:)   = ly(1:ndof_cell)
    basis_index(3,:)   = lz(1:ndof_cell)
    basis_vector(1,:)  = 1.0_r_def

  case (W1)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of Hcurl spaces
    !---------------------------------------------------------------------------

    poly_order = k+1

    do idx=1, ndof_cell
      do i=1, 3
        unit_vec(idx,i) = 0.0_r_def
      end do
    end do

    ! Compute indices of functions
    idx = 1

    ! dofs in volume
    ! u components
    do jz=2, k+1
      do jy=2, k+1
        do jx=1, k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) = tangent_to_edge(2,:)
          idx = idx + 1
        end do
      end do
    end do

    ! v components
    do jz=2, k+1
      do jy=1, k+1
        do jx=2, k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) = tangent_to_edge(1,:)
          idx = idx + 1
        end do
      end do
    end do

    ! w components
    do jz=1, k+1
      do jy=2, k+1
        do jx=2, k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) = tangent_to_edge(5,:)
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1, nfaces
      do j1=2, k+1
        do j2=1, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec(idx,:) = tangent_to_edge(edge_on_face(i,1),:)
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
      do j1=1, k+1
        do j2=2, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec(idx,:) = tangent_to_edge(edge_on_face(i,2),:)
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on edges
    do i=1, nedges
      do j1=1, k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        unit_vec(idx,:) = tangent_to_edge(i,:)
        if (i <= nedges_h )         dof_on_vert_boundary(idx,1) = 0
        if (i > nedges - nedges_h ) dof_on_vert_boundary(idx,2) = 0
        idx = idx + 1
      end do
    end do


    do i=1, ndof_cell

      nodal_coords(1,i) = abs(unit_vec(i,1))*x2(lx(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x1(lx(i))

      nodal_coords(2,i) = abs(unit_vec(i,2))*x2(ly(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x1(ly(i))

      nodal_coords(3,i) = abs(unit_vec(i,3))*x2(lz(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x1(lz(i))

      basis_order(1,i)  = poly_order - int(abs(unit_vec(i,1)))
      basis_order(2,i)  = poly_order - int(abs(unit_vec(i,2)))
      basis_order(3,i)  = poly_order - int(abs(unit_vec(i,3)))

      basis_x(:,1,i)    = abs(unit_vec(i,1))*x2(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x1(:)

      basis_x(:,2,i)    = abs(unit_vec(i,2))*x2(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x1(:)

      basis_x(:,3,i)    = abs(unit_vec(i,3))*x2(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x1(:)

      basis_vector(:,i) = unit_vec(i,:)

    end do

    basis_index(1,:) = lx(1:ndof_cell)
    basis_index(2,:) = ly(1:ndof_cell)
    basis_index(3,:) = lz(1:ndof_cell)



  case(W2)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of Hdiv spaces
    !---------------------------------------------------------------------------

    poly_order = k + 1

    do idx=1, ndof_cell
      do i=1, 3
        unit_vec(idx,i) = 0.0_r_def
      end do
    end do

    idx = 1
    ! dofs in volume
    ! u components
    do jz=1, k+1
      do jy=1, k+1
        do jx=2,k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) = normal_to_face(1,:)
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=1, k+1
      do jy=2, k+1
        do jx=1,k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) = normal_to_face(2,:)
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=2, k+1
      do jy=1, k+1
        do jx=1,k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          unit_vec(idx,:) =  normal_to_face(5,:)
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1, nfaces
      do j1=1, k+1
        do j2=1, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec(idx,:) = normal_to_face(i,:)
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    do i=1, ndof_cell

      nodal_coords(1,i) = abs(unit_vec(i,1))*x1(lx(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(lx(i))

      nodal_coords(2,i) = abs(unit_vec(i,2))*x1(ly(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(ly(i))

      nodal_coords(3,i) = abs(unit_vec(i,3))*x1(lz(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(lz(i))

      basis_order(1,i)  = poly_order - int(1 - abs(unit_vec(i,1)))
      basis_order(2,i)  = poly_order - int(1 - abs(unit_vec(i,2)))
      basis_order(3,i)  = poly_order - int(1 - abs(unit_vec(i,3)))

      basis_x(:,1,i)    = abs(unit_vec(i,1))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(:)

      basis_x(:,2,i)    = abs(unit_vec(i,2))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(:)

      basis_x(:,3,i)    = abs(unit_vec(i,3))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(:)

      basis_vector(:,i) = unit_vec(i,:)

    end do

    basis_index(1,:) = lx(1:ndof_cell)
    basis_index(2,:) = ly(1:ndof_cell)
    basis_index(3,:) = lz(1:ndof_cell)



  case(W3)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of DG spaces
    !---------------------------------------------------------------------------
    poly_order = k

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    do jz=1, k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          idx = idx + 1
        end do
      end do
    end do

    do i=1, ndof_cell
      nodal_coords(1,i) = x2(lx(i))
      nodal_coords(2,i) = x2(ly(i))
      nodal_coords(3,i) = x2(lz(i))
      basis_x(:,1,i) = x2
      basis_x(:,2,i) = x2
      basis_x(:,3,i) = x2
    end do

    basis_index(1,:)  = lx(1:ndof_cell)
    basis_index(2,:)  = ly(1:ndof_cell)
    basis_index(3,:)  = lz(1:ndof_cell)
    basis_vector(1,:) = 1.0_r_def
    basis_order(:,:)  = poly_order



  case (WTHETA)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of theta spaces
    !---------------------------------------------------------------------------
    poly_order = k + 1

    idx = 1
    ! dofs in volume - (w only)
    ! w components
    do jz=2, k+1
      do jy=1, k+1
        do jx=1, k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do


    ! dofs on faces
    do i=nfaces-1, nfaces
      do j1=1, k+1
        do j2=1, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    do i=1, ndof_cell
      nodal_coords(1,i)= x2(lx(i))
      nodal_coords(2,i)= x2(ly(i))
      nodal_coords(3,i)= x1(lz(i))

      basis_order(1,i) = poly_order - 1
      basis_order(2,i) = poly_order - 1
      basis_order(3,i) = poly_order

      basis_x(:,1,i) = x2(:)
      basis_x(:,2,i) = x2(:)
      basis_x(:,3,i) = x1(:)
    end do

    basis_index(1,:)  = lx(1:ndof_cell)
    basis_index(2,:)  = ly(1:ndof_cell)
    basis_index(3,:)  = lz(1:ndof_cell)
    basis_vector(:,:) = 1.0_r_def



  case (W2V)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of W2V space
    !---------------------------------------------------------------------------
    poly_order = k + 1

    do idx=1, ndof_cell
      do i=1, 3
        unit_vec(idx,i) = 0.0_r_def
      end do
    end do

    idx = 1
    ! dofs in volume - (w only)
    ! w components
    do jz=2, k+1
      do jy=1, k+1
        do jx=1, k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec(idx,:) = normal_to_face(5,:)
          idx = idx + 1
        end do
      end do
    end do


    ! dofs on faces
    do i=nfaces-1, nfaces
      do j1=1, k+1
        do j2=1, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec(idx,:) = normal_to_face(i,:)
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    do i=1, ndof_cell

      nodal_coords(1,i) = abs(unit_vec(i,1))*x1(lx(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(lx(i))

      nodal_coords(2,i) = abs(unit_vec(i,2))*x1(ly(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(ly(i))

      nodal_coords(3,i) = abs(unit_vec(i,3))*x1(lz(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(lz(i))

      basis_order(1,i)  = poly_order - int(1 - abs(unit_vec(i,1)))
      basis_order(2,i)  = poly_order - int(1 - abs(unit_vec(i,2)))
      basis_order(3,i)  = poly_order - int(1 - abs(unit_vec(i,3)))

      basis_x(:,1,i)    = abs(unit_vec(i,1))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(:)

      basis_x(:,2,i)    = abs(unit_vec(i,2))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(:)

      basis_x(:,3,i)    = abs(unit_vec(i,3))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(:)

      basis_vector(:,i) = unit_vec(i,:)

    end do

    basis_index(1,:) = lx(1:ndof_cell)
    basis_index(2,:) = ly(1:ndof_cell)
    basis_index(3,:) = lz(1:ndof_cell)



  case (W2H)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of W2H space
    !---------------------------------------------------------------------------
    poly_order = k + 1

    do idx=1, ndof_cell
      do i=1, 3
        unit_vec(idx,i) = 0.0_r_def
      end do
    end do

    idx = 1

    !============================================
    ! dofs in volume - (u and v only)
    !============================================
    ! u components
    do jz=1,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec(idx,:) = normal_to_face(1,:)
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=1,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec(idx,:) = normal_to_face(2,:)
          idx = idx + 1
        end do
      end do
    end do

    !============================================
    ! dofs on faces
    !============================================
    do i=1, nfaces-2
      do j1=1, k+1
        do j2=1, k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec(idx,:) = normal_to_face(i,:)
          if (i == nfaces - 1) dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do


    do i=1, ndof_cell


      nodal_coords(1,i) = abs(unit_vec(i,1))*x1(lx(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(lx(i))

      nodal_coords(2,i) = abs(unit_vec(i,2))*x1(ly(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(ly(i))

      nodal_coords(3,i) = abs(unit_vec(i,3))*x1(lz(i))                         &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(lz(i))

      basis_order(1,i)  = poly_order - int(1 - abs(unit_vec(i,1)))
      basis_order(2,i)  = poly_order - int(1 - abs(unit_vec(i,2)))
      basis_order(3,i)  = poly_order - int(1 - abs(unit_vec(i,3)))
      
      basis_x(:,1,i)    = abs(unit_vec(i,1))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,1)))*x2(:)

      basis_x(:,2,i)    = abs(unit_vec(i,2))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,2)))*x2(:)
      basis_x(:,3,i)    = abs(unit_vec(i,3))*x1(:)                             &
                        + (1.0_r_def - abs(unit_vec(i,3)))*x2(:)

      basis_vector(:,i) = unit_vec(i,:)
    end do

    basis_index(1,:) = lx(1:ndof_cell)
    basis_index(2,:) = ly(1:ndof_cell)
    basis_index(3,:) = lz(1:ndof_cell)

    case(WCHI)
    !---------------------------------------------------------------------------
    ! Section for test/trial functions of DG spaces
    !---------------------------------------------------------------------------
    poly_order = k

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    do jz=1, k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) = jx
          ly(idx) = jy
          lz(idx) = jz
          idx = idx + 1
        end do
      end do
    end do

    do i=1, ndof_cell
      nodal_coords(1,i) = x2(lx(i))
      nodal_coords(2,i) = x2(ly(i))
      nodal_coords(3,i) = x2(lz(i))
      basis_x(:,1,i) = x2
      basis_x(:,2,i) = x2
      basis_x(:,3,i) = x2
    end do

    basis_index(1,:)  = lx(1:ndof_cell)
    basis_index(2,:)  = ly(1:ndof_cell)
    basis_index(3,:)  = lz(1:ndof_cell)
    basis_vector(1,:) = 1.0_r_def
    basis_order(:,:)  = poly_order

  end select

  deallocate ( lx )
  deallocate ( ly )
  deallocate ( lz )


  ! Allocate arrays to allow on the fly evaluation of basis functions
  select case (dynamo_fs)
  case (W1, W2, W2H, W2V)
    deallocate( unit_vec )
  end select

  return
end subroutine basis_setup


  
!-----------------------------------------------------------------------------
! Setup dofmaps for the function space
!-----------------------------------------------------------------------------
!> @brief   Private routine (function_space_type construction)
!> @details Creates dofmap data for function space generation. This subroutine
!>          computes the dofmaps for the function space and stores them in a
!>          master dofmap object. The master dofmap is the same as a stencil
!>          dofmap of a single dof. This is a private
!>          routine used in the construction of function spaces.
!> @param[in] mesh Mesh to define the function space on
!> @param[in] dynamo_fs Enumeration of the function space 
!> @param[in] ncells_2d_with_ghost number of 2d cells with ghost cells
!> @param[in] ndof_vert Number of dofs on vertices
!> @param[in] ndof_edge Number of dofs on edges
!> @param[in] ndof_face Number of dofs on faces
!> @param[in] ndof_vol Number of dofs in volumes
!> @param[in] ndof_cell Number of dofs associated with a cell
!> @param[out] last_dof_owned Index of last owned dof for the partition 
!> @param[out] last_dof_annexed Index of last annexed dof for the partition
!> @param[out] last_dof_halo Index of last halo dof for the partition
!> @param[out] dofmap Array containing the dofmap indexed by cells
!> @param[out] global_dof_id Global id of of dofs
subroutine dofmap_setup( mesh, dynamo_fs, ncells_2d_with_ghost,                &
                         ndof_vert, ndof_edge, ndof_face,                      &
                         ndof_vol,  ndof_cell, last_dof_owned,                 &
                         last_dof_annexed, last_dof_halo, dofmap,              &
                         global_dof_id )

  implicit none


  type(mesh_type), intent(in), pointer :: mesh
  integer(i_def),  intent(in) :: dynamo_fs
  integer(i_def),  intent(in) :: ncells_2d_with_ghost
  integer(i_def),  intent(in) :: ndof_vert
  integer(i_def),  intent(in) :: ndof_edge
  integer(i_def),  intent(in) :: ndof_face
  integer(i_def),  intent(in) :: ndof_vol
  integer(i_def),  intent(in) :: ndof_cell
  integer(i_def), intent(out) :: last_dof_owned
  integer(i_def), intent(out) :: last_dof_annexed
  integer(i_def), intent(out) :: last_dof_halo(:)
  integer(i_def), intent(out) :: dofmap(ndof_cell,0:ncells_2d_with_ghost)

  integer(i_halo_index), intent(out) :: global_dof_id(:)


  integer(i_def) :: ncells

  ! Loop counters
  integer(i_def) :: icell, iface, iedge, ivert, idof, idepth, k

  ! Number of layers
  integer(i_def) :: nlayers

  ! Indices into the dofmap
  integer(i_def) :: id_owned, id_halo, id0, dof_idx
  integer(i_def) :: face_id, edge_id, vert_id
  integer(i_def) :: bottom_edge_id, top_edge_id, side_edge_id
  integer(i_def) :: bottom_vert_id, top_vert_id

  ! Number of entities for a single layer
  integer(i_def) :: nvert_layer, nedge_layer, nface_layer

  ! Start and end points of the cell indices to loop over
  integer(i_def) :: start,finish

  ! Entity dofmaps
  integer(i_def), allocatable :: dofmap_d0(:,:), &
                                 dofmap_d1(:,:), &
                                 dofmap_d2(:,:), &
                                 dofmap_d3(:,:)

  ! dof column heights for entities
  integer(i_def), allocatable :: dof_column_height_d0(:,:), &
                                 dof_column_height_d1(:,:), &
                                 dof_column_height_d2(:,:), &
                                 dof_column_height_d3(:,:)

  ! Cell that owns the dofs on entities
  integer(i_def), allocatable :: dof_cell_owner_d0(:,:), &
                                 dof_cell_owner_d1(:,:), &
                                 dof_cell_owner_d2(:,:), &
                                 dof_cell_owner_d3(:,:)

  ! dof column heights for whole space
  integer(i_def), allocatable :: dof_column_height(:,:)

  ! Owning cell of each entry in the dofamp
  integer(i_def), allocatable :: dof_cell_owner(:,:)

  ! Cell id in global index space
  integer(i_def) :: global_cell_id

  integer(i_def) :: dofmap_size(0:3)

  ! Number of cells in all the inner halos added together 
  integer(i_def) :: tot_num_inner

  type (select_entity_type), pointer :: select_entity => null()

  !=========================================================

  ncells = ncells_2d_with_ghost

  ! dofmaps for a 3D horizontal layer
  nlayers     =     mesh % get_nlayers()
  nvert_layer = 2 * mesh % get_nverts_2d()
  nedge_layer = 2 * mesh % get_nedges_2d() &
              +     mesh % get_nverts_2d()
  nface_layer =     mesh % get_nedges_2d() &
              + 2 * ncells

  dofmap_size(:) = 1
  dofmap_size(0) = max( dofmap_size(0), ndof_vert )
  dofmap_size(1) = max( dofmap_size(1), ndof_edge )
  dofmap_size(2) = max( dofmap_size(2), ndof_face )
  dofmap_size(3) = max( dofmap_size(3), ndof_vol  )


  allocate( dof_column_height (ndof_cell, 0:ncells))
  allocate( dof_cell_owner    (ndof_cell, 0:ncells))

  allocate( dofmap_d0             (dofmap_size(0), nvert_layer) )
  allocate( dof_column_height_d0  (dofmap_size(0), nvert_layer) )
  allocate( dof_cell_owner_d0     (dofmap_size(0), nvert_layer) )

  allocate( dofmap_d1             (dofmap_size(1), nedge_layer) )
  allocate( dof_column_height_d1  (dofmap_size(1), nedge_layer) )
  allocate( dof_cell_owner_d1     (dofmap_size(1), nedge_layer) )

  allocate( dofmap_d2             (dofmap_size(2), nface_layer) )
  allocate( dof_column_height_d2  (dofmap_size(2), nface_layer) )
  allocate( dof_cell_owner_d2     (dofmap_size(2), nface_layer) )

  allocate( dofmap_d3             (dofmap_size(3), ncells) )
  allocate( dof_column_height_d3  (dofmap_size(3), ncells) )
  allocate( dof_cell_owner_d3     (dofmap_size(3), ncells) )

  ! Initialise entity dofmaps
  dofmap_d0(:,:) = 0
  dofmap_d1(:,:) = 0
  dofmap_d2(:,:) = 0
  dofmap_d3(:,:) = 0

  dofmap_d0             (:,:) = 0
  dof_column_height_d0  (:,:) = 0
  dof_cell_owner_d0     (:,:) = 0

  dofmap_d1             (:,:) = 0
  dof_column_height_d1  (:,:) = 0
  dof_cell_owner_d1     (:,:) = 0

  dofmap_d2             (:,:) = 0
  dof_column_height_d2  (:,:) = 0
  dof_cell_owner_d2     (:,:) = 0

  dofmap_d3             (:,:) = 0
  dof_column_height_d3  (:,:) = 0
  dof_cell_owner_d3     (:,:) = 0

  ! Sum the number of cells in all the inner halos
  tot_num_inner=0
  do idepth=1,mesh%get_inner_depth()
    tot_num_inner = tot_num_inner + &
                         mesh%get_num_cells_inner(idepth)
  end do


  ! Assume we have all possible global connectivity information
  ! in practice this requires connectivity
  ! (3,2) -> faces on cells
  ! (3,1) -> edges on cells
  ! (3,0) -> vertices on cells

  id_owned = 1
  id_halo  = -1

  ! loop over 3 entities (cells) starting with core + inner halos + edge
  ! + first depth halo then proceding with further halo depths as required
  start=1
  finish=tot_num_inner + &
         mesh%get_num_cells_edge() + &
         mesh%get_num_cells_halo(1)

  select case (dynamo_fs)
  case(W0,W1,W2,W3,WCHI)
    select_entity => select_entity_all
  case(WTHETA)
    select_entity => select_entity_theta
  case(W2H)
    select_entity => select_entity_w2h
  case(W2V)
    select_entity => select_entity_w2v
  end select

  halo_loop: do idepth = 1, mesh % get_halo_depth()+1
    cell_loop: do icell = start, finish

      ! Assign dofs for connectivity (3,3) (dofs in cell)
      !---------------------------------------------------------
      if (mesh % is_cell_owned(icell)) then
        do idof=1, ndof_vol
          dofmap_d3            (idof,icell) = id_owned
          dof_column_height_d3 (idof,icell) = nlayers
          dof_cell_owner_d3    (idof,icell) = icell
          id_owned = id_owned + nlayers
        end do
      else
        do idof=1, ndof_vol
          dofmap_d3             (idof,icell) = id_halo
          dof_column_height_d3  (idof,icell) = nlayers
          dof_cell_owner_d3     (idof,icell) = icell
          id_halo = id_halo - nlayers
        end do
      end if

      ! Assign dofs for connectivity (3,2) (dofs on faces)
      !---------------------------------------------------------
      do iface=1, nfaces_h
        if (any(select_entity % faces==iface)) then
          face_id = mesh%get_face_on_cell(iface,icell)


          if (mesh%is_edge_owned(iface,icell)) then

            if ( dofmap_d2(1,face_id) == 0 ) then
              do idof=1, ndof_face
                dofmap_d2(idof,face_id) = id_owned
                dof_column_height_d2(idof,face_id) = nlayers
                dof_cell_owner_d2(idof,face_id) = &
                                     mesh%get_edge_cell_owner(iface,icell)
                id_owned = id_owned + nlayers
              end do
            end if
          else
            if ( dofmap_d2(1,face_id) == 0 ) then
              do idof=1, ndof_face
                dofmap_d2(idof,face_id) = id_halo
                dof_column_height_d2(idof,face_id) = nlayers
                dof_cell_owner_d2(idof,face_id) = &
                                     mesh%get_edge_cell_owner(iface,icell)
                id_halo = id_halo - nlayers
              end do
            end if
          end if
        end if ! select_entity
      end do

      if (mesh % is_cell_owned(icell)) then
        id0 = id_owned
        do iface=nfaces_h+1, nfaces
          if (any(select_entity % faces==iface)) then
            face_id = mesh % get_face_on_cell(iface,icell)

            if ( dofmap_d2(1,face_id) == 0 ) then
              do idof=1, ndof_face
                dofmap_d2(idof,face_id) = id_owned
                if ( iface == nfaces_h+1 ) then
                  dof_column_height_d2(idof,face_id) = nlayers + 1
                else
                  dof_column_height_d2(idof,face_id) = 0
                end if
                dof_cell_owner_d2(idof,face_id) = icell
                id_owned = id_owned + nlayers + 1
              end do
            end if

            if (iface==nfaces_h+1) then
              id_owned = id0 + 1
            else
              id_owned = id_owned - 1
            end if

          end if ! select_entity
        end do
      else
        id0 = id_halo
        do iface=nfaces_h+1, nfaces
          if (any(select_entity % faces==iface)) then
            face_id = mesh % get_face_on_cell(iface,icell)
            if ( dofmap_d2(1,face_id) == 0 ) then
              do idof=1, ndof_face
                dofmap_d2(idof,face_id) = id_halo
                if ( iface == nfaces_h+1 ) then
                  dof_column_height_d2(idof,face_id) = nlayers + 1
                else
                  dof_column_height_d2(idof,face_id) = 0
                end if
                dof_cell_owner_d2(idof,face_id) = icell
                id_halo = id_halo - nlayers - 1
              end do
            end if
            if (iface==nfaces_h+1) then
              id_halo = id0 - 1
            else
              id_halo = id_halo + 1
            end if
          end if ! select_entity
        end do
      end if ! is cell owned

      ! assign dofs for connectivity (3,1) (dofs on edges)
      do iedge=1,nedges_h
        bottom_edge_id = mesh%get_edge_on_cell(iedge,icell)
        top_edge_id    = mesh%get_edge_on_cell(iedge+nedges-nedges_h,icell)
        if (mesh%is_edge_owned(iedge,icell)) then
          if ( dofmap_d1(1,bottom_edge_id) == 0 ) then
            do idof=1,ndof_edge
              dofmap_d1(idof,bottom_edge_id)  = id_owned
              dofmap_d1(idof,top_edge_id)     = id_owned + 1
              dof_column_height_d1(idof,bottom_edge_id) = nlayers + 1
              dof_column_height_d1(idof,top_edge_id   ) = 0
              dof_cell_owner_d1(idof,bottom_edge_id) =                         &
                          mesh%get_edge_cell_owner(iedge,icell)
              dof_cell_owner_d1(idof,top_edge_id   ) =                         &
                          mesh%get_edge_cell_owner(iedge,icell)
              id_owned = id_owned + nlayers + 1
            end do
          end if
        else
          if ( dofmap_d1(1,bottom_edge_id) == 0 ) then
            do idof=1,ndof_edge
              dofmap_d1(idof,bottom_edge_id)  = id_halo
              dofmap_d1(idof,top_edge_id)     = id_halo - 1
              dof_column_height_d1(idof,bottom_edge_id) = nlayers + 1
              dof_column_height_d1(idof,top_edge_id   ) = 0
              dof_cell_owner_d1(idof,bottom_edge_id) =                         &
                          mesh%get_edge_cell_owner(iedge,icell)
              dof_cell_owner_d1(idof,top_edge_id   ) =                         &
                          mesh%get_edge_cell_owner(iedge,icell)
              id_halo = id_halo - nlayers - 1
            end do
          end if
        end if
      end do
      do iedge=nedges_h+1,nedges-nedges_h
        side_edge_id  = mesh%get_edge_on_cell(iedge,icell)
        if (mesh%is_vertex_owned(iedge-nedges_h,icell)) then
          if ( dofmap_d1(1,side_edge_id) == 0 ) then
            do idof=1,ndof_edge
              dofmap_d1(idof,side_edge_id)  = id_owned
              dof_column_height_d1(idof,side_edge_id) = nlayers
              dof_cell_owner_d1(idof,side_edge_id) =                           &
                          mesh%get_vertex_cell_owner(iedge-nedges_h,icell)
              id_owned = id_owned + nlayers
            end do
          end if
        else
          if ( dofmap_d1(1,side_edge_id) == 0 ) then
            do idof=1,ndof_edge
              dofmap_d1(idof,side_edge_id)  = id_halo
              dof_column_height_d1(idof,side_edge_id) = nlayers
              dof_cell_owner_d1(idof,side_edge_id) =                           &
                          mesh%get_vertex_cell_owner(iedge-nedges_h,icell)
              id_halo = id_halo - nlayers
            end do
          end if
        end if
      end do


      ! Assign dofs for connectivity (3,0) (dofs on verts)
      !---------------------------------------------------------
      do ivert=1, nverts_h
        bottom_vert_id  = mesh % get_vert_on_cell(ivert,icell)
        top_vert_id     = mesh % get_vert_on_cell(ivert+nverts_h,icell)

        if (mesh % is_vertex_owned(ivert,icell)) then

          if ( dofmap_d0(1,bottom_vert_id) == 0 ) then
            do idof=1, ndof_vert
              dofmap_d0(idof,bottom_vert_id)  = id_owned
              dofmap_d0(idof,top_vert_id)     = id_owned + 1
              dof_column_height_d0(idof,bottom_vert_id) = nlayers + 1
              dof_column_height_d0(idof,top_vert_id   ) = 0
              dof_cell_owner_d0(idof,bottom_vert_id) =                         &
                          mesh % get_vertex_cell_owner(ivert,icell)
              dof_cell_owner_d0(idof,top_vert_id   ) =                         &
                          mesh % get_vertex_cell_owner(ivert,icell)
              id_owned = id_owned + nlayers + 1
            end do
          end if
        else
          if ( dofmap_d0(1,bottom_vert_id) == 0 ) then
            do idof=1, ndof_vert
              dofmap_d0(idof,bottom_vert_id)  = id_halo
              dofmap_d0(idof,top_vert_id)     = id_halo - 1
              dof_column_height_d0(idof,bottom_vert_id) = nlayers + 1
              dof_column_height_d0(idof,top_vert_id   ) = 0
              dof_cell_owner_d0(idof,bottom_vert_id) =                         &
                          mesh%get_vertex_cell_owner(ivert,icell)
              dof_cell_owner_d0(idof,top_vert_id   ) =                         &
                          mesh%get_vertex_cell_owner(ivert,icell)
              id_halo = id_halo - nlayers - 1
            end do
          end if
        end if
      end do

      if(icell == tot_num_inner + mesh%get_num_cells_edge())then
        last_dof_owned = id_owned - 1
        last_dof_annexed = id_owned - id_halo - 2
      end if

    end do cell_loop

    if (idepth <= mesh%get_halo_depth())                                       &
      last_dof_halo(idepth) = id_owned - id_halo - 2

    start = finish+1
    if (idepth < mesh%get_halo_depth()) then
      finish = start + mesh % get_num_cells_halo(idepth+1)-1
    else
      finish = start + mesh % get_num_cells_ghost()-1
    end if

  end do halo_loop


  ! Copy from the dofmap_dn arrays into one dofmap array
  dof_column_height(:,:) = -999
  dof_cell_owner(:,:)    = -999
  dofmap(:,:)            = -999

  do icell=1, ncells

    dof_idx = 1

    ! dofs in volumes
    !----------------------------------------
    do idof=1, ndof_vol
      if ( dofmap_d3(idof,icell) /= 0 ) then

        if ( dofmap_d3(idof,icell) > 0 ) then
          dofmap(dof_idx,icell) = dofmap_d3(idof,icell)
        else if ( dofmap_d3(idof,icell) < 0 ) then
          dofmap(dof_idx,icell) = id_owned - (dofmap_d3(idof,icell) + 1)
        end if

        dof_column_height(dof_idx,icell) = dof_column_height_d3(idof,icell)
        dof_cell_owner(dof_idx,icell)    = dof_cell_owner_d3(idof,icell)
        dof_idx = dof_idx + 1

      end if
    end do

    ! dofs on faces
    !----------------------------------------
    do iface=1, nfaces
      face_id = mesh % get_face_on_cell(iface,icell)
      do idof=1, ndof_face
        if ( dofmap_d2(idof,face_id) /= 0 ) then
          if ( dofmap_d2(idof,face_id) > 0 ) then
            dofmap(dof_idx,icell) = dofmap_d2(idof,face_id)
          else if ( dofmap_d2(idof,face_id) < 0 ) then
            dofmap(dof_idx,icell) = id_owned - (dofmap_d2(idof,face_id) + 1)
          end if

          dof_column_height(dof_idx,icell) = dof_column_height_d2(idof,face_id)
          dof_cell_owner(dof_idx,icell)    = dof_cell_owner_d2(idof,face_id)
          dof_idx = dof_idx + 1

        end if
      end do
    end do

    ! dofs on edges
    !----------------------------------------
    do iedge=1, nedges
      edge_id = mesh % get_edge_on_cell(iedge,icell)
      do idof=1, ndof_edge
        if ( dofmap_d1(idof,edge_id) /= 0 ) then
          if ( dofmap_d1(idof,edge_id) > 0 ) then
            dofmap(dof_idx,icell) = dofmap_d1(idof,edge_id)
          else if ( dofmap_d1(idof,edge_id) < 0 ) then
            dofmap(dof_idx,icell) = id_owned - (dofmap_d1(idof,edge_id) + 1)
          end if
          dof_column_height(dof_idx,icell) = dof_column_height_d1(idof,edge_id)
          dof_cell_owner(dof_idx,icell)    = dof_cell_owner_d1(idof,edge_id)
          dof_idx = dof_idx + 1
        end if
      end do
    end do

    ! dofs on vertices
    !----------------------------------------
    do ivert=1, nverts
      vert_id = mesh % get_vert_on_cell(ivert,icell)
      do idof=1, ndof_vert
        if ( dofmap_d0(idof,vert_id) /= 0 ) then
          if ( dofmap_d0(idof,vert_id) > 0 ) then
            dofmap(dof_idx,icell) = dofmap_d0(idof,vert_id)
          else if ( dofmap_d0(idof,vert_id) < 0 ) then
            dofmap(dof_idx,icell) = id_owned - (dofmap_d0(idof,vert_id) + 1)
          end if
          dof_column_height(dof_idx,icell) = dof_column_height_d0(idof,vert_id)
          dof_cell_owner(dof_idx,icell)    = dof_cell_owner_d0(idof,vert_id)
          dof_idx = dof_idx + 1
        end if
      end do
    end do

  end do

  dofmap(:,0) = 0

  if (allocated( dofmap_d0 )) deallocate( dofmap_d0 )
  if (allocated( dofmap_d1 )) deallocate( dofmap_d1 )
  if (allocated( dofmap_d2 )) deallocate( dofmap_d2 )
  if (allocated( dofmap_d3 )) deallocate( dofmap_d3 )

  if (allocated( dof_column_height_d0 )) deallocate( dof_column_height_d0 )
  if (allocated( dof_column_height_d1 )) deallocate( dof_column_height_d1 )
  if (allocated( dof_column_height_d2 )) deallocate( dof_column_height_d2 )
  if (allocated( dof_column_height_d3 )) deallocate( dof_column_height_d3 )

  if (allocated( dof_cell_owner_d0 )) deallocate( dof_cell_owner_d0 )
  if (allocated( dof_cell_owner_d1 )) deallocate( dof_cell_owner_d1 )
  if (allocated( dof_cell_owner_d2 )) deallocate( dof_cell_owner_d2 )
  if (allocated( dof_cell_owner_d3 )) deallocate( dof_cell_owner_d3 )



  ! Calculate a globally unique id for each dof, such that each partition
  ! that needs access to that dof will calculate the same id
  global_dof_id(:) = 0_i_halo_index
  do icell=1, ncells
    global_cell_id = mesh % get_gid_from_lid(icell)
    do idof=1, ndof_cell
      if (icell == dof_cell_owner(idof,icell)) then
        do k=1, dof_column_height(idof, icell)
! The following line is very confused by the casting that is required,
!  but it is actually calculating the global id as being:
! (global_cell_id-1)*ndof_cell*(nlayers+1) + (idof-1)*(nlayers+1) + k 
          global_dof_id( dofmap(idof,icell)+k-1 ) = &
               (int(global_cell_id,i_halo_index)-1_i_halo_index)* &
      int(ndof_cell,i_halo_index)*(int(nlayers,i_halo_index)+1_i_halo_index) + &
               (int(idof,i_halo_index)-1_i_halo_index)* &
               (int(nlayers,i_halo_index)+1_i_halo_index) + int(k,i_halo_index)
        end do
      end if
    end do
  end do

  if (allocated(dof_column_height)) deallocate( dof_column_height )
  if (allocated(dof_cell_owner))    deallocate( dof_cell_owner )

  return
end subroutine dofmap_setup

!-----------------------------------------------------------------------------
! Setup the fractional levels for data on this function space
!-----------------------------------------------------------------------------
!> @brief   Creates an array of unique fractional levels
!> @details Creates an array of fractional (dof) levels for output.
!> @param[in] nlayers Number of layers
!> @param[in] fs Integer enumeration of the function space
!> @param[out] levels Array of fractional levels

subroutine levels_setup(nlayers, fs, levels)

  integer(i_def),  intent(in)          :: nlayers
  integer(i_def),  intent(in)          :: fs
  real(kind=dp_xios), dimension(:),allocatable, intent(out)   :: levels


  ! Variable to hold the number of levels we found
  integer(i_def)   :: idx
  ! working array to hold fractional levels
  real(kind=dp_xios),allocatable  :: tmp_levs(:)

  select case (fs)


  case (W0)
    ! W0 locates data on vertices

    call compute_levels(nlayers, vert_coords, select_data_entity_all%verts, tmp_levs, idx)

  case (W1)
    ! W1 locates data on edges

    call compute_levels(nlayers, edge_coords, select_data_entity_all%edges, tmp_levs, idx)

  case (W2)
    ! W2 locates data on faces 

    call compute_levels(nlayers, face_coords, select_data_entity_all%faces, tmp_levs, idx)


  case (W3)
    ! W3 locates data on cell volume

    call compute_levels(nlayers, vol_coords, select_data_entity_all%vols, tmp_levs, idx)


  case (WTHETA)
    ! WTheta locates data on selected faces
    ! (top and bottom)

    call compute_levels(nlayers, face_coords, select_data_entity_theta%faces, tmp_levs, idx)


  case (W2H)
    ! W2H locates data on selected faces
    ! (top and bottom)

    call compute_levels(nlayers, face_coords, select_data_entity_w2h%faces, tmp_levs, idx)


  case (W2V)
    ! W2V locates data on selected faces
    ! (W, S, E, N)

    call compute_levels(nlayers, face_coords, select_data_entity_w2v%faces, tmp_levs, idx)

  case default
    ! For anything else just make a dummy list

    idx = 2
    allocate(tmp_levs(1))
    tmp_levs(1) = 0

  end select

  ! Allocate the out array to be the size of the number of levels we found
  ! and copy in the data from the temp array
  allocate( levels( size(tmp_levs(1:(idx-1))) ) )
  levels=tmp_levs(1:(idx-1))

  deallocate(tmp_levs)


end subroutine levels_setup 

!-----------------------------------------------------------------------------
! Compute the fractional levels 
!-----------------------------------------------------------------------------

! A private routine to Compute fractional (dof) levels for output. Takes information
! from the data reference element combined with the number of layers to compute a unique,
! ordered list of fractional levels that fields on this function space are to be output on.


subroutine compute_levels(nlayers, coords_array, entity_array, tmp_levs, idx )

  integer(i_def),                             intent(in)    :: nlayers
  real(kind=r_def), dimension(:,:),           intent(in)    :: coords_array
  integer(i_def), dimension(:),               intent(in)    :: entity_array
  real(kind=dp_xios), dimension(:),allocatable, intent(out) :: tmp_levs
  integer(i_def),                             intent(inout) :: idx


  ! Local variables for computation
  real(kind=r_def) :: l
  integer(i_def)   :: ilayer, idof

  ! Allocate temp levels array to be the maximum possible size 
  allocate(tmp_levs(size(entity_array)*nlayers))
  tmp_levs = 999.0
  idx=1

  do ilayer=0, (nlayers - 1)
    do idof = 1, size(entity_array)
      ! Check this mesh entity is not marked as missing for this function space
      if (entity_array(idof) /= IMDI) then
        l = ilayer + coords_array(entity_array(idof),3)
        if ( .not.(any(tmp_levs == l)) ) then
          tmp_levs(idx) = l
          ! keep track of how many items we added
          idx = idx + 1
        end if
      end if
    end do
  end do

end subroutine compute_levels



end module function_space_constructor_helper_functions_mod
