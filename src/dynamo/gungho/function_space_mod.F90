!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Holds information about the function space.

!> @details A container which holds type definition of the function space and 
!> has holds a number of static copies of the function spaces require by the
!> model. It provides accessor functions (getters) to various information weld
!> in the type

module function_space_mod

use constants_mod, only : r_def
use mesh_mod,  only: mesh_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  
type, public :: function_space_type
  private
  integer              :: ndf, ncell, undf, fs, nlayers, order
  integer              :: dim_space, dim_space_diff
  !> A two dimensional, allocatable array which holds the indirection map 
  !! or dofmap for the whole function space over the bottom level of the domain.
  integer, allocatable :: dofmap(:,:)
  !> A two dimensional, allocatable array of reals which holds the coordinates
  !! of the function_space degrees of freedom
  real(kind=r_def), allocatable :: nodal_coords(:,:)
  !> A two dimensional, allocatable array which holds the local orientation of
  !! vector degrees of freedom
  integer, allocatable :: orientation(:,:)
  
  !> A two dimensional, allocatable, integer array which specifies which 
  !! dofs are on vertex boundarys
  integer, allocatable :: dof_on_vert_boundary(:,:)

  !> Arrays needed for on the fly basis evaluations
  !! integer 2 dimesional, allocatable array containing the basis order
  integer, allocatable          :: basis_order(:,:)
  !! integer 2 dimesional, allocatable array containing the basis index
  integer, allocatable          :: basis_index(:,:)
  !! real, 2 dimesional, allocatable array containing the basis vector
  real(kind=r_def), allocatable :: basis_vector(:,:)
  !! real 2 dimesional, allocatable array containing the basis x
  real(kind=r_def), allocatable :: basis_x(:,:,:)

contains
  !final :: destructor

  !> @brief Returns a pointer to a function space
  !> @detail Returns a pointer to a function space. If the required function space
  !> had not yet been created, it creates one before returning the pointer to it
  !> @param[in] mesh           The parent mesh of the function space
  !> @param[in] function_space The function space id (e.g. W0)
  procedure, nopass :: get_instance

  !> @brief Gets the total unique degrees of freedom for this space,
  !! returns an integer
  !! @param[in] self The calling function space
  procedure :: get_undf

  !> @brief Returns the number of cells in the function space
  !> @param[in] self The calling function space.
  !> @return Integer, the number of cells
  procedure :: get_ncell

  !> @brief Returns the number of cells in the function space
  !> @param[in] self The calling function space.
  !> @return Integer, the number of layers
  procedure :: get_nlayers

  !> @brief Returns a pointer to the dofmap for the cell 
  !> @param[in] self The calling function_space
  !> @param[in] cell Which cell
  !> @return The pointer which points to a slice of the dofmap
  procedure :: get_cell_dofmap

  !> @brief Obtains the number of dofs per cell
  !> @param[in] self The calling functions space
  !> @return Integer, the number of dofs per cell
  procedure :: get_ndf

  !> Gets the coordinates of the function space
  !> @param[in] self The calling function_space
  !> @return A pointer to the two dimensional array of nodal_coords, (xyz,ndf)
  procedure :: get_nodes
  
  !> @brief Returns the enumerated integer for the functions_space which
  !! is this function_space
  !> @param[in] self The calling function_space_type
  !> @return fs The enumerated integer for the functions space
  procedure :: which

  !> @brief Returns a pointer to the orientation for the cell 
  !> @param[in] self The calling function_space
  !> @param[in] cell Which cell
  !> @return The pointer which points to a slice of the orientation
  procedure :: get_cell_orientation

  !> @brief Gets the flag (0) for dofs on bottom and top faces of element
  !> @param[in] self The calling function_space
  !> @return A pointer to boundary_dofs(ndf,2) the flag for bottom (:,1) and top (:,2) boundaries
  procedure :: get_boundary_dofs

  !> @brief Evaluates the basis function at a point
  !> @param[in] self The calling function space
  !> @param[in] df The dof to compute the basis function of
  !> @param[in] xi The (x,y,z) coodinates to evaluate the basis function
  procedure evaluate_basis

  !> @brief Evaluates the differential of a basis function
  !> @param[in] self The calling function space
  !> @param[in] df The dof to compute the basis function of
  !> @param[in] xi The (x,y,z) coodinates to evaluate the basis function
  procedure evaluate_diff_basis

  !> @brief Evaluates the basis function for a given quadrature
  !> @param[in] ndf integer number of dofs
  !> @param[in] qp_h integer number of quadrature points in the horizontal
  !> @param[in] qp_v integer number of quadrature points in the vertical
  !> @param[in] x_qp real two dimensional array holding the x's horizontal
  !> @param[in] z_qp real two dimensional array holding the x's vertical
  !> @param[out] basis real 3 dimensional array holding the evaluated basis 
  !! functions
  procedure compute_basis_function

  !> @brief Evaluates the differential basis function for a given quadrature
  !> @param[in] ndf integer number of dofs
  !> @param[in] qp_h integer number of quadrature points in the horizontal
  !> @param[in] qp_v integer number of quadrature points in the vertical
  !> @param[in] x_qp real two dimensional array holding the x's horizontal
  !> @param[in] z_qp real two dimensional array holding the x's vertical
  !> @param[out] dbasis real 3 dimensional array holding the evaluated basis 
  !! functions
  procedure compute_diff_basis_function

  !> @brief Gets the size of the space 
  !!(1 is scalar 3 is vector). Returns dim
  !> @param[in] self The calling get_dim_space
  !> @return dim The size of the space
  procedure get_dim_space

  !> @brief Gets the size of the differential space 
  !! (1 is scalar 3 is vector). Returns dim
  !> @param[in] self The calling get_dim_space_diff
  !> @return dim The size of the differential space
  procedure get_dim_space_diff

end type function_space_type
!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
!> integer that defines the type of function space required
integer, public, parameter      :: W0 = 100
integer, public, parameter      :: W1 = 101
integer, public, parameter      :: W2 = 102
integer, public, parameter      :: W3 = 103

!> These are static copies of all the function spaces that will be required 
type(function_space_type), target, allocatable, save :: w0_function_space
type(function_space_type), target, allocatable, save :: w1_function_space
type(function_space_type), target, allocatable, save :: w2_function_space
type(function_space_type), target, allocatable, save :: w3_function_space

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Returns a pointer to a function space
!-------------------------------------------------------------------------------
function get_instance(mesh, function_space) result(instance)
  use basis_function_mod,      only : &
              w0_order, w1_order, w2_order, w3_order, &
              w0_nodal_coords, w1_nodal_coords, w2_nodal_coords, w3_nodal_coords, &
              w0_dof_on_vert_boundary, w1_dof_on_vert_boundary, &
              w2_dof_on_vert_boundary, w3_dof_on_vert_boundary, &
              w0_basis_order, w0_basis_index, w0_basis_vector, w0_basis_x, &
              w1_basis_order, w1_basis_index, w1_basis_vector, w1_basis_x, &
              w2_basis_order, w2_basis_index, w2_basis_vector, w2_basis_x, &
              w3_basis_order, w3_basis_index, w3_basis_vector, w3_basis_x

  use dofmap_mod,              only : &
              w0_dofmap, w1_dofmap, w2_dofmap, w3_dofmap, &
              w0_orientation, w1_orientation, w2_orientation, w3_orientation

  use slush_mod, only : w_unique_dofs


  implicit none

  type (mesh_type) :: mesh
  integer :: function_space
  type(function_space_type), pointer :: instance 
  

  select case (function_space)
  case (W0)
    if(.not.allocated(w0_function_space)) then
      allocate(w0_function_space)   
      call init_function_space(self=w0_function_space, &
         order = w0_order, &
         mesh=mesh,&
         num_dofs = w_unique_dofs(1,2), &
         num_unique_dofs = w_unique_dofs(1,1) ,  &
         dim_space = 1, dim_space_diff = 3,  &
         dofmap=w0_dofmap, &
         nodal_coords=w0_nodal_coords, &
         dof_on_vert_boundary=w0_dof_on_vert_boundary, &
         orientation=w0_orientation, fs=W0, &
         basis_order=w0_basis_order, basis_index=w0_basis_index, &
         basis_vector=w0_basis_vector, basis_x=w0_basis_x) 
    end if
    instance => w0_function_space
  case (W1)
    if(.not.allocated(w1_function_space)) then
      allocate(w1_function_space) 
      call init_function_space(self=w1_function_space, &
         order = w1_order, &
         mesh=mesh,&
         num_dofs = w_unique_dofs(2,2), &
         num_unique_dofs = w_unique_dofs(2,1) ,  &
         dim_space = 3, dim_space_diff = 3,  &
         dofmap=w1_dofmap, &
         nodal_coords=w1_nodal_coords, &
         dof_on_vert_boundary=w1_dof_on_vert_boundary, &
         orientation=w1_orientation, fs=W1, &
         basis_order=w1_basis_order, basis_index=w1_basis_index, &
         basis_vector=w1_basis_vector, basis_x=w1_basis_x )
    end if
    instance => w1_function_space
  case (W2)
    if(.not.allocated(w2_function_space)) then 
      allocate(w2_function_space)
      call init_function_space(self=w2_function_space, &
         order = w2_order, &
         mesh=mesh,&
         num_dofs = w_unique_dofs(3,2), &
         num_unique_dofs = w_unique_dofs(3,1) ,  &
         dim_space = 3, dim_space_diff = 1,  &
         dofmap=w2_dofmap, &
         nodal_coords=w2_nodal_coords, &
         dof_on_vert_boundary=w2_dof_on_vert_boundary, &
         orientation=w2_orientation, fs=W2, &
         basis_order=w2_basis_order, basis_index=w2_basis_index, &
         basis_vector=w2_basis_vector, basis_x=w2_basis_x )
    end if
    instance => w2_function_space
  case (W3)
    if(.not.allocated(w3_function_space)) then
      allocate(w3_function_space)
      call init_function_space(self=w3_function_space, &
         order = w3_order, &
         mesh=mesh,&
         num_dofs = w_unique_dofs(4,2), &
         num_unique_dofs = w_unique_dofs(4,1) ,  &
         dim_space = 1, dim_space_diff = 1,  &
         dofmap=w3_dofmap, &
         nodal_coords=w3_nodal_coords, &
         dof_on_vert_boundary=w3_dof_on_vert_boundary, &
         orientation=w3_orientation, fs=W3, &
         basis_order=w3_basis_order, basis_index=w3_basis_index, &
         basis_vector=w3_basis_vector, basis_x=w3_basis_x )
    end if
    instance => w3_function_space
  case default
    !not a recognised function space - return a null pointer
    instance => null()
  end select

  return
end function get_instance

!-----------------------------------------------------------------------------
! Initialises a function space
!-----------------------------------------------------------------------------
!> @brief Initialises a function space.
!> @param[in] mesh object to assign function space to
!> @param[in] num_dofs
!> @param[in] num_unique_dofs
!> @param[in] dim_space The dimension of this function space
!> @param[in] dim_space_diff The dimension of the differentiated function space
!> @param[in] ngp_h The number of guassian quadrature points in the horizonal
!> @param[in] ngp_v The number of guassian quadrature points in the vertical
!-----------------------------------------------------------------------------
subroutine init_function_space(self, &
                               order, &
                               mesh,&
                               num_dofs, &
                               num_unique_dofs,  &
                               dim_space, dim_space_diff,  &
                               dofmap, &
                               nodal_coords, &
                               dof_on_vert_boundary, &
                               orientation ,fs, &
                               basis_order, basis_index, &
                               basis_vector, basis_x)
  implicit none

  class(function_space_type)  :: self
  integer, intent(in)         :: order
  type(mesh_type), intent(in) :: mesh
  integer, intent(in)         :: num_dofs, num_unique_dofs
  integer, intent(in)         :: dim_space, dim_space_diff

! The following four arrays have intent inout because the move_allocs in the
! code need access to the arrays to free them in their original locations
  integer,          intent(inout), allocatable  :: dofmap(:,:)
  real(kind=r_def), intent(inout), allocatable  :: nodal_coords(:,:)
  integer,          intent(inout), allocatable  :: dof_on_vert_boundary(:,:)
  integer,          intent(inout), allocatable  :: orientation(:,:)
  integer,          intent(in)                  :: fs
  integer,          intent(inout), allocatable  :: basis_order(:,:),  basis_index(:,:)
  real(kind=r_def), intent(inout), allocatable  :: basis_vector(:,:), basis_x(:,:,:)
  self%order           =  order
  self%ncell           =  mesh%get_ncells_2d()
  self%nlayers         =  mesh%get_nlayers()
  self%ndf             =  num_dofs
  self%undf            =  num_unique_dofs
  self%dim_space       =  dim_space
  self%dim_space_diff  =  dim_space_diff
  call move_alloc(dofmap, self%dofmap)
  call move_alloc(nodal_coords , self%nodal_coords) 
  call move_alloc(dof_on_vert_boundary , self%dof_on_vert_boundary) 
  call move_alloc(orientation , self%orientation) 
  self%fs              = fs
  call move_alloc(basis_order,self%basis_order)
  call move_alloc(basis_index,self%basis_index)
  call move_alloc(basis_vector,self%basis_vector)
  call move_alloc(basis_x,self%basis_x)

  return
end subroutine init_function_space

!-----------------------------------------------------------------------------
! Gets total unique dofs for this space
!-----------------------------------------------------------------------------
integer function get_undf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_undf=self%undf

  return
end function get_undf

!-----------------------------------------------------------------------------
! Gets the number of cells for this function space
!-----------------------------------------------------------------------------
integer function get_ncell(self)
  implicit none
  class(function_space_type), intent(in) :: self
  get_ncell=self%ncell

  return
end function get_ncell

!-----------------------------------------------------------------------------
! Gets the number of layers for this functions space 
!-----------------------------------------------------------------------------
integer function get_nlayers(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_nlayers=self%nlayers

  return
end function get_nlayers

!-----------------------------------------------------------------------------
! Gets the number of dofs for a single cell 
!-----------------------------------------------------------------------------
integer function get_ndf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_ndf=self%ndf

  return
end function get_ndf

!-----------------------------------------------------------------------------
! Gets the dofmap for a single cell
!-----------------------------------------------------------------------------
function get_cell_dofmap(self,cell) result(map)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer                               :: map(:) 

  map => self%dofmap(:,cell)
  return
end function get_cell_dofmap

! ----------------------------------------------------------------
! Gets the nodal coordinates of the function_space
! ----------------------------------------------------------------
function get_nodes(self) result(nodal_coords)
  implicit none
  class(function_space_type), target, intent(in)  :: self
  real(kind=r_def),              pointer          :: nodal_coords(:,:) 
  
  nodal_coords => self%nodal_coords
  
  return
end function get_nodes

!-----------------------------------------------------------------------------
! Gets the orientation for a single cell
!-----------------------------------------------------------------------------
function get_cell_orientation(self,cell) result(cell_orientation)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer                               :: cell_orientation(:)

  cell_orientation => self%orientation(cell,:)
  return
end function get_cell_orientation

!-----------------------------------------------------------------------------
! Gets a flag for dofs on vertical boundaries
!-----------------------------------------------------------------------------
function get_boundary_dofs(self) result(boundary_dofs)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer, pointer                               :: boundary_dofs(:,:) 
  
  boundary_dofs => self%dof_on_vert_boundary(:,:) 
  return
end function get_boundary_dofs

!-----------------------------------------------------------------------------
! Gets enumerated integer for the function space
!-----------------------------------------------------------------------------
function which(self) result(fs)
  implicit none
  class(function_space_type),  intent(in) :: self
  integer :: fs
  
  fs = self%fs
  return
end function which

!-----------------------------------------------------------------------------
! Gets the size of the function space
!-----------------------------------------------------------------------------
function get_dim_space(self) result(dim)
  implicit none
  class(function_space_type), intent(in) :: self
  integer :: dim
  dim = self%dim_space
  return
end function get_dim_space

!-----------------------------------------------------------------------------
! Gets the size of the diferential function space
!-----------------------------------------------------------------------------
function get_dim_space_diff(self) result(dim)
  implicit none
  class(function_space_type), intent(in) :: self
  integer :: dim
  dim = self%dim_space_diff
  return
end function get_dim_space_diff

!-----------------------------------------------------------------------------
! Evaluates a basis function at a point
!-----------------------------------------------------------------------------
 function evaluate_basis(self, df, xi) result(p)
  use polynomial_mod, only: poly1d

  class(function_space_type), intent(in)  :: self
  integer,                    intent(in)  :: df
  real(kind=r_def),           intent(in)  :: xi(3)
  real(kind=r_def)                        :: p(self%dim_space)

  p(:) = poly1d(self%basis_order(1,df), xi(1), self%basis_x(:,1,df), self%basis_index(1,df)) &
        *poly1d(self%basis_order(2,df), xi(2), self%basis_x(:,2,df), self%basis_index(2,df)) &
        *poly1d(self%basis_order(3,df), xi(3), self%basis_x(:,3,df), self%basis_index(3,df)) &
        *self%basis_vector(:,df)

end function evaluate_basis

!-----------------------------------------------------------------------------
! Evaluates the differential of a basis function at a point
!-----------------------------------------------------------------------------
pure function evaluate_diff_basis(self, df, xi) result(dp)
  use polynomial_mod, only: poly1d, poly1d_deriv
  class(function_space_type), intent(in)  :: self
  integer,                    intent(in)  :: df
  real(kind=r_def),           intent(in)  :: xi(3)  
  real(kind=r_def)                        :: dp(self%dim_space_diff)
  real(kind=r_def)                        :: dpdx(3)

  dpdx(1) = poly1d_deriv(self%basis_order(1,df), xi(1), self%basis_x(:,1,df), self%basis_index(1,df)) &
           *poly1d      (self%basis_order(2,df), xi(2), self%basis_x(:,2,df), self%basis_index(2,df)) &
           *poly1d      (self%basis_order(3,df), xi(3), self%basis_x(:,3,df), self%basis_index(3,df))
  dpdx(2) = poly1d      (self%basis_order(1,df), xi(1), self%basis_x(:,1,df), self%basis_index(1,df)) &
           *poly1d_deriv(self%basis_order(2,df), xi(2), self%basis_x(:,2,df), self%basis_index(2,df)) &
           *poly1d      (self%basis_order(3,df), xi(3), self%basis_x(:,3,df), self%basis_index(3,df))
  dpdx(3) = poly1d      (self%basis_order(1,df), xi(1), self%basis_x(:,1,df), self%basis_index(1,df)) &
           *poly1d      (self%basis_order(2,df), xi(2), self%basis_x(:,2,df), self%basis_index(2,df)) &
           *poly1d_deriv(self%basis_order(3,df), xi(3), self%basis_x(:,3,df), self%basis_index(3,df))


  if ( self%dim_space == 1 .and. self%dim_space_diff == 3 ) then
! grad(p)
    dp(1) = dpdx(1)
    dp(2) = dpdx(2) 
    dp(3) = dpdx(3)
  elseif ( self%dim_space == 3 .and. self%dim_space_diff == 3 ) then
! curl(p)
    dp(1) = dpdx(2)*self%basis_vector(3,df) - dpdx(3)*self%basis_vector(2,df)
    dp(2) = dpdx(3)*self%basis_vector(1,df) - dpdx(1)*self%basis_vector(3,df)
    dp(3) = dpdx(1)*self%basis_vector(2,df) - dpdx(2)*self%basis_vector(1,df)
  elseif ( self%dim_space == 3 .and. self%dim_space_diff == 1 ) then
! div(p)
    dp(1) = dpdx(1)*self%basis_vector(1,df) + dpdx(2)*self%basis_vector(2,df) + dpdx(3)*self%basis_vector(3,df)
  else
    dp(:) = 0.0_r_def
  end if

end function evaluate_diff_basis

!-----------------------------------------------------------------------------
! Evaluates the basis function for a given quadrature
!-----------------------------------------------------------------------------
subroutine compute_basis_function(self, &
     basis, ndf, qp_h, qp_v, x_qp, z_qp)
  implicit none
  class(function_space_type), intent(in)  :: self
  integer,                                             intent(in)  :: ndf
  integer,                                             intent(in)  :: qp_h
  integer,                                             intent(in)  :: qp_v
  real(kind=r_def), dimension(qp_h,2),                 intent(in)  :: x_qp
  real(kind=r_def), dimension(qp_v),                   intent(in)  :: z_qp
  real(kind=r_def), dimension(self%dim_space,ndf,qp_h,qp_v), intent(out) :: basis

  ! local variables - loop counters
  integer :: df
  real(kind=r_def) :: xyz(3)
  integer :: qp1
  integer :: qp2

  do qp2 = 1, qp_v
     xyz(3)=z_qp(qp2)
     do qp1 = 1, qp_h
        xyz(1) = x_qp(qp1,1)
        xyz(2) = x_qp(qp1,2)
        do df = 1, ndf
           basis(:,df,qp1,qp2) = self%evaluate_basis(df,xyz)
        end do
     end do
  end do
  
end subroutine compute_basis_function

!-----------------------------------------------------------------------------
! Evaluates the differential basis function for a given quadrature
!-----------------------------------------------------------------------------
subroutine compute_diff_basis_function(self, &
     dbasis, ndf, qp_h, qp_v, x_qp, z_qp )
  implicit none
  class(function_space_type), intent(in)  :: self
  integer,                                             intent(in)  :: ndf
  integer,                                             intent(in)  :: qp_h
  integer,                                             intent(in)  :: qp_v
  real(kind=r_def), dimension(qp_h,2),                 intent(in)  :: x_qp
  real(kind=r_def), dimension(qp_v),                   intent(in)  :: z_qp
  real(kind=r_def), dimension(self%dim_space_diff,ndf,qp_h,qp_v), intent(out) :: dbasis

! local variables - loop counters
  integer :: df
  real(kind=r_def) :: xyz(3)
  integer :: qp1
  integer :: qp2
  
  do qp2 = 1, qp_v
     xyz(3)=z_qp(qp2)
     do qp1 = 1, qp_h
        xyz(1) = x_qp(qp1,1)
        xyz(2) = x_qp(qp1,2)
        do df = 1, ndf
           dbasis(:,df,qp1,qp2) = self%evaluate_diff_basis(df,xyz)
        end do
     end do
  end do

end subroutine compute_diff_basis_function

!-----------------------------------------------------------------------------
! Gets order for this space
!-----------------------------------------------------------------------------
!> @brief Gets the polynomial order for this space, returns an integer
!> @param[in] self the calling function space
!-----------------------------------------------------------------------------
integer function get_order(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_order=self%order

  return
end function get_order
end module function_space_mod
