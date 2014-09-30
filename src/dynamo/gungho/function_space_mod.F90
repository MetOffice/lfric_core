!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A type which holds information about the function space.

!> @details A container which holds type definition of the function space and 
!> has holds a number of static copies of the function spaces require by the
!> model. It provides accessor functions (getters) to various information weld
!> in the type

module function_space_mod

use constants_mod, only : r_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  
type, public :: function_space_type
  private
  integer              :: ndf, ncell, undf, ngp_h, ngp_v, fs, nlayers
  integer              :: dim_space, dim_space_diff
  !> A two dimensional, allocatable array which holds the indirection map 
  !! or dofmap for the whole function space over the bottom level of the domain.
  integer, allocatable :: dofmap(:,:)
  !> 4-dim allocatable array of reals which hold the values of the basis function
  real(kind=r_def), allocatable :: basis(:,:,:,:)
  !> 4-dim allocatable array of reals which hold the values of the basis function
  !! for the differential  functions space
  real(kind=r_def), allocatable :: diff_basis(:,:,:,:)

  !> A two dimensional, allocatable array of reals which holds the coordinates
  !! of the function_space degrees of freedom
  real(kind=r_def), allocatable :: nodal_coords(:,:)
  !> A two dimensional, allocatable array which holds the local orientation of
  !! vector degrees of freedom
  integer, allocatable :: orientation(:,:)
  !> 3-dim allocatable array of reals which hold the mass matrix
  real(kind=r_def), allocatable :: mass_matrix(:,:,:)  
  
  integer, allocatable :: dof_on_vert_boundary(:,:)

  !> Arrays needed for on the fly basis evaluations
  integer, allocatable          :: basis_order(:,:)
  integer, allocatable          :: basis_index(:,:)
  real(kind=r_def), allocatable :: basis_vector(:,:)
  real(kind=r_def), allocatable :: basis_x(:,:,:)

contains
  !final :: destructor

  !> Function returns a pointer to a function space. If the required function
  !> space had not yet been created, it creates one before returning the pointer
  !> to it
  procedure, nopass :: get_instance
!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
  procedure :: get_undf

!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
  procedure :: get_ncell

!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of layers
  procedure :: get_nlayers

!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
  procedure :: get_cell_dofmap

!> Function which obtains the number of dofs per cell
!! @param[in] self The calling functions space
!! return an integer, the number of dofs per cell
  procedure :: get_ndf

!>  Accessor procedure for the basis function
!! @param[in] self the calling function space
!! @return A pointer to the array to hold the values of the basis function 
  procedure :: get_basis

!> Accessor procedure for the basis function for the
!! differential function space
!! @param[in] self the calling function space
!! @return A pointer to the real array to hold the values of the
!!  basis function 
  procedure :: get_diff_basis

!> Accessor function to get the coordinates of the function space
!! @return A pointer to the two dimensional array, (xyz,ndf)
  procedure :: get_nodes
  
  !> function returns the enumerated integer for the functions_space which
  !! is this function_space
  procedure :: which

!> Subroutine Returns a pointer to the orientation for the cell 
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the orientation
  procedure :: get_cell_orientation

!> Accessor function to get the flag (0) for dofs on bottom and top faces of element
!! @return A pointer to the two dimensional array, (ndf,2)
  procedure :: get_boundary_dofs

!> Accessor function to evaluate a basis function
  procedure evaluate_basis
!> Accessor function to evaluate the differential of a basis function
  procedure evaluate_diff_basis
  
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

function get_instance(function_space) result(instance)
  use basis_function_mod,      only : &
              w0_basis, w1_basis, w2_basis, w3_basis, &
              w0_diff_basis, w1_diff_basis, w2_diff_basis, w3_diff_basis, &
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

  use gaussian_quadrature_mod, only : ngp_h, ngp_v
  use mesh_mod,                only : num_cells, w_unique_dofs, num_layers

  implicit none

  integer :: function_space
  type(function_space_type), pointer :: instance

  select case (function_space)
  case (W0)
    if(.not.allocated(w0_function_space)) then
      allocate(w0_function_space)   
      call init_function_space(self=w0_function_space, &
         num_cells = num_cells , num_layers = num_layers, &
         num_dofs = w_unique_dofs(1,2), &
         num_unique_dofs = w_unique_dofs(1,1) ,  &
         dim_space = 1, dim_space_diff = 3,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=w0_dofmap, &
         basis=w0_basis, diff_basis=w0_diff_basis, &
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
         num_cells = num_cells ,num_layers = num_layers, &
         num_dofs = w_unique_dofs(2,2), &
         num_unique_dofs = w_unique_dofs(2,1) ,  &
         dim_space = 3, dim_space_diff = 3,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=w1_dofmap, &
         basis=w1_basis, diff_basis=w1_diff_basis, &
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
         num_cells = num_cells ,num_layers = num_layers, &
         num_dofs = w_unique_dofs(3,2), &
         num_unique_dofs = w_unique_dofs(3,1) ,  &
         dim_space = 3, dim_space_diff = 1,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=w2_dofmap, &
         basis=w2_basis, diff_basis=w2_diff_basis, &
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
         num_cells = num_cells ,num_layers = num_layers, &
         num_dofs = w_unique_dofs(4,2), &
         num_unique_dofs = w_unique_dofs(4,1) ,  &
         dim_space = 1, dim_space_diff = 1,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=w3_dofmap, &
         basis=w3_basis, diff_basis=w3_diff_basis, &
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


!> Subroutine initialises a function space.
!! @param[in] num_cells
!! @param[in] num_dofs
!! @param[in] num_unique_dofs
!! @param[in] dim_space The dimension of this function space
!! @param[in] dim_space_diff The dimension of the differentiated function space
!! @param[in] ngp_h The number of guassian quadrature points in the horizonal
!! @param[in] ngp_v The number of guassian quadrature points in the vertical
subroutine init_function_space(self, &
                               num_cells,num_layers, &
                               num_dofs, &
                               num_unique_dofs,  &
                               dim_space, dim_space_diff,  &
                               ngp_h,ngp_v, &
                               dofmap, &
                               basis, diff_basis, &
                               nodal_coords, &
                               dof_on_vert_boundary, &
                               orientation ,fs, &
                               basis_order, basis_index, &
                               basis_vector, basis_x)
  implicit none

  class(function_space_type) :: self
  integer, intent(in) :: num_cells, num_layers, num_dofs, num_unique_dofs
  integer, intent(in) :: dim_space, dim_space_diff
  integer, intent(in) :: ngp_h,ngp_v
! The following four arrays have intent inout because the move_allocs in the
! code need access to the arrays to free them in their original locations
  integer,          intent(inout), allocatable  :: dofmap(:,:)
  real(kind=r_def), intent(inout), allocatable  :: basis(:,:,:,:)
  real(kind=r_def), intent(inout), allocatable  :: diff_basis(:,:,:,:)
  real(kind=r_def), intent(inout), allocatable  :: nodal_coords(:,:)
  integer,          intent(inout), allocatable  :: dof_on_vert_boundary(:,:)
  integer,          intent(inout), allocatable  :: orientation(:,:)
  integer,          intent(in)                  :: fs
  integer,          intent(inout), allocatable  :: basis_order(:,:),  basis_index(:,:)
  real(kind=r_def), intent(inout), allocatable  :: basis_vector(:,:), basis_x(:,:,:)

  self%ncell           =  num_cells
  self%nlayers         =  num_layers
  self%ndf             =  num_dofs
  self%undf            =  num_unique_dofs
  self%dim_space       =  dim_space
  self%dim_space_diff  =  dim_space_diff
  self%ngp_h           =  ngp_h
  self%ngp_v           =  ngp_v  
  call move_alloc(dofmap, self%dofmap)
  call move_alloc(basis , self%basis)
  call move_alloc(diff_basis , self%diff_basis)
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
! Get total unique dofs for this space
!-----------------------------------------------------------------------------

!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
integer function get_undf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_undf=self%undf

  return
end function get_undf

!-----------------------------------------------------------------------------
! Get the number of cells for this function space
!-----------------------------------------------------------------------------
!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
integer function get_ncell(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_ncell=self%ncell

  return
end function get_ncell

!-----------------------------------------------------------------------------
! Get the number of layers for this functions space 
!-----------------------------------------------------------------------------
integer function get_nlayers(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_nlayers=self%nlayers

  return
end function get_nlayers

!-----------------------------------------------------------------------------
! Get the number of dofs for a single cell 
!-----------------------------------------------------------------------------
integer function get_ndf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_ndf=self%ndf

  return
end function get_ndf

!-----------------------------------------------------------------------------
! Get the dofmap for a single cell
!-----------------------------------------------------------------------------
!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
function get_cell_dofmap(self,cell) result(map)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer                               :: map(:)

  map => self%dofmap(cell,:)
  return
end function get_cell_dofmap

!-----------------------------------------------------------------------------
! Get the basis function
!-----------------------------------------------------------------------------
!> Subroutine to return the basis functions for this space
!! @param[in] self The calling function_space
!! @return The pointer which points to the basis functions
function get_basis(self)  result(basis)
  implicit none
  class(function_space_type), target, intent(in)  :: self  
  real(kind=r_def),              pointer          :: basis(:,:,:,:)

  basis => self%basis

  return
end function get_basis

!-----------------------------------------------------------------------------
! Get the differential of the basis function
!-----------------------------------------------------------------------------
!> Subroutine to return the differential basis functions for this space
!! @param[in] self The calling function_space
!! @return The pointer which points to the differenrtial basis functions
function get_diff_basis(self) result(diff_basis)
  implicit none
  class(function_space_type), target, intent(in)  :: self  
  real(kind=r_def),              pointer          :: diff_basis(:,:,:,:)

  diff_basis => self%diff_basis

  return
end function get_diff_basis

! ----------------------------------------------------------------
! Get the nodal coordinates of the function_space
! ----------------------------------------------------------------
!> Subroutine to return the dof is location
!! @param[in] self The calling function_space
!! @return The pointer which points to the nodal_coords
function get_nodes(self) result(nodal_coords)
  implicit none
  class(function_space_type), target, intent(in)  :: self
  real(kind=r_def),              pointer          :: nodal_coords(:,:)
  
  nodal_coords => self%nodal_coords
  
  return
end function get_nodes

!-----------------------------------------------------------------------------
! Get the orientation for a single cell
!-----------------------------------------------------------------------------
!> Subroutine Returns a pointer to the orientation for the cell 
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the orientation
function get_cell_orientation(self,cell) result(cell_orientation)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer                               :: cell_orientation(:)

  cell_orientation => self%orientation(cell,:)
  return
end function get_cell_orientation

!-----------------------------------------------------------------------------
! Get a flag for dofs on vertical boundaries
!-----------------------------------------------------------------------------
!> Subroutine returns a pointer to flag for dofs on vertical boundaries 
!! @param[in] self The calling function_space
!! @param[in] boundary_dofs(ndf,2) the flag for bottom (:,1) and top (:,2) boundaries
function get_boundary_dofs(self) result(boundary_dofs)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer, pointer                               :: boundary_dofs(:,:)
  
  boundary_dofs => self%dof_on_vert_boundary(:,:) 
  return
end function get_boundary_dofs

function which(self) result(fs)
  implicit none
  class(function_space_type),  intent(in) :: self
  integer :: fs
  
  fs = self%fs
  return
end function which

!-----------------------------------------------------------------------------
! Evaluate a basis function at a point
!-----------------------------------------------------------------------------
!> Function to evaluate the basis function at a point
!> @param[in] self The calling function space
!> @param[in] df The dof to compute the basis function of
!> @param[in] xi The (x,y,z) coodinates to evaluate the basis function
pure function evaluate_basis(self, df, xi) result(p)
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
! Evaluate the differential of a basis function at a point
!-----------------------------------------------------------------------------
!> Function to evaluate the differential of a basis function at a point
!> @param[in] self The calling function space
!> @param[in] df The dof to compute the basis function of
!> @param[in] xi The (x,y,z) coodinates to evaluate the basis function
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


end module function_space_mod
