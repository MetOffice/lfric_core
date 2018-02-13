!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the fluxes through fitting a high order 1D
!>        upwind reconstruction
!> @details Compute the flux for a tracer density field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used
!>          for the polynomial is centred, with an upwind bias if an odd
!>          number of points are used, therefore for an upwind shceme an even
!>          ordered polynomial should be used (this means odd ordered for the
!>          tracer integral. In the vertical the order is reduced near the
!>          boundaries depending on the number of points available. This
!>          method is only valid for lowest order elements
!>
module sample_poly_flux_kernel_mod

  use argument_mod,  only : arg_type, func_type, mesh_data_type, &
                            GH_FIELD, GH_WRITE, GH_READ,         &
                            W2, W3, GH_BASIS, CELLS,             &
                            GH_EVALUATOR, STENCIL, CROSS,        &
                            reference_element_out_face_normal
  use constants_mod, only : r_def, i_def
  use kernel_mod,    only : kernel_type

  implicit none

  ! Precomputed operators, these are the same for all model columns
  real(kind=r_def), allocatable,    private :: coeff_matrix(:,:)
  integer(kind=i_def), allocatable, private :: dof_stencil(:,:)
  integer(kind=i_def), allocatable, private :: np_v(:)
  integer(kind=i_def),              private :: np
  real(kind=r_def),                 private :: x0
  real(kind=r_def), allocatable,    private :: x_to_p(:)
  real(kind=r_def), allocatable,    private :: coeff_matrix_v(:,:,:)

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> @brief The type declaration for the kernel. Contains the metadata needed
  !>        by the Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_poly_flux_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                    &
        arg_type(GH_FIELD,   GH_WRITE, W2),                &
        arg_type(GH_FIELD,   GH_READ,  W2),                &
        arg_type(GH_FIELD,   GH_READ,  W3, STENCIL(CROSS)) &
        /)
    type(func_type) :: meta_funcs(1) = (/ &
        func_type(W2, GH_BASIS)           &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_EVALUATOR
    type(mesh_data_type) :: meta_init(1) = (/             &
      mesh_data_type( reference_element_out_face_normal ) &
      /)
  contains
    procedure, nopass ::sample_poly_flux_code
  end type

  !-------------------------------------------------------------------------------
  ! Constructors
  !-------------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface sample_poly_flux_kernel_type
    module procedure sample_poly_flux_kernel_constructor
  end interface

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: sample_poly_flux_code, &
            sample_poly_flux_init, sample_poly_flux_final

contains

  type(sample_poly_flux_kernel_type) function sample_poly_flux_kernel_constructor() result(self)
    return
  end function sample_poly_flux_kernel_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Computes the fluxes for a tracer density
  !>
  !! @param[in]  nlayers Number of layers
  !! @param[out] flux Mass flux field to compute
  !! @param[in]  wind Wind field
  !! @param[in]  density Tracer density
  !! @param[in]  stencil_size Size of the stencil (number of cells)
  !! @param[in]  stencil_map Dofmaps for the stencil
  !! @param[in]  ndf_w2 Number of degrees of freedom per cell
  !! @param[in]  undf_w2 Number of unique degrees of freedom for the flux & wind fields
  !! @param[in]  map_w2 Dofmap for the cell at the base of the column
  !! @param[in]  basis_w2 Basis function array evaluated at w2 nodes
  !! @param[in]  ndf_w3 Number of degrees of freedom per cell
  !! @param[in]  undf_w3  Number of unique degrees of freedom for the density
  !!                      field.
  !! @param[in]  map_w3   Dofmap for the cell at the base of the column.
  !! @param[in]  out_face_normal  Vectors normal to the "out faces" of the
  !!                              reference element.
  !!
  subroutine sample_poly_flux_code( nlayers,              &
                                    flux,                 &
                                    wind,                 &
                                    density,              &
                                    stencil_size,         &
                                    stencil_map,          &
                                    ndf_w2,               &
                                    undf_w2,              &
                                    map_w2,               &
                                    basis_w2,             &
                                    ndf_w3,               &
                                    undf_w3,              &
                                    map_w3,               &
                                    out_face_normal )

    implicit none

    !Arguments
    integer(i_def), intent(in)  :: nlayers

    ! These three arguments must be defined out of declaratino order so they
    ! may be used by subsequent definitions.
    integer(i_def), intent(in)  :: ndf_w2
    integer(i_def), intent(in)  :: undf_w2
    integer(i_def), intent(in)  :: ndf_w3
    integer(i_def), intent(in)  :: undf_w3

    real(r_def),    intent(out) :: flux(undf_w2)
    real(r_def),    intent(in)  :: wind(undf_w2)
    real(r_def),    intent(in)  :: density(undf_w3)
    integer(i_def), intent(in)  :: stencil_size
    integer(i_def), intent(in)  :: stencil_map(ndf_w3,stencil_size)
    integer(i_def), intent(in)  :: map_w2(ndf_w2)
    real(r_def),    intent(in)  :: basis_w2(3,ndf_w2,ndf_w2)
    integer(i_def), intent(in)  :: map_w3(ndf_w3)
    real(r_def),    intent(in)  :: out_face_normal(:,:)

    !Internal variables
    integer(kind=i_def) :: k, df, p, i, nv
    real(kind=r_def)    :: direction, polynomial_density
    real(kind=r_def)    :: z0
    real(kind=r_def)    :: coeff(np), density_stencil(np)

    ! Horizontal terms
    do k = 0, nlayers - 1
      do df = 1,4
        ! Check if this is the upwind cell
        direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                                out_face_normal(:,df))

        if ( direction >= 0.0_r_def ) then
          do p = 1,np
            density_stencil(p) = density( stencil_map(1,dof_stencil(p,df)) + k )
          end do
          coeff(:) = matmul(coeff_matrix,density_stencil)
          polynomial_density = 0.0_r_def
          do p = 1,np
            polynomial_density = polynomial_density + coeff(p)*x_to_p(p)
          end do
          flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
        end if
      end do
    end do

    ! Vertical terms
    ! Boundary cells - apply zero flux on bottom boundary (df=5) of k=0
    !                  and top buondary (df=6) of k=nlayers-1

    ! fluxes on the bottom of each cell
    df = 5
    flux(map_w2(df) ) = 0.0_r_def
    do k = 1,nlayers - 1
      direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                              out_face_normal(:,df))
      if ( direction >= 0.0_r_def ) then
        nv = np_v(k)
        do p = 1,nv
          i = stencil_map(1,1) + k + (nv-1)/2 ! Location of first upwind point
          density_stencil(p) = density( i - (p-1)  )
        end do
      ! Use appropriate inverse matrix for the order of this point
        coeff(1:nv) = matmul(coeff_matrix_v(1:nv,1:nv,nv), &
                            density_stencil(1:nv))
        polynomial_density = 0.0_r_def
        z0 = real(1+nv/2,r_def)
        do p = 1,nv
          polynomial_density = polynomial_density + coeff(p)*z0**(p-1)
        end do
        flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
      end if
    end do

    ! fluxes on the top of each cell
    df = 6
    do k = 0,nlayers - 2
      direction = dot_product( wind(map_w2(df) + k )*basis_w2(:,df,df), &
                              out_face_normal(:,df))
      if ( direction >= 0.0_r_def ) then
        nv = np_v(k)
        do p = 1,nv
          i =  stencil_map(1,1) + k - (nv-1)/2 ! Location of first upwind point
          density_stencil(p) = density( i + (p-1) )
        end do
      ! Use appropriate inverse matrix for the order of this point
        coeff(1:nv) = matmul(coeff_matrix_v(1:nv,1:nv,nv), &
                            density_stencil(1:nv))
        polynomial_density = 0.0_r_def
        z0 = real(1+nv/2,r_def)
        do p = 1,nv
          polynomial_density = polynomial_density + coeff(p)*z0**(p-1)
        end do
        flux(map_w2(df) + k ) = wind(map_w2(df) + k)*polynomial_density
      end if
    end do
    flux(map_w2(df) + nlayers-1 ) = 0.0_r_def

  end subroutine sample_poly_flux_code

  !=============================================================================!
  !>@brief Initialise various quantities needed for sample_poly_flux_code
  !>@param[in] order Polynomial order for flux computations
  !>@param[in] nlayers Number of vertical layers
  subroutine sample_poly_flux_init(order, nlayers)

    use matrix_invert_mod, only: matrix_invert

    implicit none

    integer(kind=i_def), intent(in) :: order, nlayers

    integer(kind=i_def)             :: p, i, j, k, nup, ndown
    integer(kind=i_def)             :: nupwindcells, ndownwindcells
    real(kind=r_def), allocatable   :: inv_coeff_matrix(:,:)

    ! Order p polynomial has p+1 coefficients
    np = order + 1_i_def

    ! Compute the stencil table
    ! This give the cell id's in the stencilmap in a different order so that
    ! in each direction the first entry is the most upwind cell and the last
    ! the most downwind cell.
    ! e.g for a 1depth stencil map
    !     ---
    !     |5|
    !   -------
    !   |2|1|4|
    !   -------
    !     |3|
    !     ---
    !
    ! The stencil array will be:
    ! stencil = ( 4,1,2
    !             5,1,3
    !             2,1,4
    !             3,1,5)
    ! that is ( i+1, i, i-1
    !           j+1, j, j-1
    !           i-1, i, i+1
    !           j-1, j, j+1)
    allocate(dof_stencil(np,4))

    ! Index of first upwind cell is  (j + (nupwindcells-1)*4)
    ! where j = [2,3,4,5] for [W,S,E,N] directions
    nupwindcells = int(real(np,r_def)/2.0_r_def)
    i = 1
    do p = 1,nupwindcells
      dof_stencil(i,1) = int((4 + (nupwindcells-1)*4) - (p-1)*4,i_def)
      dof_stencil(i,2) = int((5 + (nupwindcells-1)*4) - (p-1)*4,i_def)
      dof_stencil(i,3) = int((2 + (nupwindcells-1)*4) - (p-1)*4,i_def)
      dof_stencil(i,4) = int((3 + (nupwindcells-1)*4) - (p-1)*4,i_def)
      i = i + 1
    end do
    !index of centre cell in stencil is always 1
    dof_stencil(i,:) = 1
    i = i + 1
    ! Index of first downwind cell is j
    ! where j = [5,4,3,2] for [W,S,E,N] directions
    ndownwindcells = int(real(np-1,r_def)/2.0_r_def)
    do p = 1,ndownwindcells
      dof_stencil(i,1) = int(2 + (p-1)*4,i_def)
      dof_stencil(i,2) = int(3 + (p-1)*4,i_def)
      dof_stencil(i,3) = int(4 + (p-1)*4,i_def)
      dof_stencil(i,4) = int(5 + (p-1)*4,i_def)
      i = i + 1
    end do

    ! Build Coefficient matrix of arbritray order
    allocate(coeff_matrix(np,np), inv_coeff_matrix(np,np))
    do i = 1,np
      do j = 1,np
        ! int(x^(j-1),x=i-1..i)
        inv_coeff_matrix(i,j) = (real(i,r_def)**j - real(i-1,r_def)**j)/real(j,r_def)
      end do
    end do
    call matrix_invert(inv_coeff_matrix,coeff_matrix,np)

    ! Find sampling point,
    ! first tracer point in stencil is at x = 0,
    ! flux points are at boundaries so first flux point is x = 1/2
    ! x0 then depends upon the number of upwind cells used
    x0 = real(1+int(real(np,r_def)/2.0_r_def),r_def)

    ! Compute monomials at sampling point
    allocate( x_to_p(np) )
    do p = 1,np
      x_to_p(p) = x0**(p-1)
    end do


    ! For vertical terms we may need to reduce orders near the boundaries
    ! due to lack of points so for each cell compute order and
    ! number of points (order+1) in stencil of vertical terms
    allocate( np_v(0:nlayers-1) )

    ! For a stencil symmetric about cell k (so that it is upwinded for both the
    ! dof on the top and bottom of the cell) the maximum number of points
    ! available is 2*min(nup,ndown)+1 where nup and ndown are the number of cells
    ! above and below the current cell
    do k = 1,nlayers-2
      ! # down = k
      ! # up   = (nlayers-1)-k
      ndown = k
      nup   = (nlayers-1)-k
      np_v(k) = min(np,2*min(ndown,nup)+1)
    end do
    ! For the boundary cells we use a constant approximation
    np_v(0)         = 1_i_def
    np_v(nlayers-1) = 1_i_def

    ! Compute inverses for vertical terms
    ! Each lower order matrix is just a subset of the high order matrix
    allocate(coeff_matrix_v(1:np,1:np,np) )
    do p = 1,np
      call matrix_invert(inv_coeff_matrix(1:p,1:p),coeff_matrix_v(1:p,1:p,p),p)
    end do

    deallocate( inv_coeff_matrix )

  end subroutine sample_poly_flux_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Releases claimed resources after kernel is finished with.
  !>
  subroutine sample_poly_flux_final

    implicit none

    if (allocated(coeff_matrix))    deallocate(coeff_matrix)
    if (allocated(dof_stencil))     deallocate(dof_stencil)
    if (allocated(np_v))            deallocate(np_v)
    if (allocated(x_to_p))          deallocate(x_to_p)
    if (allocated(coeff_matrix_v))  deallocate(coeff_matrix_v)

  end subroutine sample_poly_flux_final

end module sample_poly_flux_kernel_mod
