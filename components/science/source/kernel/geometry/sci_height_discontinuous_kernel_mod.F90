!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Computes the height above the Earth's mean radius at the degrees of
!!        freedom for a discontinuous function space.
!> @details Uses the chi coordinate fields to determine the height of nodes
!!          above the Earth's mean radius, for discontinuous function spaces.

module sci_height_discontinuous_kernel_mod

  use argument_mod,              only: arg_type, func_type,                    &
                                       GH_FIELD, GH_SCALAR,                    &
                                       GH_REAL, GH_INTEGER,                    &
                                       GH_READ, GH_WRITE,                      &
                                       ANY_DISCONTINUOUS_SPACE_1, ANY_SPACE_9, &
                                       CELL_COLUMN, GH_BASIS, GH_EVALUATOR
  use base_mesh_config_mod,      only: geometry_spherical
  use constants_mod,             only: r_def, i_def, l_def
  use finite_element_config_mod, only: coord_system_xyz
  use kernel_mod,                only: kernel_type

  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  type, public, extends(kernel_type) :: height_discontinuous_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9),               &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              &
    /)
    type(func_type) :: meta_funcs(1) = (/                                      &
        func_type(ANY_SPACE_9, GH_BASIS)                                       &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: height_discontinuous_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: height_discontinuous_code

contains

!> @brief Computes the height above the Earth's mean radius at the degrees of
!!        freedom for a discontinuous function space.
!> @param[in]     nlayers        Number of layers in the mesh
!> @param[in,out] height         The height field to compute
!> @param[in]     chi_1          1st component of the coordinate fields
!> @param[in]     chi_2          2nd component of the coordinate fields
!> @param[in]     chi_3          3rd component of the coordinate fields
!> @param[in]     geometry       The geometry of the domain
!> @param[in]     coord_system   The coordinate system of the domain
!> @param[in]     planet_radius  The planet radius
!> @param[in]     ndf_h          Num DoFs per cell for height field
!> @param[in]     undf_h         Num DoFs in this partition for height field
!> @param[in]     map_h          DoF index map for height field
!> @param[in]     ndf_chi        Num DoFs per cell for chi fields
!> @param[in]     undf_chi       Num DoFs in this partition for chi fields
!> @param[in]     map_chi        DoF index map for chi fields
!> @param[in]     basis_chi      Chi basis functions evaluated at height DoFs
subroutine height_discontinuous_code(                                          &
    nlayers,                                                                   &
    height,                                                                    &
    chi_1, chi_2, chi_3,                                                       &
    geometry, coord_system, planet_radius,                                     &
    ndf_h, undf_h, map_h,                                                      &
    ndf_chi, undf_chi, map_chi,                                                &
    basis_chi                                                                  &
)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_h, undf_h
  integer(kind=i_def), intent(in)    :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in)    :: geometry, coord_system
  integer(kind=i_def), intent(in)    :: map_h(ndf_h)
  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  real(kind=r_def),    intent(in)    :: planet_radius
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,ndf_h)
  real(kind=r_def),    intent(inout) :: height(undf_h)
  real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_3(undf_chi)

  ! Internal variables
  logical(kind=l_def) :: is_spherical_xyz
  integer(kind=i_def) :: df_chi, df_h
  integer(kind=i_def) :: h_b_idx, h_t_idx, chi_b_idx, chi_t_idx
  real(kind=r_def)    :: coord_1(nlayers), coord_2(nlayers), coord_3(nlayers)
  real(kind=r_def)    :: coord_radius(nlayers), basis_val

  ! Two cases: if the geometry is spherical with geocentric Cartesian coords,
  ! all three chi components are needed to compute the height. Otherwise, only
  ! the last chi component is needed, as this is already the vertical coordinate

  is_spherical_xyz = (                                                         &
    geometry == geometry_spherical .and. coord_system == coord_system_xyz      &
  )

  ! -------------------------------------------------------------------------- !
  ! Cartesian system and spherical
  ! -------------------------------------------------------------------------- !
  if (is_spherical_xyz) then

    do df_h = 1, ndf_h
      ! Indices for bottom and top height DoFs in this column
      h_b_idx = map_h(df_h)
      h_t_idx = map_h(df_h) + nlayers - 1

      ! Determine coordinates at this DoF using chi basis functions
      coord_1(:) = 0.0_r_def
      coord_2(:) = 0.0_r_def
      coord_3(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        ! Indices for bottom and top chi DoFs in this column
        chi_b_idx = map_chi(df_chi)
        chi_t_idx = map_chi(df_chi) + nlayers - 1
        basis_val = basis_chi(1, df_chi, df_h)

        coord_1(:) = coord_1(:) + chi_1(chi_b_idx:chi_t_idx)*basis_val
        coord_2(:) = coord_2(:) + chi_2(chi_b_idx:chi_t_idx)*basis_val
        coord_3(:) = coord_3(:) + chi_3(chi_b_idx:chi_t_idx)*basis_val
      end do

      ! Convert to radial coordinate
      coord_radius(:) = sqrt(coord_1(:)**2 + coord_2(:)**2 + coord_3(:)**2)

      ! Increment height field at this DoF
      height(h_b_idx:h_t_idx) = coord_radius(:) - planet_radius
    end do

  ! -------------------------------------------------------------------------- !
  ! Native coordinates, or planar domain
  ! -------------------------------------------------------------------------- !
  else
    ! Third coordinate already gives the height above the Earth's mean radius
    do df_h = 1, ndf_h
      ! Indices for bottom and top height DoFs in this column
      h_b_idx = map_h(df_h)
      h_t_idx = map_h(df_h) + nlayers - 1
      height(h_b_idx:h_t_idx) = 0.0_r_def

      do df_chi = 1, ndf_chi
        ! Indices for bottom and top chi DoFs in this column
        chi_b_idx = map_chi(df_chi)
        chi_t_idx = map_chi(df_chi) + nlayers - 1
        basis_val = basis_chi(1, df_chi, df_h)

        ! Increment height field at this DoF
        height(h_b_idx:h_t_idx) = (                                            &
            height(h_b_idx:h_t_idx) + chi_3(chi_b_idx:chi_t_idx)*basis_val     &
        )
      end do
    end do
  end if

end subroutine height_discontinuous_code

end module sci_height_discontinuous_kernel_mod
