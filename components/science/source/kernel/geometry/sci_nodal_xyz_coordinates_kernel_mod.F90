!-------------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @details Kernel to compute the (X,Y,Z) coordinates fields at nodal points of
!>          another function space. This computes the coordinates at the nodes
!>          in the native coordinate system before converting to (X,Y,Z).

module sci_nodal_xyz_coordinates_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_SCALAR, GH_FIELD,       &
                                    GH_WRITE, GH_READ,         &
                                    GH_REAL, GH_INTEGER,       &
                                    ANY_SPACE_9, ANY_SPACE_1,  &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_BASIS, CELL_COLUMN,     &
                                    GH_EVALUATOR
use constants_mod,           only : r_def, i_def
use sci_chi_transform_mod,   only : chi2xyz

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: nodal_xyz_coordinates_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                     &
       arg_type(GH_FIELD*3, GH_REAL, GH_WRITE, ANY_SPACE_1),              & ! nodal_x, nodal_y, nodal_z
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               & ! chi_1, chi_2, chi_3
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                         & ! geometry
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                         & ! topology
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                         & ! coord_system
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                          & ! scaled_radius
       /)
  type(func_type) :: meta_funcs(1) = (/                                  &
       func_type(ANY_SPACE_9, GH_BASIS)                                  &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: nodal_xyz_coordinates_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public nodal_xyz_coordinates_code
contains

!> @brief   Compute the coordinates fields at nodal points of another
!>          function space
!> @param[in] nlayers Number of layers
!> @param[in,out] nodal_x Nodal X coordinate
!> @param[in,out] nodal_y Nodal Y coordinate
!> @param[in,out] nodal_z Nodal Z coordinate
!> @param[in] chi1 Coordinates in the first direction
!> @param[in] chi2 Coordinates in the second direction
!> @param[in] chi3 Coordinates in the third direction
!> @param[in] panel_id A field giving the ID for mesh panels.
!> @param[in] geometry      Mesh geometry enumeration value
!> @param[in] topology      Mesh topology enumeration value
!> @param[in] coord_system  Finite-Element coord-system enumeration value
!> @param[in] scaled_radius Planet scaled radius
!> @param[in] ndf_x Number of degrees of freedom per cell for the output field
!> @param[in] undf_x Number of unique degrees of freedom for the output field
!> @param[in] map_x Dofmap for the output field
!> @param[in] ndf_chi Number of degrees of freedom per cell for the input field
!> @param[in] undf_chi Number of unique degrees of freedom for the input field
!> @param[in] map_chi Dofmap for the input field
!> @param[in] basis_chi Basis functions of the chi function space evaluated at
!>                      the nodal points of the x function space
!> @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!> @param[in] undf_pid Number of unique degrees of freedom for panel_id
!> @param[in] map_pid  Dofmap for the panel_id field
subroutine nodal_xyz_coordinates_code( nlayers,                     &
                                       nodal_x, nodal_y, nodal_z,   &
                                       chi1, chi2, chi3,            &
                                       panel_id,                    &
                                       geometry, topology,          &
                                       coord_system, scaled_radius, &
                                       ndf_x, undf_x, map_x,        &
                                       ndf_chi, undf_chi, map_chi,  &
                                       basis_chi,                   &
                                       ndf_pid, undf_pid, map_pid )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_x, ndf_chi, undf_x, undf_chi, ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_x),        intent(in)    :: map_x
  integer(kind=i_def), dimension(ndf_chi),      intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),      intent(in)    :: map_pid
  real(kind=r_def), dimension(undf_x),          intent(inout) :: nodal_x, nodal_y, nodal_z
  real(kind=r_def), dimension(undf_chi),        intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),        intent(in)    :: panel_id
  real(kind=r_def), dimension(1,ndf_chi,ndf_x), intent(in)    :: basis_chi

  integer(i_def), intent(in) :: geometry
  integer(i_def), intent(in) :: topology
  integer(i_def), intent(in) :: coord_system
  real(r_def),    intent(in) :: scaled_radius

  ! Internal variables
  integer(kind=i_def) :: df_x, df_chi, k, ipanel
  real(kind=r_def)    :: coords(3)

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df_x = 1,ndf_x
      coords(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        coords(1) = coords(1) + chi1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        coords(2) = coords(2) + chi2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        coords(3) = coords(3) + chi3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do

      call chi2xyz( coords(1), coords(2),        &
                    coords(3), ipanel,           &
                    geometry, topology,          &
                    coord_system, scaled_radius, &
                    nodal_x(map_x(df_x)+k),      &
                    nodal_y(map_x(df_x)+k),      &
                    nodal_z(map_x(df_x)+k) )
    end do
  end do

end subroutine nodal_xyz_coordinates_code

end module sci_nodal_xyz_coordinates_kernel_mod
