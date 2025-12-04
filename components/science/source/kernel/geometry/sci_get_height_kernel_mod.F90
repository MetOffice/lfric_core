!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Returns a height field (r or z) from the chi array.
!>
!> @details Returns a height field (r or z) from the chi array
!>
module sci_get_height_kernel_mod

  use argument_mod,              only: arg_type, func_type,       &
                                       GH_FIELD, GH_REAL,         &
                                       GH_SCALAR,                 &
                                       GH_WRITE, GH_READ, GH_INC, &
                                       ANY_DISCONTINUOUS_SPACE_1, &
                                       ANY_SPACE_9, GH_BASIS,     &
                                       CELL_COLUMN, GH_EVALUATOR
  use base_mesh_config_mod,      only: geometry, &
                                       geometry_spherical
  use constants_mod,             only: r_def, i_def
  use finite_element_config_mod, only: coord_system, &
                                       coord_system_xyz
  use kernel_mod,                only: kernel_type

  implicit none
  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: get_height_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS)                                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: get_height_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: get_height_code

contains

!> @brief Returns a height field (r or z) from the chi array
!>        Will only work at lowest order for now
!! @param[in] nlayers Number of layers
!! @param[in,out] height The height field
!! @param[in] chi_1 1st component of the coordinate
!! @param[in] chi_2 2nd component of the coordinate
!! @param[in] chi_3 3rd component of the coordinate
!! @param[in] planet_radius The planet radius
!! @param[in] ndf_x Number of degrees of freedom per cell for height
!! @param[in] undf_x Number of unique degrees of freedom for height
!! @param[in] map_x Dofmap for the cell at the base of the column for height
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] basis_chi Basis functions evaluated at nodal points for height
subroutine get_height_code(nlayers,                         &
                           height,                          &
                           chi_1, chi_2, chi_3,             &
                           planet_radius,                   &
                           ndf_x, undf_x, map_x,            &
                           ndf_chi, undf_chi, map_chi,      &
                           basis_chi                        &
                           )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def),                         intent(in) :: ndf_x, undf_x
  integer(kind=i_def),                         intent(in) :: ndf_chi, undf_chi
  real(kind=r_def),    dimension(undf_x),   intent(inout) :: height
  real(kind=r_def),    dimension(undf_chi),    intent(in) :: chi_1, chi_2, chi_3
  real(kind=r_def),                            intent(in) :: planet_radius

  integer(kind=i_def), dimension(ndf_x),           intent(in) :: map_x
  integer(kind=i_def), dimension(ndf_chi),         intent(in) :: map_chi
  real(kind=r_def),    dimension(1,ndf_chi,ndf_x), intent(in) :: basis_chi

  ! Internal variables
  integer(kind=i_def) :: df_chi, df_x, k
  real(kind=r_def)    :: coord(3), coord_radius, height_at_dof

  if ( (geometry == geometry_spherical) .and. &
       (coord_system == coord_system_xyz) ) then
    ! NB This will result in the height above
    ! the spherical representation of the planet
    ! but not necessarily the height above the bottom
    ! of the mesh
    ! This should be reviewed with ticket #562

    do k = 0, nlayers-1
      do df_x = 1, ndf_x
        coord(:) = 0.0_r_def
        do df_chi = 1, ndf_chi
          coord(1) = coord(1) + chi_1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
          coord(2) = coord(2) + chi_2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
          coord(3) = coord(3) + chi_3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        end do

        coord_radius = sqrt(coord(1)**2 + coord(2)**2 + coord(3)**2)

        height( map_x(df_x) + k ) = coord_radius - planet_radius

      end do
    end do

  else

    ! Either the domain is Cartesian or we are using a non-Cartesian spherical
    ! coordinate system. In both these cases, chi_3 will be the height.
    do k = 0, nlayers-1
      do df_x = 1, ndf_x
        height_at_dof = 0.0_r_def
        do df_chi = 1, ndf_chi
          height_at_dof = height_at_dof + &
                          chi_3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        end do

        height( map_x(df_x) + k ) = height_at_dof

      end do
    end do
  end if

end subroutine get_height_code

end module sci_get_height_kernel_mod
