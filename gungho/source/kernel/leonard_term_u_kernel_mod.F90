!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates Leonard term vertical fluxes and increments for any
!>        field defined at W3-points.
!>
module leonard_term_u_kernel_mod

  use argument_mod,          only : arg_type,                     &
                                    GH_FIELD, GH_SCALAR, GH_REAL, &
                                    GH_READ, GH_WRITE, GH_INC,    &
                                    CELL_COLUMN, STENCIL, REGION
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta, W2, W1
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: leonard_term_u_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                   &
         arg_type(GH_FIELD,  GH_REAL, GH_INC,   W2),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2, STENCIL(REGION)),      &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta, STENCIL(REGION)),  &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W1),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                       &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: leonard_term_u_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: leonard_term_u_code

contains

!> @brief Calculates increment due to vertical Leonard term flux
!>        for variables held in wt space
!! @param[in] nlayers  Number  of layers in the mesh
!! @param[in,out] u_inc  Leonard term increment on w2
!! @param[in] u_n Wind on w2
!! @param[in] map_w2_stencil_size  Number of cells in the stencil at the base
!!                                 of the column for w2
!! @param[in] map_w2_stencil  Array holding the dofmap for the stencil at the
!!                            base of the column for w2
!! @param[in] velocity_w2v  Velocity normal to the top face of the cell
!! @param[in] map_wt_stencil_size  Number of cells in the stencil at the base
!!                                 of the column for Wtheta
!! @param[in] map_wt_stencil  Array holding the dofmap for the stencil at the
!!                            base of the column for Wtheta
!! @param[in] vel_w2v_inc  Leonard term increment of velocity_w2v
!! @param[in] dtrdz_fd2  Array of dt/(r*dz) at FD2 points
!! @param[in] height_w1  Height of w1 space levels above the surface
!! @param[in] height_w2  Height of w2 space levels above the surface
!! @param[in] wetrho_in_w2  Density on w2 space levels
!! @param[in] planet_radius  The planet radius
!! @param[in] leonard_kl  The user-specified Leonard term parameter
!! @param[in] dt  The model timestep length
!! @param[in] ndf_w2  Number of degrees of freedom per cell for w2 space
!! @param[in] undf_w2  Number of unique degrees of freedom for w2 space
!! @param[in] map_w2  Cell dofmap for w2 space
!! @param[in] ndf_wt  Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] map_wt  Cell dofmap for theta space
!! @param[in] ndf_w1  Number of degrees of freedom per cell for w1 space
!! @param[in] undf_w1  Number of unique degrees of freedom for w1 space
!! @param[in] map_w1  Cell dofmap for w1 space
subroutine leonard_term_u_code( nlayers,                                &
                                 u_inc,                                 &
                                 u_n,                                   &
                                 map_w2_stencil_size, map_w2_stencil,   &
                                 velocity_w2v,                          &
                                 map_wt_stencil_size, map_wt_stencil,   &
                                 vel_w2v_inc,                           &
                                 dtrdz_fd2,                             &
                                 height_w1,                             &
                                 height_w2,                             &
                                 wetrho_in_w2,                          &
                                 planet_radius,                         &
                                 leonard_kl,                            &
                                 dt,                                    &
                                 ndf_w2, undf_w2, map_w2,               &
                                 ndf_wt, undf_wt, map_wt,               &
                                 ndf_w1, undf_w1, map_w1                &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: map_w2_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_stencil_size), intent(in)  :: map_w2_stencil
  integer(kind=i_def), intent(in) :: map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_w1),  intent(in)  :: map_w1
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: velocity_w2v
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: vel_w2v_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: dtrdz_fd2
  real(kind=r_def), dimension(undf_w1),  intent(in)    :: height_w1
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: height_w2
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: wetrho_in_w2
  real(kind=r_def),                      intent(in)    :: planet_radius
  real(kind=r_def),                      intent(in)    :: leonard_kl
  real(kind=r_def),                      intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k, km, df
  integer(kind=i_def) :: df2, df2p1, df2p2, df2m1, df2m2, dfp2x2

  ! interpolation weights
  real(kind=r_def)    :: weight_pl, weight_min
  ! density at FD1 points
  real(kind=r_def)    :: rho_fd1
  ! Leonard term parameter at FD1 points
  real(kind=r_def), dimension(1:nlayers-1,4) :: kl_fd1
  ! Leonard term vertical flux at FD1 points
  real(kind=r_def), dimension(0:nlayers-1,4) :: flux
  ! density * r^2 at FD1
  real(kind=r_def) , dimension(1:nlayers-1,4):: rho_rsq_fd1
  ! timestep / ( dz * rho * r^2 ) at FD2
  real(kind=r_def) , dimension(0:nlayers-2,4):: dtrdzrho_fd2

  ! If the full stencil isn't available, we must be at the domain edge.
  ! The increment is already 0, so we just exit the routine.
  if (map_w2_stencil_size < 9_i_def) then
    return
  end if

  ! Loop over horizontal faces of the cell
  do df = 1,4

    if (u_inc(map_w2(df)+1) == 0.0_r_def) then

      ! Calculate indices for finite differences between faces of cells
      df2 = 2_i_def * df
      df2m1 = df2 - 1_i_def
      df2m2 = df2 - 2_i_def
      df2p1 = df2 + 1_i_def
      df2p2 = df2 + 2_i_def
      dfp2x2 = 2_i_def * (df + 2_i_def)

      ! If index is less than 2 add 8
      if (df2m1 < 2_i_def) then
        df2m1 = df2m1 + 8_i_def
      end if
      if (df2m2 < 2_i_def) then
        df2m2 = df2m2 + 8_i_def
      end if
      ! If index is greater than 9 subtract 8
      if (df2p2 > 9_i_def) then
        df2p2 = df2p2 - 8_i_def
      end if
      if (dfp2x2 > 9_i_def) then
        dfp2x2 = dfp2x2 - 8_i_def
      end if

      ! Calculate kl at FD1 points,
      ! accounting for stability limit
      do k = 1, nlayers - 1
        kl_fd1(k,df) = MIN( leonard_kl,                                &
                    6.0_r_def * ( height_w2(map_w2(df) + k) -          &
                                  height_w2(map_w2(df) + k-1) )        &
                    / (dt * MAX(                                       &
                    ! Difference normal to face
                    ! (point lies exactly between these):
                    ABS( velocity_w2v(map_wt_stencil(1,df2) + k) -     &
                         velocity_w2v(map_wt_stencil(1,1) + k) ),      &
                    ! Terms from gradient parallel to face...
                    ! from difference to the left (when looking at face):
                    ABS( velocity_w2v(map_wt_stencil(1,1) + k) -       &
                         velocity_w2v(map_wt_stencil(1,df2m2) + k) ),  &
                    ABS( velocity_w2v(map_wt_stencil(1,df2) + k) -     &
                         velocity_w2v(map_wt_stencil(1,df2m1) + k) ),  &
                    ! from difference to the right (when looking at face):
                    ABS( velocity_w2v(map_wt_stencil(1,df2p2) + k) -   &
                         velocity_w2v(map_wt_stencil(1,1) + k) ),      &
                    ABS( velocity_w2v(map_wt_stencil(1,df2p1) + k) -   &
                         velocity_w2v(map_wt_stencil(1,df2) + k) ),    &
                    EPSILON( leonard_kl )                              &
                    ) ) )
      end do

      ! Calculate rho * r^2 at FD1
      do k = 1, nlayers - 1
        km = k - 1

        ! Vertical interpolation weights:
        weight_pl = (height_w1(map_w1(df) + k) - height_w2(map_w2(df) + km)) /  &
                    (height_w2(map_w2(df) + k) - height_w2(map_w2(df) + km))
        weight_min = (height_w2(map_w2(df) + k) - height_w1(map_w1(df) + k)) /  &
                     (height_w2(map_w2(df) + k) - height_w2(map_w2(df) + km))

        ! Vertical interpolation of wetrho_in_w2 from w2 to FD1 points
        rho_fd1 = ( weight_min * wetrho_in_w2(map_w2(df) + km) ) +  &
                  ( weight_pl * wetrho_in_w2(map_w2(df) + k) )

        ! Calculate rho * r^2 at FD1
        rho_rsq_fd1(k,df) = rho_fd1 *                                       &
                           ( height_w1(map_w1(df) + k) + planet_radius ) *  &
                           ( height_w1(map_w1(df) + k) + planet_radius )

      end do

      ! Calculate timestep / ( dz * rho * r^2 ) at w2
      do k = 0, nlayers - 2
        dtrdzrho_fd2(k,df) = dtrdz_fd2(map_w2(df) + k) / &
                             wetrho_in_w2(map_w2(df) + k)
      end do

      ! Calculate vertical flux at FD1 points:
      ! Leonard term sub-grid vertical flux is proportional to the
      ! product of the grid-scale horizontal differences in u and w.

      ! Set flux and increment to zero at k=0
      k = 0
      flux(k,df) = 0.0_r_def

      do k = 1, nlayers - 1
        km = k - 1
        flux(k,df) = ( kl_fd1(k,df) / 12.0_r_def )                    &
              ! 8 terms contribute to each direction, so scale by 1/8
              * ( 1.0_r_def / 8.0_r_def ) * (                         &
              ! Terms from gradient normal to face...
              ! ( point is half-way between 2 w-points,
              !   so only one dw_x contributes ):
              2.0_r_def * ( ( u_n(map_w2_stencil(df,dfp2x2) + km) -   &
                              u_n(map_w2_stencil(df,df2) + km) )      &
                          + ( u_n(map_w2_stencil(df,dfp2x2) + k) -    &
                              u_n(map_w2_stencil(df,df2) + k) ) )     &
              * ( velocity_w2v(map_wt_stencil(1,1) + k) -             &
                  velocity_w2v(map_wt_stencil(1,df2) + k) )           &
              ! Terms from gradient parallel to face...
              ! from difference to the left (when looking at face):
              + ( ( u_n(map_w2_stencil(df,1) + km) -                  &
                    u_n(map_w2_stencil(df,df2m2) + km) )              &
                + ( u_n(map_w2_stencil(df,1) + k) -                   &
                    u_n(map_w2_stencil(df,df2m2) + k) ) )             &
              * ( ( velocity_w2v(map_wt_stencil(1,1) + k) -           &
                    velocity_w2v(map_wt_stencil(1,df2m2) + k) )       &
                + ( velocity_w2v(map_wt_stencil(1,df2) + k) -         &
                    velocity_w2v(map_wt_stencil(1,df2m1) + k) ) )     &
              ! from difference to the right (when looking at face):
              + ( ( u_n(map_w2_stencil(df,df2p2) + km) -              &
                    u_n(map_w2_stencil(df,1) + km) )                  &
                + ( u_n(map_w2_stencil(df,df2p2) + k) -               &
                    u_n(map_w2_stencil(df,1) + k) ) )                 &
              * ( ( velocity_w2v(map_wt_stencil(1,df2p2) + k) -       &
                    velocity_w2v(map_wt_stencil(1,1) + k) )           &
                + ( velocity_w2v(map_wt_stencil(1,df2p1) + k) -       &
                    velocity_w2v(map_wt_stencil(1,df2) + k) ) )       &
              )                                                       &
              * rho_rsq_fd1(k,df)
              ! Flux now includes density * r^2 factor
      end do

      ! Difference the flux in the vertical to get increment on w2
      do k = 0, nlayers - 2
        u_inc(map_w2(df) + k) = -(flux(k+1,df) - flux(k,df))           &
                                * dtrdzrho_fd2(k,df)
      end do

    end if

  end do

  ! Add vertical velocity increment to vertical DoFs (5/6) of u_inc
  ! to give the total vector increment in a single field
  do k = 0, nlayers - 1
    u_inc(map_w2(5) + k) = vel_w2v_inc(map_wt(1) + k)
  end do

end subroutine leonard_term_u_code

end module leonard_term_u_kernel_mod
