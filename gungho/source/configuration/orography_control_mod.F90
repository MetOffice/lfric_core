!-----------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------
!> @brief Sets analytic orography parameters for the model. 
!> 
!> @details Contains subroutines for reading analytic orography parameters and 
!>          constructing analytic orography types. The parameters are stored in 
!>          relevant namelist files. Each namelist sets up the selected 
!>          orography profile.
!-------------------------------------------------------------------------------
module orography_control_mod

  use constants_mod,          only : i_def, str_short, str_max_filename
  use base_mesh_config_mod,   only : geometry, &
                                     base_mesh_geometry_spherical
  use orography_config_mod,   only : profile,                  &
                                     orography_profile_schar,  & 
                                     orography_profile_agnesi
  use log_mod,                only : log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_INFO
  use analytic_orography_mod, only : orography_profile
 
  implicit none

  private 
 
  public :: set_orography_option

contains

  !=============================================================================
  !> @brief Sets analytic orography options. 
  !>
  !> @details Reads namelists with the parameters required for setting up  
  !>          analytic orography profiles and initialises corresponding types.
  !>          There are currently two profiles to choose from:
  !>          1) Schar mountain, 
  !>          2) Witch-of-Agnesi mountain.
  !>          Each profile is available for both Cartesian and spherical 
  !>          coordinates. 
  !>          The default option in no orography (flat planet surface).
  !=============================================================================
  subroutine set_orography_option()
 
    implicit none
 
    ! ----------- Deallocate abstract orography type if allocated -------------!
    if ( allocated (orography_profile) ) deallocate(orography_profile)

    ! ----------- Read orography namelist and set orography type --------------!
    select case( profile )
      ! Witch-of-Agnesi orography
      case( orography_profile_agnesi )    
        if ( geometry == base_mesh_geometry_spherical ) then
          ! Read parameters for Witch-of-Agnesi mountain in spherical 
          ! coordinates and initialise the corresponding type
           call set_orography_agnesi_spherical() 
        else 
          ! Read parameters for Witch-of-Agnesi mountain in Cartesian 
          ! coordinates and initialise the corresponding type
          call set_orography_agnesi_cartesian()
        end if 
      ! Schar orography 
      case( orography_profile_schar )   
        if ( geometry == base_mesh_geometry_spherical ) then
          ! Read parameters for Schar mountain in spherical coordinates and 
          ! initialise the corresponding type
           call set_orography_schar_spherical()
        else 
          ! Read parameters for Schar mountain in Cartesian coordinates and 
          ! initialise the corresponding type
           call set_orography_schar_cartesian()
        end if
      ! No orography (default) 
      case default   
        write(log_scratch_space,'(A,A)') &
              "set_orography_option: No analytic orography set."
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
    end select
   
    return 
  end subroutine set_orography_option

  !=============================================================================
  !> @brief Initialises analytic orography type for Witch-of-Agnesi mountain in 
  !>        spherical coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_agnesi_spherical()
   
    use agnesi_orography_spherical_mod,        only : agnesi_spherical_type
    use orography_agnesi_spherical_config_mod, only : mountain_height, &
                                                      half_width,      &
                                                      lambda_centre,   &
                                                      phi_centre,      &
                                                      lambda_focus,    &
                                                      phi_focus

    implicit none

    ! ----------- Initialise Witch-of-Agnesi spherical orography type ---------!
    allocate( orography_profile, &
              source = agnesi_spherical_type( mountain_height, &
                                              half_width,      &
                                              lambda_centre,   &
                                              phi_centre,      &
                                              lambda_focus,    &
                                              phi_focus ) )

    write(log_scratch_space,'(A,A)') "set_orography_agnesi_spherical: "// &  
          "Set analytic orography type to spherical Witch-of-Agnesi mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return 
  end subroutine set_orography_agnesi_spherical

  !=============================================================================
  !> @brief Initialises analytic orography type for Witch-of-Agnesi mountain in 
  !>        Cartesian coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_agnesi_cartesian()

    use agnesi_orography_cartesian_mod,        only : agnesi_cartesian_type
    use orography_agnesi_cartesian_config_mod, only : mountain_height, &
                                                      half_width_x,    &
                                                      half_width_y,    &
                                                      x_centre,        &
                                                      y_centre,        &
                                                      direction

    implicit none

    ! Internal variables
    integer(kind=i_def) :: direction_cart 

    ! Convert i_native configuration namelist parameter to i_def
    direction_cart = int(direction, i_def)

    ! ----------- Initialise Witch-of-Agnesi Cartesian orography type ---------!
    allocate( orography_profile,                               &
              source = agnesi_cartesian_type( mountain_height, &
                                              half_width_x,    &
                                              half_width_y,    &
                                              x_centre,        &
                                              y_centre,        &
                                              direction_cart ) )

    write(log_scratch_space,'(A,A)') "set_orography_agnesi_cartesian: "// &  
          "Set analytic orography type to Cartesian Witch-of-Agnesi mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return 
  end subroutine set_orography_agnesi_cartesian

  !=============================================================================
  !> @brief Initialises analytic orography type for Schar mountain in 
  !>        spherical coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_schar_spherical()
   
    use schar_orography_spherical_mod,        only : schar_spherical_type
    use orography_schar_spherical_config_mod, only : mountain_height, &
                                                     half_width,      &
                                                     wavelength,      &
                                                     lambda_centre,   &
                                                     phi_centre

    implicit none

    ! ----------- Initialise Schar spherical orography type -------------------!
    allocate( orography_profile,                              &
              source = schar_spherical_type( mountain_height, &
                                             half_width,      &
                                             wavelength,      &
                                             lambda_centre,   &
                                             phi_centre ) )

    write(log_scratch_space,'(A,A)') "set_orography_schar_spherical: "// &  
          "Set analytic orography type to spherical Schar mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return 
  end subroutine set_orography_schar_spherical

  !=============================================================================
  !> @brief Initialises analytic orography type for Schar mountain in 
  !>        Cartesian coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_schar_cartesian()

    use schar_orography_cartesian_mod,        only : schar_cartesian_type
    use orography_schar_cartesian_config_mod, only : mountain_height, &
                                                     half_width_x,    &
                                                     half_width_y,    &
                                                     wavelength,      &
                                                     x_centre,        &
                                                     y_centre,        &
                                                     direction

    implicit none

    ! Internal variables
    integer(kind=i_def) :: direction_cart 

    ! Convert i_native configuration namelist parameter to i_def
    direction_cart = int(direction, i_def)

    ! ----------- Initialise Schar Cartesian orography type -------------------!
    allocate( orography_profile,                              &
              source = schar_cartesian_type( mountain_height, &
                                             half_width_x,    &
                                             half_width_y,    &
                                             wavelength,      &
                                             x_centre,        &
                                             y_centre,        &
                                             direction_cart ) )

    write(log_scratch_space,'(A,A)') "set_orography_schar_cartesian: "// &  
          "Set analytic orography type to Cartesian Schar mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return 
  end subroutine set_orography_schar_cartesian

end module orography_control_mod

