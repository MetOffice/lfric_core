!-----------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------
!> @brief Module to assign the values of the surface height to model  
!> coordinates using analytic orography function.
!-------------------------------------------------------------------------------
module assign_orography_field_mod

  use constants_mod,                  only : r_def, i_def
  use base_mesh_config_mod,           only : geometry, &
                                             base_mesh_geometry_spherical
  use mesh_collection_mod,            only : mesh_collection
  use coord_transform_mod,            only : xyz2llr, llr2xyz
  use orography_helper_functions_mod, only : z2eta_linear, eta2z_linear
  use analytic_orography_mod,         only : orography_profile
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR

  implicit none

  private

  public :: assign_orography_field
  public :: assign_orography_spherical
  public :: assign_orography_cartesian

  interface

    subroutine assign_orography_interface(nlayers, ndf, undf, map,    &
                                          domain_surface, domain_top, &
                                          chi_1, chi_2, chi_3)
      import :: r_def, i_def
      implicit none
      integer(kind=i_def), intent(in) :: nlayers, undf
      integer(kind=i_def), intent(in) :: ndf
      integer(kind=i_def), intent(in) :: map(ndf)
      real(kind=r_def), intent(in)    :: domain_surface, domain_top
      real(kind=r_def), intent(inout) :: chi_1(undf), chi_2(undf), chi_3(undf)
    end subroutine assign_orography_interface

  end interface

contains

  !=============================================================================
  !> @brief Updates model vertical coordinate using selected analytic orography.  
  !>
  !> @details Model coordinate array of size 3 for the type field is passed in 
  !> to be updated. The field proxy is used to break encapsulation and access 
  !> the function space and the data attributes of the field so that its values
  !> can be updated. Model coordinates are updated by calling single column
  !> subroutines, one for spherical and the other for Cartesian domain. These
  !> routines calculate analytic orography from horizontal coordinates and then 
  !> update the vertical coordinate.  
  !>
  !> @param[in,out] chi     Model coordinate array of size 3 (x,y,z) of fields
  !> @param[in]     mesh_id Id of mesh on which this field is attached
  !============================================================================= 
  subroutine assign_orography_field(chi, mesh_id)

    use field_mod,                      only : field_type, field_proxy_type
    use mesh_mod,                       only : mesh_type
    use mesh_constructor_helper_functions_mod, &
                                        only : domain_size_type
    use orography_helper_functions_mod, only : calc_domain_size_horizontal                                             
                                             
    implicit none

    ! Arguments
    type( field_type ), intent( inout ) :: chi(3)
    integer(kind=i_def),     intent(in) :: mesh_id  
    ! Local variables
    type( field_proxy_type )     :: chi_proxy(3)
    type( mesh_type), pointer    :: mesh => null()
    type( domain_size_type )     :: domain_size 
    real(kind=r_def), pointer    :: dof_coords(:,:) => null()
    real(kind=r_def)             :: domain_top, domain_surface
    integer(kind=i_def)          :: cell
    integer(kind=i_def)          :: undf, ndf, nlayers
    integer(kind=i_def), pointer :: map(:) => null()
    ! Procedure pointer
    procedure(assign_orography_interface), pointer :: assign_orography => null()

    if ( allocated(orography_profile) ) then
      call log_event( "assign_orography_field: "// &
                      "Calculating analytic orography.", LOG_LEVEL_INFO )

      ! Get mesh object
      mesh => mesh_collection%get_mesh(mesh_id)

      ! Get domain size
      domain_size = mesh%get_domain_size()
      ! Calculate horizontal domain size from the domain_size object
      call calc_domain_size_horizontal(domain_size%minimum%x, &
                                       domain_size%maximum%x, &
                                       domain_size%minimum%y, &
                                       domain_size%maximum%y)

      ! Get physical height of flat domain surface from the domain_size object
      domain_surface = domain_size%base_height

      ! Get domain top from the mesh object and domain_surface
      domain_top = mesh%get_domain_top() + domain_surface

      ! Point to appropriate procedure to assign orography
      if ( geometry == base_mesh_geometry_spherical ) then
        assign_orography => assign_orography_spherical
      else
        assign_orography => assign_orography_cartesian
      end if

      ! Break encapsulation and get the proxy
      chi_proxy(1) = chi(1)%get_proxy()
      chi_proxy(2) = chi(2)%get_proxy()
      chi_proxy(3) = chi(3)%get_proxy()
      undf    = chi_proxy(1)%vspace%get_undf()
      ndf     = chi_proxy(1)%vspace%get_ndf( )
      nlayers = chi_proxy(1)%vspace%get_nlayers()
      ! Get DoF coordinates
      dof_coords => chi_proxy(1)%vspace%get_nodes( )

      ! Call column procedure
      do cell = 1,chi_proxy(1)%vspace%get_ncell()
        map => chi_proxy(1)%vspace%get_cell_dofmap( cell )

        call assign_orography(nlayers, &
                              ndf,               &
                              undf,              &
                              map,               &
                              domain_surface,    &
                              domain_top,        &
                              chi_proxy(1)%data, &
                              chi_proxy(2)%data, &
                              chi_proxy(3)%data )
      end do

    else

      call log_event( "assign_orography_field: No orography set "// &
                      "(flat planet surface).", LOG_LEVEL_INFO )

    end if

    return
  end subroutine assign_orography_field

  !=============================================================================
  !> @brief Updates spherical vertical coordinate for a single column using 
  !>        selected analytic orography.  
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. As model coordinates for 
  !>          spherical domain are currently (x,y,z) form they first need to be
  !>          converted to (long,lat,r) to assign orography to the model surface.
  !>          After evaluation of the new surface height chi_3 is updated using  
  !>          its nondimensional eta coordinate and then transformed back to 
  !>          (x,y,z) form.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf            Array size and loop bound for map
  !> @param[in]     undf           Column coordinates' array size and loop bound 
  !> @param[in]     map            Indirection map
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          Size undf x coord 
  !> @param[in,out] chi_2          Size undf y coord
  !> @param[in,out] chi_3          Size undf z coord
  !=============================================================================  
  subroutine assign_orography_spherical(nlayers, ndf, undf, map,    &
                                        domain_surface, domain_top, &
                                        chi_1, chi_2, chi_3)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf
    integer(kind=i_def), intent(in)    :: ndf 
    integer(kind=i_def), intent(in)    :: map(ndf) 
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf), chi_2(undf), chi_3(undf)
    ! Internal variables
    integer(kind=i_def) :: k, df, dfk
    real(kind=r_def)    :: chi_3_r
    real(kind=r_def)    :: chi_oro, surface_height
    real(kind=r_def)    :: eta
    real(kind=r_def)    :: longitude, latitude, r

    ! Calculate orography and update chi_3
    do df = 1, ndf
      do k = 0, nlayers-1
        dfk = map(df)+k 

        ! Model coordinates for spherical domain are in (x,y,z) form so they need
        ! to be converted to (long,lat,r) first   
        call xyz2llr(chi_1(dfk), chi_2(dfk), chi_3(dfk), longitude, latitude, r)

        ! Calculate surface height for each DoF using selected analytic orography
        chi_oro = orography_profile%analytic_orography(longitude, latitude)

        ! Calculate nondimensional coordinate from current height coordinate 
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(r, domain_surface, domain_top)

        ! Calculate new surface_height from flat domain_surface and orography
        surface_height = domain_surface + chi_oro

        ! Calculate new height spherical coordinate (chi_3_r) from its  
        ! nondimensional coordinate eta and surface_height
        chi_3_r = eta2z_linear(eta, surface_height, domain_top)

        ! Convert spherical coordinates back to model (x,y,z) form
        call llr2xyz(longitude, latitude, chi_3_r, &
                     chi_1(dfk), chi_2(dfk), chi_3(dfk))
      end do
    end do

    return
  end subroutine assign_orography_spherical

  !=============================================================================
  !> @brief Updates Cartesian vertical coordinate for a single column using 
  !>        selected analytic orography.  
  !>
  !> @details Calculates analytic orography from chi_1 and chi_2 horizontal
  !>          coordinates and then updates chi_3. After evaluation of the new 
  !>          surface height chi_3 is updated using its nondimensional eta 
  !>          coordinate.
  !>
  !> @param[in]     nlayers        Number of vertical layers
  !> @param[in]     ndf            Array size and loop bound for map
  !> @param[in]     undf           Column coordinates' array size and loop bound 
  !> @param[in]     map            Indirection map
  !> @param[in]     domain_surface Physical height of flat domain surface (m)
  !> @param[in]     domain_top     Physical height of domain top (m)
  !> @param[in,out] chi_1          Size undf x coord 
  !> @param[in,out] chi_2          Size undf y coord
  !> @param[in,out] chi_3          Size undf z coord
  !=============================================================================  
  subroutine assign_orography_cartesian(nlayers, ndf, undf, map,   &
                                       domain_surface, domain_top, &
                                       chi_1, chi_2, chi_3)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers, undf
    integer(kind=i_def), intent(in)    :: ndf 
    integer(kind=i_def), intent(in)    :: map(ndf) 
    real(kind=r_def),    intent(in)    :: domain_surface, domain_top
    real(kind=r_def),    intent(inout) :: chi_1(undf), chi_2(undf), chi_3(undf)
    ! Internal variables
    integer(kind=i_def) :: k, df, dfk
    real(kind=r_def)    :: chi_oro, surface_height
    real(kind=r_def)    :: eta

    ! Calculate orography and update chi_3
    do df = 1, ndf
      do k = 0, nlayers-1
        dfk = map(df)+k 

        ! Calculate surface height for each DoF using selected analytic orography
        chi_oro = orography_profile%analytic_orography(chi_1(dfk), chi_2(dfk))

        ! Calculate nondimensional coordinate from current height coordinate 
        ! (chi_3) with flat domain_surface
        eta = z2eta_linear(chi_3(dfk), domain_surface, domain_top)

        ! Calculate new surface_height from flat domain_surface and orography
        surface_height = domain_surface + chi_oro

        ! Calculate new height coordinate from its nondimensional coordinate
        ! eta and surface_height
        chi_3(dfk) = eta2z_linear(eta, surface_height, domain_top)
      end do 
    end do

    return
  end subroutine assign_orography_cartesian

end module assign_orography_field_mod