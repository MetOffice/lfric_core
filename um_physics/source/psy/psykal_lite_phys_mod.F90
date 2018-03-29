!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer for physics

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_phys_mod

  use field_mod,             only : field_type, field_proxy_type 
  use mesh_mod,              only : mesh_type

  implicit none
  public

contains

  !---------------------------------------------------------------------
  !> For some unknown reason the bl kernel code accesses uninitialized values when extending
  !> into the halos.  As a result we can't use psyclone for this kernel as it will
  !> loop into the halos with the automatic openmp transformations
  !> This will be investigated in ticket #<open up a ticket for this>
  subroutine invoke_bl_kernel(theta_in_wth, rho_in_w3, rho_in_wth,           &
                              exner_in_w3, exner_in_wth, u1_in_w3, u2_in_w3, &
                              height_w3, height_wth, tstar_2d, zh_2d,        &
                              z0msea_2d, theta_inc, m_v, m_cl )

      use bl_kernel_mod, only: bl_code
      use mesh_mod, only: mesh_type
      type(field_type), intent(inout) :: theta_inc, tstar_2d, zh_2d,    &
           z0msea_2d, m_v, m_cl
      type(field_type), intent(in) :: theta_in_wth, rho_in_w3, rho_in_wth, &
           exner_in_w3, exner_in_wth, u1_in_w3, u2_in_w3, height_w3, height_wth
      integer :: cell
      integer :: ndf_wtheta, undf_wtheta, ndf_w3, undf_w3
      integer :: ndf_2d, undf_2d
      type(mesh_type), pointer :: mesh => null()
      integer :: nlayers, nlayers_2d
      type(field_proxy_type) theta_in_wth_proxy, rho_in_w3_proxy,       &
           rho_in_wth_proxy, exner_in_w3_proxy, exner_in_wth_proxy,     &
           u1_in_w3_proxy, u2_in_w3_proxy, height_w3_proxy,             &
           height_wth_proxy, tstar_2d_proxy, zh_2d_proxy, z0msea_2d_proxy, &
           theta_inc_proxy, m_v_proxy, m_cl_proxy
      integer, pointer :: map_w3(:,:) => null()
      integer, pointer :: map_wtheta(:,:) => null()
      integer, pointer :: map_2d(:,:) => null()
      !
      ! initialise field and/or operator proxies
      !
      theta_in_wth_proxy = theta_in_wth%get_proxy()
      rho_in_w3_proxy = rho_in_w3%get_proxy()
      rho_in_wth_proxy = rho_in_wth%get_proxy()
      exner_in_w3_proxy = exner_in_w3%get_proxy()
      exner_in_wth_proxy = exner_in_wth%get_proxy()
      u1_in_w3_proxy = u1_in_w3%get_proxy()
      u2_in_w3_proxy = u2_in_w3%get_proxy()
      height_w3_proxy = height_w3%get_proxy()
      height_wth_proxy = height_wth%get_proxy()
      tstar_2d_proxy = tstar_2d%get_proxy()
      zh_2d_proxy = zh_2d%get_proxy()
      z0msea_2d_proxy = z0msea_2d%get_proxy()
      theta_inc_proxy = theta_inc%get_proxy()
      m_v_proxy = m_v%get_proxy()
      m_cl_proxy = m_cl%get_proxy()      
      !
      ! initialise number of layers
      !
      nlayers    = theta_in_wth_proxy%vspace%get_nlayers()
      nlayers_2d = tstar_2d_proxy%vspace%get_nlayers()
      !
      ! create a mesh object
      !
      mesh => theta_in_wth%get_mesh()
      !
      ! look-up dofmaps for each function space
      !
      map_w3 => rho_in_w3_proxy%vspace%get_whole_dofmap()
      map_wtheta => theta_in_wth_proxy%vspace%get_whole_dofmap()
      map_2d => tstar_2d_proxy%vspace%get_whole_dofmap()
      !
      ! initialise sizes and allocate any basis arrays for wtheta
      !
      ndf_wtheta = theta_in_wth_proxy%vspace%get_ndf()
      undf_wtheta = theta_in_wth_proxy%vspace%get_undf()
      !
      ! initialise sizes and allocate any basis arrays for w3
      !
      ndf_w3 = rho_in_w3_proxy%vspace%get_ndf()
      undf_w3 = rho_in_w3_proxy%vspace%get_undf()
      !
      ! initialise sizes and allocate any basis arrays for w3 2d fields
      !
      ndf_2d = tstar_2d_proxy%vspace%get_ndf()
      undf_2d = tstar_2d_proxy%vspace%get_undf()
      !
      ! call kernels and communication routines
      !

      do cell=1,mesh%get_last_edge_cell()

        call bl_code(nlayers, nlayers_2d, theta_in_wth_proxy%data,      &
                     rho_in_w3_proxy%data, rho_in_wth_proxy%data,       &
                     exner_in_w3_proxy%data, exner_in_wth_proxy%data,   &
                     u1_in_w3_proxy%data, u2_in_w3_proxy%data,          &
                     height_w3_proxy%data, height_wth_proxy%data,       &
                     tstar_2d_proxy%data, zh_2d_proxy%data,             &
                     z0msea_2d_proxy%data, theta_inc_proxy%data,        &
                     m_v_proxy%data, m_cl_proxy%data,                   &
                     ndf_wtheta, undf_wtheta,                           &
                     map_wtheta(:,cell), ndf_w3,                        &
                     undf_w3, map_w3(:,cell),                           &
                     ndf_2d, undf_2d, map_2d(:,cell))
      end do 

    end subroutine invoke_bl_kernel

end module psykal_lite_phys_mod
