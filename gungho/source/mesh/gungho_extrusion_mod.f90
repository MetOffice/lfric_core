!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!>
!> @brief Sets up vertical extrusion of mesh
!> @details This code contains two functions, create_extrusion() and
!>          create_shifted_extrusion(). The function create_extrusion() generates
!>          the standard vertical extrusion with different options, including
!>          uniform, quadratic, geometric and dcmip spacing.
!>          create_shifted_extrusion() creates a vertical extrusion with the same
!>          options as create_extrusion() but the top and bottom layers are
!>          half the normal cell height.
!>          There are technical infrastructure limitations which mean two different
!>          functions have been used to create the normal vertical mesh and the shifted
!>          vertical mesh. Tickets #1645 and #1659 deal with this issue of multiple
!>          instances of the mesh with different vertical extrusion.

module gungho_extrusion_mod

  use base_mesh_config_mod, only : geometry,          &
                                   key_from_geometry, &
                                   geometry_planar,   &
                                   geometry_spherical
  use constants_mod,        only : r_def, i_def
  use extrusion_mod,        only : extrusion_type,                             &
                                   uniform_extrusion_type,                     &
                                   quadratic_extrusion_type,                   &
                                   geometric_extrusion_type,                   &
                                   dcmip_extrusion_type,                       &
                                   shifted_uniform_extrusion_type,             &
                                   shifted_quadratic_extrusion_type,           &
                                   shifted_geometric_extrusion_type,           &
                                   shifted_dcmip_extrusion_type,               &
                                   um_L38_29t_9s_40km_extrusion_type
  use extrusion_config_mod, only : method,                    &
                                   key_from_method,           &
                                   method_uniform,            &
                                   method_quadratic,          &
                                   method_geometric,          &
                                   method_dcmip,              &
                                   method_um_L38_29t_9s_40km, &
                                   domain_top,                &
                                   number_of_layers
  use log_mod,              only : log_event,       &
                                   log_level_error, &
                                   log_scratch_space
  use planet_config_mod,    only : scaled_radius

  implicit none

  private
  public create_extrusion, create_shifted_extrusion

  character(*), parameter :: module_name = 'gungho_extrusion_mod'

contains

  !> @brief Creates vertical mesh extrusion
  !> @details Creates vertical mesh with nlayers.
  !> @return new     Extrusion class
  function create_extrusion() result(new)

    implicit none

    class(extrusion_type), allocatable :: new

    real(r_def) :: atmosphere_bottom

    if (allocated(new)) deallocate(new)

    select case (geometry)
      case (geometry_planar)
        atmosphere_bottom = 0.0_r_def
      case (geometry_spherical)
        atmosphere_bottom = scaled_radius
      case default
        write( log_scratch_space,                      &
               '(A, ": Unrecognised geometry: ", A)' ) &
             module_name, key_from_geometry( geometry )
        call log_event( log_scratch_space, log_level_error )
    end select

    select case (method)
      case (method_uniform)
        allocate( new, source=uniform_extrusion_type( atmosphere_bottom, &
                                                      domain_top,        &
                                                      number_of_layers ) )
      case (method_um_L38_29t_9s_40km)
        allocate( new, source=um_L38_29t_9s_40km_extrusion_type( atmosphere_bottom, &
                                                                 domain_top,        & 
                                                                 number_of_layers ) )
      case (method_quadratic)
        allocate( new, source=quadratic_extrusion_type( atmosphere_bottom, &
                                                        domain_top,        &
                                                        number_of_layers ) )
      case (method_geometric)
        allocate( new, source=geometric_extrusion_type( atmosphere_bottom, &
                                                        domain_top,        &
                                                        number_of_layers ) )
      case (method_dcmip)
        allocate( new, source=dcmip_extrusion_type( atmosphere_bottom, &
                                                    domain_top,        &
                                                    number_of_layers ) )
      case default
        write( log_scratch_space,                         &
               '(A, ": Unrecognised extrusion method: ", A)' ) &
             module_name, key_from_method( method )
        call log_event( log_scratch_space, log_level_error )
    end select

  end function create_extrusion

  !> @brief Creates vertical mesh extrusion for vertically shifted mesh.
  !> @details Creates vertically shifted mesh with nlayers+1 with the top and 
  !>          bottom levels having half the cell height of the normal extrusion.
  !> @return new     Extrusion class
  function create_shifted_extrusion() result(new)

    implicit none

    class(extrusion_type), allocatable :: new

    real(r_def) :: atmosphere_bottom
    integer(i_def) :: nlayers_plus_one

    nlayers_plus_one = number_of_layers + 1

    if (allocated(new)) deallocate(new)

    select case (geometry)
      case (geometry_planar)
        atmosphere_bottom = 0.0_r_def
      case (geometry_spherical)
        atmosphere_bottom = scaled_radius
      case default
        write( log_scratch_space,                      &
               '(A, ": Unrecognised geometry: ", A)' ) &
             module_name, key_from_geometry( geometry )
        call log_event( log_scratch_space, log_level_error )
    end select

    select case (method)
      case (method_uniform)
        allocate( new, source=shifted_uniform_extrusion_type( atmosphere_bottom, &
                                                              domain_top,        &
                                                              nlayers_plus_one ) )
      case (method_quadratic)
        allocate( new, source=shifted_quadratic_extrusion_type( atmosphere_bottom, &
                                                                domain_top,        &
                                                                nlayers_plus_one ) )
      case (method_geometric)
        allocate( new, source=shifted_geometric_extrusion_type( atmosphere_bottom, &
                                                                domain_top,        &
                                                                nlayers_plus_one ) )
      case (method_dcmip)
        allocate( new, source=shifted_dcmip_extrusion_type( atmosphere_bottom, &
                                                            domain_top,        &
                                                            nlayers_plus_one ) )
      case default
        write( log_scratch_space,                         &
               '(A, ": Unrecognised shifted extrusion method: ", A)' ) &
             module_name, key_from_method( method )
        call log_event( log_scratch_space, log_level_error )
    end select

  end function create_shifted_extrusion

end module gungho_extrusion_mod
