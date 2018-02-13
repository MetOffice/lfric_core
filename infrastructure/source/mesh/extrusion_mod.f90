!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Provides extrusion methods for converting a 2D mesh to a unitless
!>        3D mesh.
!>
!> An abstract extrusion_type holds the data and provides getters to access
!> it. Concrete classes derived from it then implement the specific extrusion
!> routine.
!>
!> The result of this design is that once an extrusion object is created it
!> may be passed around as its abstract base class. Thus the ultimate point
!> of use need know nothing about which extrusion is being used.
!>
module extrusion_mod

  use constants_mod,         only : i_def, r_def
  use global_mesh_mod,       only : global_mesh_type
  use log_mod,               only : log_scratch_space, log_event, &
                                    log_level_error
  use reference_element_mod, only : reference_element_type

  implicit none

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief All extrusion implementations inherit from this class.
  !>
  type, public, abstract :: extrusion_type

    private

    real(r_def)    :: atmosphere_bottom
    real(r_def)    :: atmosphere_top
    integer(i_def) :: number_of_layers

  contains

    private

    procedure, public :: get_atmosphere_bottom
    procedure, public :: get_atmosphere_top
    procedure, public :: get_number_of_layers
    procedure, public :: get_reference_element
    procedure(extrude_method), public, deferred :: extrude

    procedure :: extrusion_constructor

  end type extrusion_type

  interface
    subroutine extrude_method( this, eta )
      import extrusion_type, r_def
      class(extrusion_type), intent(in)  :: this
      real(r_def),           intent(out) :: eta(0:this%number_of_layers)
    end subroutine extrude_method
  end interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with equal distribution of layers.
  !>
  type, public, extends(extrusion_type) :: uniform_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => uniform_extrude
  end type uniform_extrusion_type

  interface uniform_extrusion_type
    module procedure uniform_extrusion_constructor
  end interface uniform_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with a @f$\left(\frac{layer}{n_{layers}}\right)^2@f$
  !>        distribution of layers.
  !>
  type, public, extends(extrusion_type) :: quadratic_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => quadratic_extrude
  end type quadratic_extrusion_type

  interface quadratic_extrusion_type
    module procedure quadratic_extrusion_constructor
  end interface quadratic_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with "geometric" layers.
  !>
  type, public, extends(extrusion_type) :: geometric_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => geometric_extrude
  end type geometric_extrusion_type

  interface geometric_extrusion_type
    module procedure geometric_extrusion_constructor
  end interface geometric_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes using DCMIP scheme.
  !>
  type, public, extends(extrusion_type) :: dcmip_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => dcmip_extrude
  end type dcmip_extrusion_type

  interface dcmip_extrusion_type
    module procedure dcmip_extrusion_constructor
  end interface dcmip_extrusion_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a uniform_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function uniform_extrusion_constructor( atmosphere_bottom, &
                                          atmosphere_top,    &
                                          number_of_layers ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers

    type(uniform_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers )

  end function uniform_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh to give constant delta between layers.
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine uniform_extrude( this, eta )

    implicit none

    class(uniform_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:this%number_of_layers)

    integer(i_def) :: k

    do k = 0, this%number_of_layers
      eta(k) = real(k,r_def)/real(this%number_of_layers,r_def)
    end do

  end subroutine uniform_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a quadratic_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !>
  !> @return New quadratic_extrusion_type object.
  !>
  function quadratic_extrusion_constructor( atmosphere_bottom, &
                                            atmosphere_top,    &
                                            number_of_layers ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers

    type(quadratic_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers )

  end function quadratic_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh to give layer boundaries
  !>        @f$\frac{l}{n_{layers}}^2@f$.
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine quadratic_extrude( this, eta )

    implicit none

    class(quadratic_extrusion_type), intent(in)  :: this
    real(r_def),                     intent(out) :: eta(0:this%number_of_layers)

    integer(i_def) :: k

    do k = 0, this%number_of_layers
      eta(k) = ( real(k,r_def)/real(this%number_of_layers,r_def) )**2_i_def
    end do

  end subroutine quadratic_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a geometric_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !>
  !> @return New geometric_extrusion_type object.
  !>
  function geometric_extrusion_constructor( atmosphere_bottom, &
                                            atmosphere_top,    &
                                            number_of_layers ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers

    type(geometric_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers )

  end function geometric_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh to give a John Thuburn ENDGame non-staggered grid.
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine geometric_extrude( this, eta )

    implicit none

    class(geometric_extrusion_type), intent(in)  :: this
    real(r_def),                     intent(out) :: eta(0:this%number_of_layers)

    real(r_def), parameter :: stretching_factor = 1.03_r_def

    integer(i_def) :: k
    real(r_def)    :: delta_eta

    delta_eta = (stretching_factor - 1.0_r_def) &
                / (stretching_factor**(this%number_of_layers) - 1.0_r_def)

    eta(0) = 0.0_r_def
    do k = 1, this%number_of_layers
      eta(k) = eta(k-1) + delta_eta
      delta_eta = delta_eta*stretching_factor
    end do

  end subroutine geometric_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a dcmip_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !>
  !> @return New dcmip_extrusion_type object.
  !>
  function dcmip_extrusion_constructor( atmosphere_bottom, &
                                        atmosphere_top,    &
                                        number_of_layers ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers

    type(dcmip_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers )

  end function dcmip_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh using the DCMIP scheme.
  !>
  !> For more information see DCMIP-TestCaseDocument_v1.7.pdf,
  !> Appendix F.2. - Eq. 229.
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine dcmip_extrude( this, eta )

    implicit none

    class(dcmip_extrusion_type), intent(in)  :: this
    real(r_def),                 intent(out) :: eta(0:this%number_of_layers)

    real(r_def), parameter :: phi_flatten = 15.0_r_def

    integer(i_def) :: k
    real(r_def)    :: eta_uni

    do k = 0, this%number_of_layers
      eta_uni = real(k,r_def)/real(this%number_of_layers,r_def)
      eta(k) = ( sqrt(phi_flatten*(eta_uni**2_i_def) + 1.0_r_def) &
                      - 1.0_r_def ) / &
                    ( sqrt(phi_flatten + 1.0_r_def) - 1.0_r_def )
    end do

  end subroutine dcmip_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Initialises the extrusion base class.
  !>
  !> This method should be called from child method constructors in order to
  !> populate the parent fields.
  !>
  !> @param [in] atmosphere_bottom Bottom of the atmosphere (planet surface)
  !>                               in meters.
  !> @param [in] atmosphere_top Top of the atmosphere in meters.
  !> @param [in] number_of_layers Number of layers to split atmosphere into.
  !>
  subroutine extrusion_constructor( this,              &
                                    atmosphere_bottom, &
                                    atmosphere_top,    &
                                    number_of_layers )

    implicit none

    class(extrusion_type), intent(inout) :: this
    real(r_def),           intent(in) :: atmosphere_bottom
    real(r_def),           intent(in) :: atmosphere_top
    integer(i_def),        intent(in) :: number_of_layers

    this%atmosphere_bottom = atmosphere_bottom
    this%atmosphere_top    = atmosphere_top
    this%number_of_layers  = number_of_layers

  end subroutine extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Gets the reference element for this extrusion given a particular
  !>        base mesh.
  !>
  !> @param[in] mesh Base mesh object.
  !> @param[out] reference_element Shape of a 3D element given the extrusion.
  !>
  subroutine get_reference_element( this, mesh, reference_element )

    use reference_element_mod, only : reference_prism_type, &
                                      reference_cube_type
    implicit none

    class(extrusion_type),   intent(in) :: this
    class(global_mesh_type), intent(in) :: mesh
    class(reference_element_type), &
                             intent(out), allocatable :: reference_element

    type(reference_prism_type) :: reference_prism
    type(reference_cube_type)  :: reference_cube

    select case (mesh%get_nverts_per_cell())
      case (3)
        allocate( reference_element, source=reference_prism_type() )
      case (4)
        allocate( reference_element, source=reference_cube_type() )
      case default
        write( log_scratch_space, &
              '("Base mesh with ", I0, " vertices per cell not supported.")' &
             ) mesh%get_nverts_per_cell()
        call log_event( log_scratch_space, log_level_error )
    end select

  end subroutine get_reference_element

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Gets the bottom of the atmosphere or the surface of the planet.
  !>
  !> @return Bottom of the atmosphere in meters.
  !>
  function get_atmosphere_bottom( this ) result(bottom)

    implicit none

    class(extrusion_type), intent(in) :: this
    real(r_def) :: bottom

    bottom = this%atmosphere_bottom

  end function get_atmosphere_bottom

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Gets the top of the atmosphere.
  !>
  !> @return Top of the atmosphere in meters.
  !>
  function get_atmosphere_top( this ) result(top)

    implicit none

    class(extrusion_type), intent(in) :: this
    real(r_def) :: top

    top = this%atmosphere_top

  end function get_atmosphere_top

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Gets the number of layers in the atmosphere.
  !>
  !> @return Number of layers.
  !>
  function get_number_of_layers( this ) result(layers)

    implicit none

    class(extrusion_type), intent(in) :: this
    integer(i_def) :: layers

    layers = this%number_of_layers

  end function get_number_of_layers

end module extrusion_mod
