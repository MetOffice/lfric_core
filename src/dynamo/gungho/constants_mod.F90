!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Define various constants for the application.
!>
!> @details Various computational, physical and geometrical constants are
!>          defined in this module. Their values are also set here.
module constants_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, int32

  implicit none

  ! Define default application-defined kinds for all intrinsic data types

  !> @name Set up default kinds for real and double-precision variables.
  !> @{
  real,             private :: r_val  !< A native real used to compute kind of native real.
  double precision, private :: dp_val !< A native double-precision used to compute kind of native dp.

  integer, parameter :: r_def     = real64 !< Default real kind for application.
  integer, parameter :: r_single  = real32 !< Default single precision real kind for application.
  integer, parameter :: r_double  = real64 !< Default double precision real kind for application.

  integer, parameter :: r_native  = kind(r_val)  !< Native kind for real.
  integer, parameter :: dp_native = kind(dp_val) !< Native kind for double precision.
  !> @}

  !> @name Complex
  !> @{
  !> @}

  !> @name Set up default kinds for integers.
  !> @{
  integer, private   :: i_val     !< A native integer used to compute kind of native integer.

  integer, parameter :: i_def     = int32       !< Default integer kind for application.
  integer, parameter :: i_native  = kind(i_val) !< Native kind for integer.
  !> @}

  !> @name Set up default kinds for logicals.
  !> @{
  logical, private   :: l_val     !< A native logical used to compute kind of native logical.

  integer, parameter :: l_def     = kind(l_val) !< Default logical kind for application.
  integer, parameter :: l_native  = kind(l_val) !< Native kind for logical.
  !> @}

  !> @name Set up default kinds for character variables.
  !> @{
  character, private :: c_val     !< A native character used to compute kind of native character.

  integer, parameter :: c_def     = kind(c_val) !< Default character kind for application.
  integer, parameter :: c_native  = kind(c_val) !< Native kind for character.
  !> @}

  !> @name Set up default lengths for string variables.
  !> @{
  integer, parameter :: str_def          = 128 !< Default string length for normal strings.
  integer, parameter :: str_long         = 255 !< Default length of long string.
  integer, parameter :: str_max_filename = 255 !< Default maximum length of a file-name.
  !> @}

  !> @name Platform constants
  !> @{
  real(kind=r_def), parameter :: LARGE_REAL = huge(0.0_r_def) !< The largest number of kind r_def that is not an infinity.
  !> @}

  !> @name Numerical constants
  !> @{
  real(kind=r_def), parameter :: EPS = 3.0e-15_r_def !< Relative precision: if (abs(x-y) > EPS) then assume x==y.
  !> @}

  !> @name Mathematical constants
  !> @{
  real(kind=r_def), parameter :: PI  = 3.141592654_r_def !< Value of pi.
  !> @}

  !> @name Physical constants
  !> @{
  real(kind=r_def), parameter :: GRAVITY = 9.80665_r_def !< Surface equatorial value of gravity [m/s^2].
  real(kind=r_def), parameter :: Rd = 287.05_r_def       !< Gas constant for dry air [J/(kg K)].
  real(kind=r_def), parameter :: Cp = 1005.0_r_def       !< Specific heat of dry air at constant pressure [J/(kg K)].
  real(kind=r_def), parameter :: Cv = Cp - Rd            !< Specific heat of dry air at constant volume [J/(kg K)].
  real(kind=r_def), parameter :: KAPPA = Rd/Cp           !< Ratio of Rd and Cp [dimensionless].
  real(kind=r_def), parameter :: P_ZERO = 100000.0_r_def !< Reference surface pressure [Pa].
  real(kind=r_def), parameter :: N_SQ = 0.0001_r_def     !< The square of Brunt-Vaisala frequency [1/s^2].
  real(kind=r_def), parameter :: OMEGA_UNSCALED = 7.292116E-5_r_def !< Angular velocity [rad/s].
  real(kind=r_def), parameter :: EARTH_RADIUS_UNSCALED = 6371229.0_r_def !< Radius of the Earth [m].
  !> @}

  !> @name Small Earth scalings
  !> @{
  real(kind=r_def), parameter :: EARTH_SCALING = 125.0_r_def !< Scaling factor to modify Earth parameters.
  real(kind=r_def)            :: omega = OMEGA_UNSCALED*EARTH_SCALING !< Scale up rotation [rad/s].
  real(kind=r_def)            :: earth_radius = EARTH_RADIUS_UNSCALED/EARTH_SCALING !< Scale down Earth radius [m].
  !> @}

  !> @name Enumeration of the available choices for the linear solver 
  !> @{
  integer (kind=i_def), parameter :: CG_SOLVER     = 1 !< Conjugate Gradient solver option.
  integer (kind=i_def), parameter :: BICG_SOLVER   = 2 !< BiCGSTAB solver option.
  integer (kind=i_def), parameter :: JACOBI_SOLVER = 3 !< Jacobi solver option.
  integer (kind=i_def), parameter :: GMRES_SOLVER  = 4 !< Generalized Minimal RESidual solver option.
  integer (kind=i_def), parameter :: GCR_SOLVER    = 5 !< Generalized Conjugate Residual solver option.
  !> @}
  !> @name Linear solver constants
  !> @{
  integer (kind=i_def), parameter :: MAX_ITER = 99 !< Maximum iteration number for solver.
  integer (kind=i_def), parameter :: SOLVER_OPTION = BICG_SOLVER !< Choice of solver from the above list
                                                                 !< (currently hard-wired).
  integer (kind=i_def), parameter :: NO_PRE_COND       = -1 !< No preconditioner option.
  integer (kind=i_def), parameter :: DIAGONAL_PRE_COND = 1  !< Diagonal preconditioner option.
  real(kind=r_def),     parameter :: SOLVER_TOL = 1.0e-4_r_def !< Relative tolerance of solver.
  integer (kind=i_def), parameter :: GCRK  = 32                !< Dimension of the approximate Krylov subspace.
                                                               !< In other words, it is the number of potential 
                                                               !< residual vectors to calculate at each iteration 
                                                               !< of the solver

  integer (kind=i_def), parameter :: TRI  = 1                  !< For triangular reference elements
  integer (kind=i_def), parameter :: QUAD = 2                  !< For quadrilateral reference elements

  ! Missing data indicators
  real    (r_def), parameter :: RMDI = -huge(0.0_r_def)        !< Missing data indicator value for real numbers
  integer (i_def), parameter :: IMDI = -huge(0_i_def)          !< Missing data indicator value for integer numbers

  ! Grid Types
  integer(i_def), parameter :: PLANE             = 1
  integer(i_def), parameter :: PLANE_BI_PERIODIC = 2
  
  !> @}
  !> @name Formulation switches
  logical, parameter :: L_NONLINEAR     = .true.  !< Solve the full nonlinear equation set
  logical, parameter :: L_SEMI_IMPLICIT = .true.  !< Use the iterative timestepping method or runge kutta method
  logical, parameter :: L_ROTATING      = .true.  !< Turn on/off Coriolis terms
  logical, parameter :: L_NEWTON_KRYLOV = .false. !< Use Newton-Krylov method to compute lhs
  logical, parameter :: L_SUPG          = .false. !< Use Streamline-Upwind-Petrov-Galerkin method for stabilisation of CG advection
  !> @}

  !> @}
  !> @name Runtime options
  real(kind=r_def),    parameter :: DT = 10.0_r_def !< Timestep in seconds
  integer(kind=i_def), parameter :: NT = 3          !< Number of timesteps to run for
  !> @}
  !> @}
  !> @name Idealised test switches
  logical, parameter :: L_COLD_BUBBLE  = .false. ! Straka density current test (planer domain only)
  logical, parameter :: L_GRAVITY_WAVE = .true.  ! Gravity wave test (either planer or spherical)
  !> @}
  !> @}
  !> @name Iterative timestepping options
  integer, parameter :: N_OUTER_ITER = 2 !< Number of outer (advection) iterations to do
  integer, parameter :: N_INNER_ITER = 2 !< Number of inner (Newton) iterations to do
  real(kind=r_def), parameter :: ALPHA = 0.5_r_def !< Time off-centering parameter
  real(kind=r_def), parameter :: BETA = (1.0_r_def - ALPHA) !< 1 - Time off-centering parameter
  !> @}
 !> @name Enumeration of the available choices for the initialisation profiles
 !> @{
 integer, parameter :: ZERO_WIND                = 0
 integer, parameter :: SOLID_BODY_ROTATION_WIND = 1
 integer, parameter :: CONSTANT_UV_WIND         = 2
 integer, parameter :: CONSTANT_SHEAR_UV_WIND   = 3
 !> @{

 !> @{
 !> @name Initialisation options
 integer(kind=i_def), parameter :: INITIAL_U_PROFILE = zero_wind
 real(kind=r_def),    parameter :: U0 = 0.0_r_def
 real(kind=r_def),    parameter :: V0 = 0.0_r_def
 real(kind=r_def),    parameter :: ROTATION_ANGLE = 0.0_r_def
 !> @{
end module constants_mod

