!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Main program for running the CMA operator tests in cma_test_mod.x90.
!>@details To run a particular test, provide the appropriate command line
!>         argument of the form --test_XXX where a list of possible values
!>         for XXX is given below. Running with --test_all will carry out all
!>         tests.
!>
!>         Any configuration data is read from the namelist file
!>         cma_test_configuration.nml.
!>
!>         Note that the test of the diag_DhMDhT term currently only works
!>         reliably sequentially.

program cma_test

  use, intrinsic :: iso_fortran_env,  only : real64
  use base_mesh_config_mod,           only : geometry,           &
                                             geometry_spherical, &
                                             prime_mesh_name
  use cma_test_algorithm_mod,         only : cma_test_init,                  &
                                             test_cma_apply_mass_p,          &
                                             test_cma_apply_mass_v,          &
                                             test_cma_apply_div_v,           &
                                             test_cma_multiply_div_v_mass_v, &
                                             test_cma_multiply_grad_v_div_v, &
                                             test_cma_add,                   &
                                             test_cma_apply_inv,             &
                                             test_cma_diag_DhMDhT
  use constants_mod,                  only : i_def, r_def, pi
  use derived_config_mod,             only : set_derived_config
  use yaxt,                           only : xt_initialize, xt_finalize
  use mpi_mod,                        only : initialise_comm, store_comm, &
                                             finalise_comm, &
                                             get_comm_size, get_comm_rank, &
                                             global_sum
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use fs_continuity_mod,              only : W0,W1,W2,W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use configuration_mod,              only : read_configuration, &
                                             ensure_configuration
  use create_mesh_mod,                only : init_mesh
  use local_mesh_collection_mod,      only : local_mesh_collection, &
                                             local_mesh_collection_type
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             initialise_logging, &
                                             finalise_logging, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO
  use mesh_collection_mod,            only : mesh_collection, &
                                             mesh_collection_type
  use mesh_mod,                       only : mesh_type
  use mesh_collection_mod,            only : mesh_collection
  use planet_config_mod,              only : radius

  implicit none

  ! MPI communicator
  integer(kind=i_def) :: comm

  ! Meshes index
  integer(kind=i_def) :: mesh_id

  ! Number of processes and local rank
  integer(kind=i_def) :: total_ranks, local_rank
  ! Filename to read namelist from
  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag
  ! Usage message to print
  character(len=256) :: usage_message
  ! The following variables are required for checking that the mesh
  ! has a suitable size (i.e. \f$dx \approx dz\f$)
  ! Pointer to mesh
  type(mesh_type), pointer :: mesh
  ! Number of cells of 2d mesh (local and global),
  integer(kind=i_def) :: ncells_2d_local, ncells_2d
  ! Number of vertical layers
  integer(kind=i_def) :: nlayers
  ! vertical domain size
  real   (kind=r_def) :: domain_top
  ! Grid spacing in horizontal and vertical
  real   (kind=r_def) :: dx, dz
  ! Variables for reading configuration from namelist file
  character(*), parameter ::       &
       required_configuration(7) = &
      (/'finite_element      ',    &
        'base_mesh           ',    &
        'formulation         ',    &
        'planet              ',    &
        'extrusion           ',    &
        'partitioning        ',    &
        'domain_size         '/)
  logical              :: okay
  logical, allocatable :: success_map(:)
  integer :: i

  ! Flags which determine the tests that will be carried out
  logical :: do_test_apply_mass_p = .false.
  logical :: do_test_apply_mass_v = .false.
  logical :: do_test_apply_div_v = .false.
  logical :: do_test_multiply_div_v_mass_v = .false.
  logical :: do_test_multiply_grad_v_div_v = .false.
  logical :: do_test_add = .false.
  logical :: do_test_apply_inv = .false.
  logical :: do_test_diag_dhmdht = .false.

  ! Error tolerance for tests
  ! Note: tolerance is for r_def = real64
  !       Tolerance at r_def = real32: generated using spacing() in
  !           ./gungho/source/algorithm/cma_test_algorithm_mod.x90
  real(kind=r_def) :: tolerance
  tolerance = 1.0E-12_r_def

  ! Initialise MPI communicatios and get a valid communicator
  call initialise_comm(comm)

  ! Save the communicator for later use
  call store_comm(comm)

  ! Initialise YAXT
  call xt_initialize(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, 'cma_test')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call log_event( 'CMA functional testing running ...', LOG_LEVEL_INFO )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parse command line parameters
  !
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()
  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
     write(usage_message,*) "Usage: ",trim(program_name), &
          " <namelist filename> " // &
          " test_XXX with XXX in { " // &
          "apply_mass_p, "             // &
          "apply_mass_v, "             // &
          "apply_div_v, "              // &
          "multiply_div_v_mass_v, "    // &
          "multiply_grad_v_div_v, "    // &
          "add, "                      // &
          "apply_inv, "                // &
          "diag_dhmdht, "              // &
          "all"                        // &
          " } "
     call log_event( trim(usage_message), LOG_LEVEL_ERROR )
  end if

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call get_command_argument( 2, dummy, length, status )
  allocate(character(length)::test_flag)
  call get_command_argument( 2, test_flag, length, status )

  ! Choose test case depending on flag provided in the first command
  ! line argument
  select case (trim(test_flag))
  case ("test_apply_mass_p")
     do_test_apply_mass_p = .true.
  case ("test_apply_mass_v")
     do_test_apply_mass_v = .true.
  case ("test_apply_div_v")
     do_test_apply_div_v = .true.
  case ("test_multiply_div_v_mass_v")
     do_test_multiply_div_v_mass_v = .true.
  case ("test_multiply_grad_v_div_v")
     do_test_multiply_grad_v_div_v = .true.
  case ("test_add")
     do_test_add = .true.
  case ("test_apply_inv")
     do_test_apply_inv = .true.
  case ("test_diag_dhmdht")
     do_test_diag_dhmdht = .true.
  case ("test_all")
     do_test_apply_mass_p = .true.
     do_test_apply_mass_v = .true.
     do_test_apply_div_v = .true.
     do_test_multiply_div_v_mass_v = .true.
     do_test_multiply_grad_v_div_v = .true.
     do_test_add = .true.
     do_test_apply_inv = .true.
     do_test_diag_dhmdht = .true.
  case default
     call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  deallocate(program_name)
  deallocate(test_flag)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Load configuration from file
  !
  write( log_scratch_space, '("Reading namelist file: ", A)' ) trim(filename)
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  allocate( success_map(size(required_configuration)) )

  call read_configuration( filename )

  okay = ensure_configuration( required_configuration, success_map )
  if (.not. okay) then
     write( log_scratch_space, '(A)' ) &
            'The following required namelists were not loaded:'
     do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
             log_scratch_space = trim(log_scratch_space) // ' ' &
                                 // required_configuration(i)
     end do
     call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if
  deallocate( success_map )

  call set_derived_config(.false.)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise

  call log_event( 'Initialising harness', LOG_LEVEL_INFO )

  allocate( local_mesh_collection, &
            source = local_mesh_collection_type() )
  allocate( mesh_collection, &
            source=mesh_collection_type() )

  call init_mesh( local_rank, total_ranks, mesh_id )

  ! Work out grid spacing, which should be of order 1
  mesh => mesh_collection%get_mesh( mesh_id )
  ncells_2d_local = mesh%get_ncells_2d()

  ! Ensure that a spherical geometry is used (otherwise tests are too simple)
  if (geometry /= geometry_spherical) then
     call log_event( "Geometry has to be spherical", &
                     LOG_LEVEL_ERROR )
  end if

  ! Work out total number of cells
  call global_sum(ncells_2d_local, ncells_2d)

  ! Check that the grid spacings are of order 1 (between 0.1 and 10).
  ! Otherwise the derivative and mass terms in the
  ! Helmholtz-operator are not well balanced and errors could go unnoticed
  nlayers = mesh%get_nlayers()
  domain_top = mesh%get_domain_top()
  dx = sqrt(4.*pi/ncells_2d)*radius
  dz = domain_top / nlayers

  if ( ( dx < 0.1_r_def) .or. ( dx > 10.0_r_def) ) then
     call log_event( "Average dx has to be in range 0.1 ... 10.0", &
                     LOG_LEVEL_ERROR )
  end if
  if ( ( dz < 0.1_r_def) .or. ( dz > 10.0_r_def) ) then
     call log_event( "Average dz has to be in range 0.1 ... 10.0", &
                     LOG_LEVEL_ERROR )
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Do tests
  !
  call log_event( 'Initialising test', LOG_LEVEL_INFO )

  ! Initialise CMA test module
  call cma_test_init(mesh_id)

  ! Run all requested tests
  if (do_test_apply_mass_p)             call test_cma_apply_mass_p(tolerance)
  if (do_test_apply_mass_v)             call test_cma_apply_mass_v(tolerance)
  if (do_test_apply_div_v)              call test_cma_apply_div_v(tolerance)
  if (do_test_multiply_div_v_mass_v)    call test_cma_multiply_div_v_mass_v(tolerance)
  if (do_test_multiply_grad_v_div_v)    call test_cma_multiply_grad_v_div_v(tolerance)
  if (do_test_add)                      call test_cma_add(tolerance)
  if (do_test_apply_inv)                call test_cma_apply_inv(tolerance)
  if (do_test_diag_dhmdht)              call test_cma_diag_DhMDhT(tolerance)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down
  !
  call log_event( ' CMA functional testing completed ...', LOG_LEVEL_INFO )

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise MPI communications
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end program cma_test
