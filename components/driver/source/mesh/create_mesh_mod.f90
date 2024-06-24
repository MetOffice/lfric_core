!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Functions/Routines related to creating a <mesh_object_type>
module create_mesh_mod

  use constants_mod, only: i_def, str_def, r_def, imdi, &
                           str_max_filename
  use log_mod,       only: log_event,         &
                           log_scratch_space, &
                           LOG_LEVEL_DEBUG,   &
                           LOG_LEVEL_ERROR,   &
                           LOG_LEVEL_INFO

  use extrusion_mod,       only: extrusion_type,           &
                                 uniform_extrusion_type,   &
                                 geometric_extrusion_type, &
                                 quadratic_extrusion_type
  use local_mesh_mod,      only: local_mesh_type
  use mesh_mod,            only: mesh_type
  use ugrid_mesh_data_mod, only: ugrid_mesh_data_type

  use local_mesh_collection_mod,  only: local_mesh_collection
  use mesh_collection_mod,        only: mesh_collection

  ! Configuration modules
  use extrusion_config_mod, only: METHOD_UNIFORM,   &
                                  METHOD_GEOMETRIC, &
                                  METHOD_QUADRATIC

  implicit none

  private
  public :: create_extrusion, create_mesh

  interface create_mesh
    module procedure create_mesh_single
    module procedure create_mesh_multiple
  end interface create_mesh

contains


!> @brief  Creates vertical mesh extrusion.
!> @return  Resulting extrusion object
function create_extrusion( extrusion_method, &
                           domain_top,       &
                           domain_bottom,    &
                           n_layers,         &
                           extrusion_id ) result(new)

  implicit none

  class(extrusion_type), allocatable :: new

  integer(i_def), intent(in) :: extrusion_method
  real(r_def),    intent(in) :: domain_top
  real(r_def),    intent(in) :: domain_bottom
  integer(i_def), intent(in) :: n_layers
  integer(i_def), intent(in) :: extrusion_id

  if (allocated(new)) deallocate(new)

  select case (extrusion_method)
    case (METHOD_UNIFORM)
      allocate( new, source=uniform_extrusion_type(        &
                                domain_bottom, domain_top, &
                                n_layers, extrusion_id ) )
    case (METHOD_QUADRATIC)
      allocate( new, source=quadratic_extrusion_type(      &
                                domain_bottom, domain_top, &
                                n_layers, extrusion_id ) )
    case (METHOD_GEOMETRIC)
      allocate( new, source=geometric_extrusion_type(      &
                                domain_bottom, domain_top, &
                                n_layers, extrusion_id ) )
    case default
      call log_event("Invalid method for simple extrusion", LOG_LEVEL_ERROR)
  end select

end function create_extrusion


!> @brief    Generates mesh object types.
!> @details  Creates multiple mesh objects from local_mesh_type objects
!!           held in the application local_mesh_collection and the
!!           specified extrusion.
!!
!> @param[in]  local_mesh_names  Names of the local_mesh_types to extrude.
!> @param[in]  extrusion         Extrusion to employ.
!> @param[in]  alt_name          Optional, Alternative names for the
!!                               extruded meshes, defaults to local_mesh_names
!!                               if absent.
subroutine create_mesh_multiple( local_mesh_names, &
                                 extrusion,        &
                                 alt_name )
  implicit none

  character(str_def),    intent(in) :: local_mesh_names(:)
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: alt_name(:)

  ! Local variables
  integer(i_def) :: i
  character(str_def), allocatable :: names(:)

  if (present(alt_name)) then

    if ( size(alt_name) /= size(local_mesh_names) ) then
      write(log_scratch_space, '(A)')                          &
          'Number of alternative mesh names does not match '// &
          'number of requested meshes.'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    allocate(names, source=alt_name)

  else

    allocate(names, source=local_mesh_names)

  end if

  do i=1, size(local_mesh_names)
    call create_mesh_single( local_mesh_names(i), &
                             extrusion,           &
                             alt_name=names(i) )
  end do

  deallocate(names)

end subroutine create_mesh_multiple

!> @brief    Generates a single mesh_type object.
!> @details  Instantiates a mesh_type object and adds it to the applications
!!           mesh collection. Multiple meshes may be generated in the model
!!           based on the same global mesh but with differing extrusions.
!!
!> @param[in]  local_mesh_name  Name of local_mesh_type object in
!!                              application local_mesh_collection.
!> @param[in]  extrusion        Extrusion to employ for this mesh_type object
!> @param[in]  alt_name         Optional, Alternative name for the
!!                              extruded mesh, defaults to local_mesh_name
!!                              if absent.
subroutine create_mesh_single( local_mesh_name, &
                               extrusion,       &
                               alt_name )

  implicit none

  character(str_def),    intent(in) :: local_mesh_name
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: alt_name

  type(local_mesh_type), pointer :: local_mesh_ptr => null()

  type(mesh_type)        :: mesh
  integer(kind=i_def)    :: mesh_id
  character(len=str_def) :: name

  if ( .not. present(alt_name) ) then
    name = local_mesh_name
  else
    name = alt_name
  end if


  ! 1.0 Check if mesh_type already exists.
  !===============================================
  if ( mesh_collection%check_for(name) ) then
    write(log_scratch_space,'(A)')                          &
        'No action taken: Mesh '//trim(name)//' already '// &
        'exists in the program mesh_collection object.'
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    return
  end if


  ! 2.0 Extrude the local_mesh_object.
  !===============================================
  local_mesh_ptr => local_mesh_collection%get_local_mesh(local_mesh_name)

  if ( .not. associated(local_mesh_ptr) ) then
    write(log_scratch_space,'(A)')                                &
        'Specified local mesh object ('//trim(local_mesh_name)//  &
        ') was not found in the program local_mesh_collection '// &
        'object.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  mesh = mesh_type( local_mesh_ptr, extrusion, mesh_name=name )

  mesh_id = mesh_collection%add_new_mesh( mesh )
  call mesh%clear()


  ! 3.0 Report on mesh_type creation.
  !===============================================
  write(log_scratch_space,'(A,I0,A)')                 &
      '   ... "'//trim(name)//'"(id:', mesh_id,') '// &
      'based on mesh "'//trim(local_mesh_name)//'"'

  if (mesh_id /= imdi) then
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else
    write(log_scratch_space,'(A,I0,A)') &
        trim(log_scratch_space)//' (FAILED)'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine create_mesh_single

end module create_mesh_mod
