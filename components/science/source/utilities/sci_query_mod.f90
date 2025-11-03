!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief Module container for query functions related to science.
module sci_query_mod

  use constants_mod,   only: i_def, str_def
  use global_mesh_mod, only: global_mesh_type
  use local_mesh_mod,  only: local_mesh_type
  use log_mod,         only: log_event, log_scratch_space, &
                             log_level_error

  implicit none

  private
  public  :: valid_for_global_model, &
             check_lbc,              &
             is_lbc

  interface is_lbc
    module procedure is_global_lbc
    module procedure is_local_lbc
  end interface is_lbc

contains

!> @brief  Queries whether a global mesh is valid for a global model domain.
!>         A negative result does not imply that the mesh is valid for a
!>         regional model.
!> @param[in]  global_mesh
!> @return     answwer
function valid_for_global_model( global_mesh ) result ( answer )

   implicit none

   type(global_mesh_type), intent(in) :: global_mesh
   logical :: answer

   if ( global_mesh%is_topology_periodic() .and. &
        global_mesh%is_coord_sys_ll()      .and. &
        global_mesh%is_geometry_spherical() ) then

      ! Note: if these conditions where satisfied from the
      !       planar mesh generator then the model would be
      !       a torus, though as this is not supported, this
      !       is only .true. for a sphere.
      answer = .true.
  else
      answer = .false.
  end if

end function valid_for_global_model

!---------------------------------------------------------------------------
!> @brief  Queries if the global mesh object is an LBC mesh
!>
!> @return answer .true. if mesh is for storing Lateral Boundary Conditions.
!>
function is_global_lbc( mesh ) result ( answer )

  implicit none

  type(global_mesh_type), intent(in) :: mesh

  character(str_def) :: mesh_name
  logical            :: answer

  mesh_name = mesh%get_mesh_name()
  answer    = check_lbc( mesh_name )

end function is_global_lbc

!---------------------------------------------------------------------------
!> @brief  Queries if the local mesh object is an LBC mesh
!>
!> @return answer .true. if mesh is for storing Lateral Boundary Conditions.
!>
function is_local_lbc( mesh ) result ( answer )

  implicit none

  type(local_mesh_type), intent(in) :: mesh

  character(str_def) :: mesh_name
  logical            :: answer

  mesh_name = mesh%get_mesh_name()
  answer    = check_lbc( mesh_name )

end function is_local_lbc


!---------------------------------------------------------------------------
!> @brief  Queries if the mesh_name provided contains the lbc-suffix
!>
!> @return answer .true. if mesh is for storing Lateral Boundary Conditions.
!>
function check_lbc( mesh_name ) result ( answer )

  implicit none

  character(str_def), intent(in) :: mesh_name

  logical :: answer

  character(4), parameter  :: lbc_tag = '-lbc'

  integer(i_def) :: mesh_name_length
! character(4)   :: tmp_str

  ! Check to see if the loaded mesh is an LBC mesh
  ! where the lbc_suffix would be appended.
  mesh_name_length = len(trim(mesh_name))
  answer = .false.

  if (mesh_name_length > 4) then
    if ( index( mesh_name, lbc_tag ) /= 0 ) then
      answer = .true.
    end if
  end if

end function check_lbc

end module sci_query_mod
