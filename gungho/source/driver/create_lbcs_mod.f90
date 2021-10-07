!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create LBC fields.
!> @details Create LBC field collection and add fields.
module create_lbcs_mod

  use constants_mod,              only : i_def, l_def, str_def
  use log_mod,                    only : log_event,             &
                                         log_scratch_space,     &
                                         LOG_LEVEL_INFO
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use fs_continuity_mod,          only : W0, W2, W3, Wtheta, W2h
  use function_space_mod,         only : function_space_type
  use mr_indices_mod,             only : nummr,                      &
                                         mr_names
  use linked_list_mod,            only : linked_list_type
  use lfric_xios_time_axis_mod,   only : time_axis_type,        &
                                         update_interface
  use lfric_xios_read_mod,        only : read_field_time_var
  use init_time_axis_mod,         only : init_time_axis,        &
                                         setup_field
  use initialization_config_mod,  only : lbc_option,            &
                                         lbc_option_analytic,   &
                                         lbc_option_file

  implicit none

  public  :: create_lbc_fields

  contains

  !> @brief   Create and add LBC fields.
  !> @details Create the lateral boundary condition field collection.
  !!          On every timestep these fields will be updated and used by the
  !!          limited area model.
  !> @param[in]     mesh_id           The current 3d mesh identifier.
  !> @param[in,out] depository        Main collection of all fields in memory.
  !> @param[in,out] prognostic_fields The prognostic variables in the model.
  !> @param[in,out] lbc_fields        The lbc_fields used on every timestep.
  subroutine create_lbc_fields( mesh_id, depository, prognostic_fields, &
                                lbc_fields, lbc_times_list )

    implicit none

    integer(i_def),              intent(in)    :: mesh_id
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(inout) :: lbc_fields
    type(linked_list_type),      intent(out)   :: lbc_times_list

    logical(l_def)                             :: checkpoint_restart_flag
    procedure(update_interface), pointer       :: tmp_update_ptr => null()

    type(time_axis_type), save                 :: lbc_time_axis
    logical(l_def),   parameter                :: cyclic=.false.
    logical(l_def),   parameter                :: interp_flag=.true.
    character(len=*), parameter                :: axis_id="lbc_axis"
    character(str_def)                         :: name
    integer(i_def)                             :: imr

    call log_event( 'Create LBC fields', LOG_LEVEL_INFO )

    lbc_fields = field_collection_type( name='lbc_fields' )

    ! Checkpoint all necessary LBC fields.
    ! At the moment, LBC fields will be set analytically and so all fields
    ! need to be checkpointed.
    ! When we add an option to read in LBC fields from file,  some
    ! LBC fields will no longer need to be checkpointed and this flag can be
    ! set appropriately for each field.

    if ( lbc_option == lbc_option_analytic ) then
      checkpoint_restart_flag = .true.

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_theta", Wtheta, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_u", W2, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_rho", W3, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_exner", W3, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_diff", W2, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_driving", W2, mesh_id, checkpoint_restart_flag )

      do imr = 1, nummr
        name = trim('lbc_' // adjustl(mr_names(imr)) )
        call setup_field( lbc_fields, depository, prognostic_fields, &
           name, wtheta, mesh_id, checkpoint_restart_flag )
      enddo

   else if ( lbc_option == lbc_option_file ) then
      checkpoint_restart_flag = .false.
      ! Set pointer to time axis read behaviour
      tmp_update_ptr => read_field_time_var

      call init_time_axis( "lbc_time", lbc_time_axis, axis_id, &
                           cyclic, interp_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_theta", Wtheta, mesh_id, checkpoint_restart_flag,       &
          time_axis=lbc_time_axis  )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_rho", W3, mesh_id, checkpoint_restart_flag,             &
          time_axis=lbc_time_axis  )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_exner", W3, mesh_id, checkpoint_restart_flag,           &
          time_axis=lbc_time_axis  )

      call lbc_time_axis%set_update_behaviour(tmp_update_ptr)
      call lbc_times_list%insert_item(lbc_time_axis)

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_h_u", W2h, mesh_id, checkpoint_restart_flag,            &
          time_axis=lbc_time_axis  )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_v_u", Wtheta, mesh_id, checkpoint_restart_flag,         &
          time_axis=lbc_time_axis  )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "lbc_u", W2, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_diff", W2, mesh_id, checkpoint_restart_flag )

      call setup_field( lbc_fields, depository, prognostic_fields, &
         "boundary_u_driving", W2, mesh_id, checkpoint_restart_flag )

    end if

  end subroutine create_lbc_fields

end module create_lbcs_mod
