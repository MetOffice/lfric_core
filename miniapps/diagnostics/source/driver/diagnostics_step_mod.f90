!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the stepping of the diagnostics miniapp
!>
module diagnostics_step_mod

    use clock_mod,                      only : clock_type
    use colours__prognostics__meta_mod, only : colours__prognostics__meta_type
    use colours__diagnostics__meta_mod, only : colours__diagnostics__meta_type
    use colours__non_spatial__meta_mod, only : colours__non_spatial__meta_type
    use constants_mod,                  only : i_def, str_def, r_def
    use field_mod,                      only : field_type
    use field_parent_mod,               only : field_parent_type
    use field_collection_mod,           only : field_collection_type
    use hex_alg_mod,                    only : hex_alg
    use io_config_mod,                  only : write_diag
    use diagnostics_model_data_mod,     only : model_data_type
    use log_mod,                        only : log_event, LOG_LEVEL_INFO
    use mesh_mod,                       only : mesh_type
    use non_spatial_alg_mod,            only : non_spatial_alg
    use spread_alg_mod,                 only : spread_alg

    implicit none

    private
    public diagnostics_step


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Performs time steps.
    !>
    subroutine diagnostics_step( mesh,       &
                                 twod_mesh,  &
                                 model_data, &
                                 clock )

        implicit none

        ! Model run working data set
        type(mesh_type),         pointer, intent(in)    :: mesh      ! included for consistency with gungho
        type(mesh_type),         pointer, intent(in)    :: twod_mesh ! included for consistency with gungho
        type( model_data_type ), target,  intent(inout) :: model_data
        class(clock_type),                intent(in)    :: clock     ! included for consistency with gungho

        type(field_type), pointer :: red => null()
        type(field_type), pointer :: green => null()
        type(field_type), pointer :: blue => null()
        type(field_type), pointer :: hex => null()
        type(field_type), pointer :: mutable_categories => null()
        type(field_type), pointer :: mutable_numbers => null()
        type(field_type), pointer :: immutable_both => null()

        type(colours__prognostics__meta_type) :: prognostics_meta
        type(colours__diagnostics__meta_type) :: diagnostics_meta
        type(colours__non_spatial__meta_type) :: non_spatial_meta
        type(field_collection_type), pointer :: prognostic_fields => null()
        type(field_collection_type), pointer :: non_spatial_fields => null()
        character(str_def) :: hex_id

        ! Demonstrate here that fields can be obtained from their field
        ! collection in model_data, or directly from the depository

        prognostics_meta = colours__prognostics__meta_type()
        diagnostics_meta = colours__diagnostics__meta_type()
        ! Use the metadata to get field collection from model_data
        ! Get individual fields using their unique IDs
        prognostic_fields => model_data%get_field_collection( prognostics_meta%name )

        if (prognostic_fields%field_exists( &
                prognostics_meta%red%get_unique_id() )) then
            red => prognostic_fields%get_field( prognostics_meta%red%get_unique_id() )
            call spread_alg(red)
        end if
        if (prognostic_fields%field_exists( &
                prognostics_meta%green%get_unique_id() )) then
            green => prognostic_fields%get_field( prognostics_meta%green%get_unique_id() )
            call spread_alg(green)
        end if
        if (prognostic_fields%field_exists( &
                prognostics_meta%blue%get_unique_id() )) then
            blue => prognostic_fields%get_field( prognostics_meta%blue%get_unique_id() )
            call spread_alg(blue)
        end if

        ! Get an individual field from depository
        ! hex only exists if red, green and blue are all on so can assume
        ! they are present
        hex_id = diagnostics_meta%hex%get_unique_id()
        if (model_data%depository%field_exists(hex_id)) then
            hex => model_data%depository%get_field( hex_id )
          call hex_alg(red, green, blue, hex)
        end if

        ! Write the fields
        if (write_diag) then
            if(associated(red)) then
              call log_event("Writing " // red%get_name(), LOG_LEVEL_INFO)
              call red%write_field(red%get_name())
            end if
            if(associated(green)) then
              call log_event("Writing " // green%get_name(), LOG_LEVEL_INFO)
              call green%write_field(green%get_name())
            end if
            if(associated(blue)) then
              call log_event("Writing " // blue%get_name(), LOG_LEVEL_INFO)
              call blue%write_field(blue%get_name())
            end if
            if(associated(hex)) then
              call log_event("Writing " // hex%get_name(), LOG_LEVEL_INFO)
              call hex%write_field(hex%get_name())
            end if
        end if

        ! Check non-spatial fields are present
        non_spatial_meta = colours__non_spatial__meta_type()
        if (model_data%field_collection_exists(non_spatial_meta%name)) then
            non_spatial_fields => model_data%get_field_collection( non_spatial_meta%name )

            ! Get fields
            if (non_spatial_fields%field_exists(non_spatial_meta%mutable_categories%get_unique_id()) .and. &
                    non_spatial_fields%field_exists(non_spatial_meta%mutable_numbers%get_unique_id()) .and. &
                    non_spatial_fields%field_exists(non_spatial_meta%immutable_both%get_unique_id())) then
                mutable_categories => non_spatial_fields%get_field( non_spatial_meta%mutable_categories%get_unique_id() )
                mutable_numbers => non_spatial_fields%get_field( non_spatial_meta%mutable_numbers%get_unique_id() )
                immutable_both => non_spatial_fields%get_field( non_spatial_meta%immutable_both%get_unique_id() )

                ! Call the algorithm
                call non_spatial_alg(&
                    mutable_categories, &
                    mutable_numbers, &
                    immutable_both)

                ! Write the fields
                call log_event("Writing " // mutable_categories%get_name(), LOG_LEVEL_INFO)
                call mutable_categories%write_field(mutable_categories%get_name())
                call log_event("Writing " // mutable_numbers%get_name(), LOG_LEVEL_INFO)
                call mutable_numbers%write_field(mutable_numbers%get_name())
                call log_event("Writing " // immutable_both%get_name(), LOG_LEVEL_INFO)
                call immutable_both%write_field(immutable_both%get_name())
            end if
        end if

    end subroutine diagnostics_step

end module diagnostics_step_mod
