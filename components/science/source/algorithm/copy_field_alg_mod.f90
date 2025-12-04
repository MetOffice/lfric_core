!-----------------------------------------------------------------------------
! Copyright (c) 2024,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which
! you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

!> @brief Provides an implementation of mixed precisions field copy for
!  real32 and real64 fields.


  module copy_field_alg_mod

    use, intrinsic :: iso_fortran_env, only: real32, real64

    use sci_psykal_builtin_light_mod, only : &
                            invoke_copy_field_32_64, &
                            invoke_copy_field_64_32, &
                            invoke_copy_field_32_32, &
                            invoke_copy_field_64_64

    use field_real32_mod,      only : field_real32_type
    use field_real64_mod,      only : field_real64_type

    implicit none

    interface copy_field
       module procedure &
          copy_field_32_32, copy_field_64_64, &
          copy_field_32_64, copy_field_64_32
    end interface copy_field

  contains
    !--------------------------------------------------------------
    !> Field copy routines for all exisiting combinations of field
    !> real types. Algorithms use native Fortran reals, essentially
    !> seeing through the rdef, rtran, etc. types.
    !> The copy routines are prepared for parallelism.
    subroutine copy_field_32_64(fsrce_32, fdest_64)
       implicit none
       type(field_real32_type), intent(in)     :: fsrce_32
       type(field_real64_type), intent(inout)  :: fdest_64
       call invoke_copy_field_32_64(fsrce_32, fdest_64)
    end subroutine copy_field_32_64

    subroutine copy_field_64_32(fsrce_64, fdest_32)
       implicit none
       type(field_real64_type), intent(in)     :: fsrce_64
       type(field_real32_type), intent(inout)  :: fdest_32
       call invoke_copy_field_64_32(fsrce_64, fdest_32)
    end subroutine copy_field_64_32

    subroutine copy_field_32_32(fsrce_32, fdest_32)
       implicit none
       type(field_real32_type), intent(in)     :: fsrce_32
       type(field_real32_type), intent(inout)  :: fdest_32
       call invoke_copy_field_32_32(fsrce_32, fdest_32)
    end subroutine copy_field_32_32

    subroutine copy_field_64_64(fsrce_64, fdest_64)
       implicit none
       type(field_real64_type), intent(in)     :: fsrce_64
       type(field_real64_type), intent(inout)  :: fdest_64
       call invoke_copy_field_64_64(fsrce_64, fdest_64)
    end subroutine copy_field_64_64

  end module copy_field_alg_mod
