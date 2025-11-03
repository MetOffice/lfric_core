!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> lbc_demo knows what configuration it needs.
!>
module lbc_demo_mod

  implicit none

  private

  character(*), public, parameter ::               &
      required_namelists(5) =  [ 'lbc_demo      ', &
                                 'base_mesh     ', &
                                 'extrusion     ', &
                                 'finite_element', &
                                 'planet        ']

end module lbc_demo_mod
