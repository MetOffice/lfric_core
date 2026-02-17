!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   Temporary module to redirect code that calls generated configuration
!>          module.
!> @details This isolates the change to new namelist access pattern to core repository.
!>          It allows inplementation in lfric apps repo to be done piecemeal
!>
!----------------------------------------------------------------------------
module configuration_mod

  use config_loader_mod, only: read_configuration,   &
                               ensure_configuration, &
                               final_configuration

  implicit none

  private
  public :: read_configuration, ensure_configuration, final_configuration

end module configuration_mod
