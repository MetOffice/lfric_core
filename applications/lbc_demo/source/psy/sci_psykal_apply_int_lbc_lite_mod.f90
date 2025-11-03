!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Provides an implementation of the PSy layer
!>
!> @details Contains hand-rolled versions of the PSy layer that can be used for
!>          simple testing and development of the scientific code

module psykal_apply_int_lbc_lite_mod

  use constants_mod, only: r_def, i_def
  use log_mod,       only: log_event, log_level_error, log_scratch_space
  use mesh_mod,      only: mesh_type
  use mesh_map_mod,  only: mesh_map_type
  use mesh_mod,      only: mesh_type

  use integer_field_mod,  only: integer_field_type, &
                                integer_field_proxy_type
  use function_space_mod, only: function_space_type

  implicit none

  private

  public :: invoke_apply_int_lbc_kernel_type

contains

  !> @brief Kernel to update field with LBC data.
  !> @description Loops over LBC field cells and identifies
  !>              the linked cell on the LAM field to be passed
  !>              to the kernel.
  !>
  !> While this lite code is written for the case of updating
  !> LAM fields with LBC fields on an LBC mesh with appropriate
  !> 1-to-1 mapping. The "LBC mesh" could be any mesh as long as
  !> it has a 1-to-1 mapping to the mesh used by the "LAM field".
  !>
  !> Assumptions:
  !> * LBC-LAM pairing will only have a single 1-to-1 mapping from
  !>   LBC to LAM mesh cells.
  !> * Once read/initialised, LBC-fields will be read only. Currently on
  !>   WCHI label, will LBC fields need their own LABEL?
  !> * Not sure if LBCs require a stencil or not? Don't think
  !>   so, though would be best check with Christine Johnson.
  !>
  !! @param[in,out] lam_field  Integer LAM field to be updated
  !! @param[in]     lbc_field  Integer LBC field to apply to integer LAM field
  subroutine invoke_apply_int_lbc_kernel_type( lam_field, lbc_field )

    use apply_int_lbc_kernel_mod, only: apply_int_lbc_kernel_code

    implicit none

    type(integer_field_type), intent(inout) :: lam_field ! Field to have LBCs updated
    type(integer_field_type), intent(in)    :: lbc_field ! Driver model LBC field


    integer(i_def) :: lam_id   ! LAM field local LAM mesh cell id
    integer(i_def) :: lbc_id   ! LBC field local LBC mesh cell id
    integer(i_def) :: nlayers
    integer(i_def) :: ndata
    logical        :: ndata_first

    type(integer_field_proxy_type) :: lbc_field_proxy
    type(integer_field_proxy_type) :: lam_field_proxy

    integer(i_def), pointer :: map_adspc1_lam_field(:,:) => null() ! LAM Cell DoF-map
    integer(i_def), pointer :: map_adspc2_lbc_field(:,:) => null() ! LAM Cell DoF-map
    integer(i_def), pointer :: cell_map_lbc_field(:,:,:) => null() ! Inter-Grid cell map
                                                                   ! LBC as source

    type(function_space_type), pointer :: fs

    integer(i_def) :: ndf_adspc1_lam_field, undf_adspc1_lam_field
    integer(i_def) :: ndf_adspc2_lbc_field, undf_adspc2_lbc_field

    type(mesh_type),     pointer :: mesh_lam_field           ! LAM field mesh
    type(mesh_type),     pointer :: mesh_lbc_field           ! LBC field mesh
    type(mesh_map_type), pointer :: mmap_lbc_field_lam_field ! Mesh map for
                                                             ! LBC mesh -> LAM mesh

    integer(i_def)              :: icolour, ncolour, cell
    integer(i_def), pointer     :: cmap(:,:)
    integer(i_def), allocatable :: last_edge_cell_all_colours(:)

    nullify(mesh_lam_field)
    nullify(mesh_lbc_field)
    nullify(mmap_lbc_field_lam_field)
    nullify(map_adspc1_lam_field)
    nullify(map_adspc2_lbc_field)
    nullify(cell_map_lbc_field)
    nullify(fs)

    !
    ! Initialise field and/or operator proxies
    !
    lam_field_proxy = lam_field%get_proxy()
    lbc_field_proxy = lbc_field%get_proxy()

    !
    ! Initialise number of layers
    !
    nlayers = lam_field_proxy%vspace%get_nlayers()

    !
    ! Account for multidata fields
    ! LBC and LAM field should have same
    ! data layout and ndata points
    !
    fs => lam_field%get_function_space()
    ndata = fs%get_ndata()
    ndata_first = fs%is_ndata_first()
    fs => lbc_field%get_function_space()

    if ( ndata /= fs%get_ndata() ) then
      write(log_scratch_space, '(A)') &
          'LBC/LAM fields have differing ndata points.'
      call log_event(log_scratch_space, log_level_error)
    end if

    if (ndata_first .neqv. fs%is_ndata_first() ) then
      write(log_scratch_space, '(A)') &
          'LBC/LAM multidata fields have different multi-data layouts.'
      call log_event(log_scratch_space, log_level_error)
    end if

    !
    ! Clean halos for any LAM field. If the LAM is being nudged
    ! the value will need to be correct in order to set the increments.
    !
    call lam_field_proxy%set_dirty()

    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_lam_field => lam_field_proxy%vspace%get_mesh()
    mesh_lbc_field => lbc_field_proxy%vspace%get_mesh()

    !
    ! Get the colourmap
    !
    ncolour =  mesh_lbc_field%get_ncolours()
    cmap    => mesh_lbc_field%get_colour_map()

    !
    ! Access the inter-grid cell map
    !
    mmap_lbc_field_lam_field => mesh_lbc_field%get_mesh_map(mesh_lam_field)
    cell_map_lbc_field       => mmap_lbc_field_lam_field%get_whole_cell_map()

    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_lam_field => lam_field_proxy%vspace%get_whole_dofmap()
    map_adspc2_lbc_field => lbc_field_proxy%vspace%get_whole_dofmap()

    !
    ! Initialise number of DoFs for adspc1_target_field
    !
    ndf_adspc1_lam_field  = lam_field_proxy%vspace%get_ndf()
    undf_adspc1_lam_field = lam_field_proxy%vspace%get_undf()

    !
    ! Initialise number of DoFs for adspc2_source_field
    !
    ndf_adspc2_lbc_field  = lbc_field_proxy%vspace%get_ndf()
    undf_adspc2_lbc_field = lbc_field_proxy%vspace%get_undf()

    !
    ! Get colours on lbc mesh to avoid write contention on lam mesh.
    ! Use edge cells as LBCs have no halos.
    !
    last_edge_cell_all_colours = mesh_lbc_field%get_last_edge_cell_all_colours()

    !
    ! Loop over LBC mesh local cell ids on partition
    ! Colouring from LBC mesh should avoid contention on LAM mesh
    !
    do icolour=1, ncolour
      !$omp parallel default(shared), private(cell, lbc_id, lam_id)
      !$omp do schedule(static)
      do cell=1, last_edge_cell_all_colours(icolour)

        lbc_id = cmap(icolour, cell)
        lam_id = cell_map_lbc_field( 1, 1, lbc_id )

        call apply_int_lbc_kernel_code( nlayers,                        &
                                        lam_field_proxy%data,           &
                                        lbc_field_proxy%data,           &
                                        undf_adspc1_lam_field,          &
                                        map_adspc1_lam_field(:,lam_id), &
                                        ndf_adspc1_lam_field,           &
                                        ndata,                          &
                                        ndata_first,                    &
                                        undf_adspc2_lbc_field,          &
                                        map_adspc2_lbc_field(:,lbc_id) )

      end do
      !$omp end do
      !$omp end parallel
    end do

  end subroutine invoke_apply_int_lbc_kernel_type

end module psykal_apply_int_lbc_lite_mod
