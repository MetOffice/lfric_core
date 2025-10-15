!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module sci_psykal_light_mod

  use constants_mod,      only : i_def, r_def, r_double, r_solver
  use field_mod,          only : field_type, field_proxy_type
  use integer_field_mod,  only : integer_field_type, integer_field_proxy_type
  use r_solver_field_mod, only : r_solver_field_type, r_solver_field_proxy_type
  use r_tran_field_mod,   only : r_tran_field_type, r_tran_field_proxy_type

  implicit none

  public

contains

  !-------------------------------------------------------------------------------
  !> This PSyKAl-lite code is required because, currently, PSYclone does not support
  !> the output of scalar variables from kernels.(See PSyclone issue #1818)
  !> This subroutine recovers a scalar value from a field. This is required
  !> as scalars can't currently be written to checkpoint files. The workaround is to
  !> copy the scalar to a field, which may then be checkpointed. On a restart the
  !> scalar value needs to be recovered from the checkpointed field.
    subroutine invoke_getvalue(field, val)
      implicit none
      real(r_def), intent(out)     :: val
      type(field_type), intent(in) :: field
      type(field_proxy_type)       :: field_proxy
      field_proxy = field%get_proxy()
      val = field_proxy%data(1)
    end subroutine invoke_getvalue


  !-----------------------------------------------------------------------------
  !> @brief Performs a halo exchange on a field
  !> @details This PSyKAl-lite code is required to allow an interface to the
  !!          field's underlying halo exchange routine. Accessing a field's
  !!          proxy type directly from an algorithm is not allowed in the API,
  !!          so this subroutine provides an acceptable way to do this.
  !!          NOTE: there is not an associated PSyclone issue for this, as this
  !!          does not relate to a kernel or its metadata. The intention is that
  !!          this routine is used in the algorithm layer as an optimisation to
  !!          prevent a field from performing a small and then large halo
  !!          exchange, instead of just the large exchange. It is not really
  !!          conceivable that PSyclone could ever detect such a situation when
  !!          it is separated across different modules.
  !> @param[in,out] field  The field to perform a halo exchange on
  !> @param[in]     depth  Depth of the halo exchange to be performed
  subroutine invoke_rtran_halo_exchange(field, depth)

    implicit none

    type(r_tran_field_type), intent(inout) :: field
    integer(kind=i_def),     intent(in)    :: depth
    type(r_tran_field_proxy_type)          :: field_proxy

    field_proxy = field%get_proxy()
    if (field_proxy%is_dirty(depth=depth)) then
      call field_proxy%halo_exchange(depth=depth)
    end if

  end subroutine invoke_rtran_halo_exchange

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Perform innerproduct of a r_solver precision field in r_double precision
    subroutine invoke_rdouble_X_innerproduct_X(field_norm, field)

      use scalar_mod,         only: scalar_type
      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type

      implicit none

      real(kind=r_def), intent(out) :: field_norm
      type(r_solver_field_type), intent(in) :: field

      type(scalar_type)                           :: global_sum
      integer(kind=i_def)                         :: df
      real(kind=r_double), allocatable, dimension(:) :: l_field_norm
      integer(kind=i_def)                         :: th_idx
      integer(kind=i_def)                         :: loop0_start, loop0_stop
      integer(kind=i_def)                         :: nthreads
      type(r_solver_field_proxy_type)             :: field_proxy
      integer(kind=i_def)                         :: max_halo_depth_mesh
      type(mesh_type), pointer                    :: mesh => null()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise field and/or operator proxies
      !
      field_proxy = field%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => field_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = field_proxy%vspace%get_last_dof_owned()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      field_norm = 0.0_r_def
      ALLOCATE (l_field_norm(nthreads))
      l_field_norm = 0.0_r_double
      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        l_field_norm(th_idx) = l_field_norm(th_idx) + real(field_proxy%data(df),r_double)**2
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        field_norm = field_norm+real(l_field_norm(th_idx),r_def)
      END DO
      DEALLOCATE (l_field_norm)
      global_sum%value = field_norm
      field_norm = global_sum%get_sum()
      !
    end subroutine invoke_rdouble_X_innerproduct_X


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Perform innerproduct of a r_solver precision field in r_def precision
    subroutine invoke_rdouble_X_innerproduct_Y(field_norm, field1, field2)

      use scalar_mod,         only: scalar_type
      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type

      implicit none

      real(kind=r_def), intent(out) :: field_norm
      type(r_solver_field_type), intent(in) :: field1, field2

      type(scalar_type)                           :: global_sum
      integer(kind=i_def)                         :: df
      real(kind=r_double), allocatable, dimension(:) :: l_field_norm
      integer(kind=i_def)                         :: th_idx
      integer(kind=i_def)                         :: loop0_start, loop0_stop
      integer(kind=i_def)                         :: nthreads
      type(r_solver_field_proxy_type)             :: field1_proxy, field2_proxy
      integer(kind=i_def)                         :: max_halo_depth_mesh
      type(mesh_type), pointer                    :: mesh => null()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise field and/or operator proxies
      !
      field1_proxy = field1%get_proxy()
      field2_proxy = field2%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => field1_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = field1_proxy%vspace%get_last_dof_owned()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      field_norm = 0.0_r_def
      ALLOCATE (l_field_norm(nthreads))
      l_field_norm = 0.0_r_double
      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        l_field_norm(th_idx) = l_field_norm(th_idx) + real(field1_proxy%data(df),r_double)*real(field2_proxy%data(df),r_double)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        field_norm = field_norm+real(l_field_norm(th_idx),r_def)
      END DO
      DEALLOCATE (l_field_norm)
      global_sum%value = field_norm
      field_norm = global_sum%get_sum()
      !
    end subroutine invoke_rdouble_X_innerproduct_Y


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine invoke_inc_rdefX_plus_rsolverY(X, Y)

      use mesh_mod, only: mesh_type

      implicit none

      type(field_type),          intent(inout) :: X
      type(r_solver_field_type), intent(in)    :: Y
      integer(kind=i_def) :: df
      integer(kind=i_def) :: loop0_start, loop0_stop
      type(field_proxy_type) :: X_proxy
      type(r_solver_field_proxy_type) :: Y_proxy
      integer(kind=i_def) :: max_halo_depth_mesh
      type(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      X_proxy = X%get_proxy()
      Y_proxy = Y%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => X_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = X_proxy%vspace%get_last_dof_annexed()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        X_proxy%data(df) = X_proxy%data(df) + real(Y_proxy%data(df),r_def)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop(s)
      !
      CALL X_proxy%set_dirty()
      !
      ! End of set dirty/clean section for above loop(s)
      !
      !
    end subroutine invoke_inc_rdefX_plus_rsolverY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_coarse of the dummy_coarse field of the intermesh kernel
  ! alongside that of that of the fine field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_weights_scalar_inter_element_order_kernel_type(dummy_fine, dummy_coarse, weights_high_low, &
&weights_low_high, qr)
      USE sci_weights_scalar_inter_element_order_kernel_mod, ONLY: weights_scalar_inter_element_order_code
      USE quadrature_xyoz_mod, ONLY: quadrature_xyoz_type, quadrature_xyoz_proxy_type
      USE function_space_mod, ONLY: BASIS, DIFF_BASIS
      USE mesh_map_mod, ONLY: mesh_map_type
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(in) :: dummy_fine, dummy_coarse, weights_high_low, weights_low_high
      TYPE(quadrature_xyoz_type), intent(in) :: qr
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      REAL(KIND=r_def), allocatable :: basis_adspc2_dummy_coarse_qr(:,:,:,:)
      INTEGER(KIND=i_def) dim_adspc2_dummy_coarse
      REAL(KIND=r_def), pointer :: weights_xy_qr(:) => null(), weights_z_qr(:) => null()
      INTEGER(KIND=i_def) np_xy_qr, np_z_qr
      INTEGER(KIND=i_def) nlayers_dummy_fine
      REAL(KIND=r_def), pointer, dimension(:) :: weights_low_high_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: weights_high_low_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: dummy_coarse_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: dummy_fine_data => null()
      TYPE(field_proxy_type) dummy_fine_proxy, dummy_coarse_proxy, weights_high_low_proxy, weights_low_high_proxy
      TYPE(quadrature_xyoz_proxy_type) qr_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_dummy_fine(:,:) => null(), map_adspc2_dummy_coarse(:,:) => null(), &
&map_adspc3_weights_high_low(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_dummy_fine, undf_adspc1_dummy_fine, ndf_adspc2_dummy_coarse, undf_adspc2_dummy_coarse, &
&ndf_adspc3_weights_high_low, undf_adspc3_weights_high_low
      INTEGER(KIND=i_def) ncell_dummy_fine, ncpc_dummy_fine_dummy_coarse_x, ncpc_dummy_fine_dummy_coarse_y
      INTEGER(KIND=i_def), pointer :: cell_map_dummy_coarse(:,:,:) => null()
      TYPE(mesh_map_type), pointer :: mmap_dummy_fine_dummy_coarse => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_dummy_fine
      TYPE(mesh_type), pointer :: mesh_dummy_fine => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_dummy_coarse
      TYPE(mesh_type), pointer :: mesh_dummy_coarse => null()
      !
      ! Initialise field and/or operator proxies
      !
      dummy_fine_proxy = dummy_fine%get_proxy()
      dummy_fine_data => dummy_fine_proxy%data
      dummy_coarse_proxy = dummy_coarse%get_proxy()
      dummy_coarse_data => dummy_coarse_proxy%data
      weights_high_low_proxy = weights_high_low%get_proxy()
      weights_high_low_data => weights_high_low_proxy%data
      weights_low_high_proxy = weights_low_high%get_proxy()
      weights_low_high_data => weights_low_high_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_dummy_fine = dummy_fine_proxy%vspace%get_nlayers()
      !
      ! Look-up mesh objects and loop limits for inter-grid kernels
      !
      mesh_dummy_fine => dummy_fine_proxy%vspace%get_mesh()
      max_halo_depth_mesh_dummy_fine = mesh_dummy_fine%get_halo_depth()
      mesh_dummy_coarse => dummy_coarse_proxy%vspace%get_mesh()
      max_halo_depth_mesh_dummy_coarse = mesh_dummy_coarse%get_halo_depth()
      mmap_dummy_fine_dummy_coarse => mesh_dummy_coarse%get_mesh_map(mesh_dummy_fine)
      cell_map_dummy_coarse => mmap_dummy_fine_dummy_coarse%get_whole_cell_map()
      ncell_dummy_fine = mesh_dummy_fine%get_last_halo_cell(depth=2)
      ncpc_dummy_fine_dummy_coarse_x = mmap_dummy_fine_dummy_coarse%get_ntarget_cells_per_source_x()
      ncpc_dummy_fine_dummy_coarse_y = mmap_dummy_fine_dummy_coarse%get_ntarget_cells_per_source_y()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_dummy_fine => dummy_fine_proxy%vspace%get_whole_dofmap()
      map_adspc2_dummy_coarse => dummy_coarse_proxy%vspace%get_whole_dofmap()
      map_adspc3_weights_high_low => weights_high_low_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_dummy_fine
      !
      ndf_adspc1_dummy_fine = dummy_fine_proxy%vspace%get_ndf()
      undf_adspc1_dummy_fine = dummy_fine_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_dummy_coarse
      !
      ndf_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_ndf()
      undf_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_weights_high_low
      !
      ndf_adspc3_weights_high_low = weights_high_low_proxy%vspace%get_ndf()
      undf_adspc3_weights_high_low = weights_high_low_proxy%vspace%get_undf()
      !
      ! Look-up quadrature variables
      !
      qr_proxy = qr%get_quadrature_proxy()
      np_xy_qr = qr_proxy%np_xy
      np_z_qr = qr_proxy%np_z
      weights_xy_qr => qr_proxy%weights_xy
      weights_z_qr => qr_proxy%weights_z
      !
      ! Allocate basis/diff-basis arrays
      !
      dim_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_dim_space()
      ALLOCATE (basis_adspc2_dummy_coarse_qr(dim_adspc2_dummy_coarse, ndf_adspc2_dummy_coarse, np_xy_qr, np_z_qr))
      !
      ! Compute basis/diff-basis arrays
      !
      CALL qr%compute_function(BASIS, dummy_coarse_proxy%vspace, dim_adspc2_dummy_coarse, ndf_adspc2_dummy_coarse, &
&basis_adspc2_dummy_coarse_qr)
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh_dummy_coarse%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL weights_scalar_inter_element_order_code(nlayers_dummy_fine, cell_map_dummy_coarse(:,:,cell), &
&ncpc_dummy_fine_dummy_coarse_x, ncpc_dummy_fine_dummy_coarse_y, ncell_dummy_fine, dummy_fine_data, dummy_coarse_data, &
&weights_high_low_data, weights_low_high_data, ndf_adspc1_dummy_fine, undf_adspc1_dummy_fine, map_adspc1_dummy_fine, &
&ndf_adspc2_dummy_coarse, undf_adspc2_dummy_coarse, map_adspc2_dummy_coarse(:,cell), basis_adspc2_dummy_coarse_qr, &
&undf_adspc3_weights_high_low, map_adspc3_weights_high_low(:,cell), np_xy_qr, np_z_qr, weights_xy_qr, weights_z_qr)
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL weights_high_low_proxy%set_dirty()
      CALL weights_low_high_proxy%set_dirty()
      !
      !
      ! Deallocate basis arrays
      !
      DEALLOCATE (basis_adspc2_dummy_coarse_qr)
      !
    END SUBROUTINE invoke_weights_scalar_inter_element_order_kernel_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_source of the source field of the intermesh kernel
  ! alongside that of that of the target field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_map_scalar_fe_to_fv_kernel_type(fine_field, coarse_field, weights_high_to_low)
      USE sci_map_scalar_fe_to_fv_kernel_mod, ONLY: map_scalar_fe_to_fv_code
      USE mesh_map_mod, ONLY: mesh_map_type
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(in) :: fine_field, coarse_field, weights_high_to_low
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_fine_field
      REAL(KIND=r_def), pointer, dimension(:) :: weights_high_to_low_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: coarse_field_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: fine_field_data => null()
      TYPE(field_proxy_type) fine_field_proxy, coarse_field_proxy, weights_high_to_low_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_fine_field(:,:) => null(), map_adspc2_coarse_field(:,:) => null(), &
&map_adspc3_weights_high_to_low(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_fine_field, undf_adspc1_fine_field, ndf_adspc2_coarse_field, undf_adspc2_coarse_field, &
&ndf_adspc3_weights_high_to_low, undf_adspc3_weights_high_to_low
      INTEGER(KIND=i_def) ncell_fine_field, ncpc_fine_field_coarse_field_x, ncpc_fine_field_coarse_field_y
      INTEGER(KIND=i_def), pointer :: cell_map_coarse_field(:,:,:) => null()
      TYPE(mesh_map_type), pointer :: mmap_fine_field_coarse_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_fine_field
      TYPE(mesh_type), pointer :: mesh_fine_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_coarse_field
      TYPE(mesh_type), pointer :: mesh_coarse_field => null()
      !
      ! Initialise field and/or operator proxies
      !
      fine_field_proxy = fine_field%get_proxy()
      fine_field_data => fine_field_proxy%data
      coarse_field_proxy = coarse_field%get_proxy()
      coarse_field_data => coarse_field_proxy%data
      weights_high_to_low_proxy = weights_high_to_low%get_proxy()
      weights_high_to_low_data => weights_high_to_low_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_fine_field = fine_field_proxy%vspace%get_nlayers()
      !
      ! Look-up mesh objects and loop limits for inter-grid kernels
      !
      mesh_fine_field => fine_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_fine_field = mesh_fine_field%get_halo_depth()
      mesh_coarse_field => coarse_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_coarse_field = mesh_coarse_field%get_halo_depth()
      mmap_fine_field_coarse_field => mesh_coarse_field%get_mesh_map(mesh_fine_field)
      cell_map_coarse_field => mmap_fine_field_coarse_field%get_whole_cell_map()
      ncell_fine_field = mesh_fine_field%get_last_halo_cell(depth=2)
      ncpc_fine_field_coarse_field_x = mmap_fine_field_coarse_field%get_ntarget_cells_per_source_x()
      ncpc_fine_field_coarse_field_y = mmap_fine_field_coarse_field%get_ntarget_cells_per_source_y()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_fine_field => fine_field_proxy%vspace%get_whole_dofmap()
      map_adspc2_coarse_field => coarse_field_proxy%vspace%get_whole_dofmap()
      map_adspc3_weights_high_to_low => weights_high_to_low_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_fine_field
      !
      ndf_adspc1_fine_field = fine_field_proxy%vspace%get_ndf()
      undf_adspc1_fine_field = fine_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_coarse_field
      !
      ndf_adspc2_coarse_field = coarse_field_proxy%vspace%get_ndf()
      undf_adspc2_coarse_field = coarse_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_weights_high_to_low
      !
      ndf_adspc3_weights_high_to_low = weights_high_to_low_proxy%vspace%get_ndf()
      undf_adspc3_weights_high_to_low = weights_high_to_low_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh_coarse_field%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL map_scalar_fe_to_fv_code(nlayers_fine_field, cell_map_coarse_field(:,:,cell), &
&ncpc_fine_field_coarse_field_x, ncpc_fine_field_coarse_field_y, ncell_fine_field, &
&fine_field_data, coarse_field_data, weights_high_to_low_data, ndf_adspc1_fine_field, &
&undf_adspc1_fine_field, map_adspc1_fine_field,ndf_adspc2_coarse_field, &
&undf_adspc2_coarse_field, map_adspc2_coarse_field(:,cell), undf_adspc3_weights_high_to_low, &
&map_adspc3_weights_high_to_low(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL fine_field_proxy%set_dirty()
      !
      !
  END SUBROUTINE invoke_map_scalar_fe_to_fv_kernel_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_source of the source field of the intermesh kernel
  ! alongside that of that of the target field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_map_scalar_fv_to_fe_kernel_type(coarse_field, fine_field, weights_low_to_high)
      USE sci_map_scalar_fv_to_fe_kernel_mod, ONLY: map_scalar_fv_to_fe_code
      USE mesh_map_mod, ONLY: mesh_map_type
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(in) :: coarse_field, fine_field, weights_low_to_high
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_coarse_field
      REAL(KIND=r_def), pointer, dimension(:) :: weights_low_to_high_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: fine_field_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: coarse_field_data => null()
      TYPE(field_proxy_type) coarse_field_proxy, fine_field_proxy, weights_low_to_high_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_coarse_field(:,:) => null(), map_adspc2_fine_field(:,:) => null(), &
&map_adspc3_weights_low_to_high(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_coarse_field, undf_adspc1_coarse_field, ndf_adspc2_fine_field, undf_adspc2_fine_field, &
&ndf_adspc3_weights_low_to_high, undf_adspc3_weights_low_to_high
      INTEGER(KIND=i_def) ncell_fine_field, ncpc_fine_field_coarse_field_x, ncpc_fine_field_coarse_field_y
      INTEGER(KIND=i_def), pointer :: cell_map_coarse_field(:,:,:) => null()
      TYPE(mesh_map_type), pointer :: mmap_fine_field_coarse_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_fine_field
      TYPE(mesh_type), pointer :: mesh_fine_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_coarse_field
      TYPE(mesh_type), pointer :: mesh_coarse_field => null()
      !
      ! Initialise field and/or operator proxies
      !
      coarse_field_proxy = coarse_field%get_proxy()
      coarse_field_data => coarse_field_proxy%data
      fine_field_proxy = fine_field%get_proxy()
      fine_field_data => fine_field_proxy%data
      weights_low_to_high_proxy = weights_low_to_high%get_proxy()
      weights_low_to_high_data => weights_low_to_high_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_coarse_field = coarse_field_proxy%vspace%get_nlayers()
      !
      ! Look-up mesh objects and loop limits for inter-grid kernels
      !
      mesh_fine_field => fine_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_fine_field = mesh_fine_field%get_halo_depth()
      mesh_coarse_field => coarse_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_coarse_field = mesh_coarse_field%get_halo_depth()
      mmap_fine_field_coarse_field => mesh_coarse_field%get_mesh_map(mesh_fine_field)
      cell_map_coarse_field => mmap_fine_field_coarse_field%get_whole_cell_map()
      ncell_fine_field = mesh_fine_field%get_last_halo_cell(depth=2)
      ncpc_fine_field_coarse_field_x = mmap_fine_field_coarse_field%get_ntarget_cells_per_source_x()
      ncpc_fine_field_coarse_field_y = mmap_fine_field_coarse_field%get_ntarget_cells_per_source_y()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_coarse_field => coarse_field_proxy%vspace%get_whole_dofmap()
      map_adspc2_fine_field => fine_field_proxy%vspace%get_whole_dofmap()
      map_adspc3_weights_low_to_high => weights_low_to_high_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_coarse_field
      !
      ndf_adspc1_coarse_field = coarse_field_proxy%vspace%get_ndf()
      undf_adspc1_coarse_field = coarse_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_fine_field
      !
      ndf_adspc2_fine_field = fine_field_proxy%vspace%get_ndf()
      undf_adspc2_fine_field = fine_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_weights_low_to_high
      !
      ndf_adspc3_weights_low_to_high = weights_low_to_high_proxy%vspace%get_ndf()
      undf_adspc3_weights_low_to_high = weights_low_to_high_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh_coarse_field%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL map_scalar_fv_to_fe_code(nlayers_coarse_field, cell_map_coarse_field(:,:,cell), &
&ncpc_fine_field_coarse_field_x, ncpc_fine_field_coarse_field_y, ncell_fine_field, &
&coarse_field_data, fine_field_data, weights_low_to_high_data, ndf_adspc1_coarse_field, &
&undf_adspc1_coarse_field, map_adspc1_coarse_field(:,cell), ndf_adspc2_fine_field, &
&undf_adspc2_fine_field, map_adspc2_fine_field, undf_adspc3_weights_low_to_high, &
&map_adspc3_weights_low_to_high(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL coarse_field_proxy%set_dirty()
      !
      !
  END SUBROUTINE invoke_map_scalar_fv_to_fe_kernel_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_source of the source field of the intermesh kernel
  ! alongside that of that of the target field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_weights_w2_inter_element_order_kernel_type(dummy_fine, dummy_coarse, weights_high_low, weights_low_high, qr)
    USE sci_weights_w2_inter_element_order_kernel_mod, ONLY: weights_w2_inter_element_order_code
    USE quadrature_xyoz_mod, ONLY: quadrature_xyoz_type, quadrature_xyoz_proxy_type
    USE function_space_mod, ONLY: BASIS, DIFF_BASIS
    USE mesh_map_mod, ONLY: mesh_map_type
    USE mesh_mod, ONLY: mesh_type
    TYPE(field_type), intent(in) :: dummy_fine, dummy_coarse, weights_high_low, weights_low_high
    TYPE(quadrature_xyoz_type), intent(in) :: qr
    INTEGER(KIND=i_def) cell
    INTEGER(KIND=i_def) loop0_start, loop0_stop
    REAL(KIND=r_def), allocatable :: basis_adspc2_dummy_coarse_qr(:,:,:,:)
    INTEGER(KIND=i_def) dim_adspc2_dummy_coarse
    REAL(KIND=r_def), pointer :: weights_xy_qr(:) => null(), weights_z_qr(:) => null()
    INTEGER(KIND=i_def) np_xy_qr, np_z_qr
    INTEGER(KIND=i_def) nlayers_dummy_fine
    REAL(KIND=r_def), pointer, dimension(:) :: weights_low_high_data => null()
    REAL(KIND=r_def), pointer, dimension(:) :: weights_high_low_data => null()
    REAL(KIND=r_def), pointer, dimension(:) :: dummy_coarse_data => null()
    REAL(KIND=r_def), pointer, dimension(:) :: dummy_fine_data => null()
    TYPE(field_proxy_type) dummy_fine_proxy, dummy_coarse_proxy, weights_high_low_proxy, weights_low_high_proxy
    TYPE(quadrature_xyoz_proxy_type) qr_proxy
    INTEGER(KIND=i_def), pointer :: map_adspc1_dummy_fine(:,:) => null(), map_adspc2_dummy_coarse(:,:) => null(), &
&map_adspc3_weights_high_low(:,:) => null()
    INTEGER(KIND=i_def) ndf_adspc1_dummy_fine, undf_adspc1_dummy_fine, ndf_adspc2_dummy_coarse, undf_adspc2_dummy_coarse, &
&ndf_adspc3_weights_high_low, undf_adspc3_weights_high_low
    INTEGER(KIND=i_def) ncell_dummy_fine, ncpc_dummy_fine_dummy_coarse_x, ncpc_dummy_fine_dummy_coarse_y
    INTEGER(KIND=i_def), pointer :: cell_map_dummy_coarse(:,:,:) => null()
    TYPE(mesh_map_type), pointer :: mmap_dummy_fine_dummy_coarse => null()
    INTEGER(KIND=i_def) max_halo_depth_mesh_dummy_fine
    TYPE(mesh_type), pointer :: mesh_dummy_fine => null()
    INTEGER(KIND=i_def) max_halo_depth_mesh_dummy_coarse
    TYPE(mesh_type), pointer :: mesh_dummy_coarse => null()
    !
    ! Initialise field and/or operator proxies
    !
    dummy_fine_proxy = dummy_fine%get_proxy()
    dummy_fine_data => dummy_fine_proxy%data
    dummy_coarse_proxy = dummy_coarse%get_proxy()
    dummy_coarse_data => dummy_coarse_proxy%data
    weights_high_low_proxy = weights_high_low%get_proxy()
    weights_high_low_data => weights_high_low_proxy%data
    weights_low_high_proxy = weights_low_high%get_proxy()
    weights_low_high_data => weights_low_high_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_dummy_fine = dummy_fine_proxy%vspace%get_nlayers()
    !
    ! Look-up mesh objects and loop limits for inter-grid kernels
    !
    mesh_dummy_fine => dummy_fine_proxy%vspace%get_mesh()
    max_halo_depth_mesh_dummy_fine = mesh_dummy_fine%get_halo_depth()
    mesh_dummy_coarse => dummy_coarse_proxy%vspace%get_mesh()
    max_halo_depth_mesh_dummy_coarse = mesh_dummy_coarse%get_halo_depth()
    mmap_dummy_fine_dummy_coarse => mesh_dummy_coarse%get_mesh_map(mesh_dummy_fine)
    cell_map_dummy_coarse => mmap_dummy_fine_dummy_coarse%get_whole_cell_map()
    ncell_dummy_fine = mesh_dummy_fine%get_last_halo_cell(depth=2)
    ncpc_dummy_fine_dummy_coarse_x = mmap_dummy_fine_dummy_coarse%get_ntarget_cells_per_source_x()
    ncpc_dummy_fine_dummy_coarse_y = mmap_dummy_fine_dummy_coarse%get_ntarget_cells_per_source_y()
    !
    ! Look-up dofmaps for each function space
    !
    map_adspc1_dummy_fine => dummy_fine_proxy%vspace%get_whole_dofmap()
    map_adspc2_dummy_coarse => dummy_coarse_proxy%vspace%get_whole_dofmap()
    map_adspc3_weights_high_low => weights_high_low_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for adspc1_dummy_fine
    !
    ndf_adspc1_dummy_fine = dummy_fine_proxy%vspace%get_ndf()
    undf_adspc1_dummy_fine = dummy_fine_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_dummy_coarse
    !
    ndf_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_ndf()
    undf_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_weights_high_low
    !
    ndf_adspc3_weights_high_low = weights_high_low_proxy%vspace%get_ndf()
    undf_adspc3_weights_high_low = weights_high_low_proxy%vspace%get_undf()
    !
    ! Look-up quadrature variables
    !
    qr_proxy = qr%get_quadrature_proxy()
    np_xy_qr = qr_proxy%np_xy
    np_z_qr = qr_proxy%np_z
    weights_xy_qr => qr_proxy%weights_xy
    weights_z_qr => qr_proxy%weights_z
    !
    ! Allocate basis/diff-basis arrays
    !
    dim_adspc2_dummy_coarse = dummy_coarse_proxy%vspace%get_dim_space()
    ALLOCATE (basis_adspc2_dummy_coarse_qr(dim_adspc2_dummy_coarse, ndf_adspc2_dummy_coarse, np_xy_qr, np_z_qr))
    !
    ! Compute basis/diff-basis arrays
    !
    CALL qr%compute_function(BASIS, dummy_coarse_proxy%vspace, dim_adspc2_dummy_coarse, ndf_adspc2_dummy_coarse, &
&basis_adspc2_dummy_coarse_qr)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh_dummy_coarse%get_last_edge_cell()
    !
    ! Call kernels and communication routines
    !
    DO cell = loop0_start, loop0_stop, 1
      CALL weights_w2_inter_element_order_code(nlayers_dummy_fine, cell_map_dummy_coarse(:,:,cell), &
&ncpc_dummy_fine_dummy_coarse_x, ncpc_dummy_fine_dummy_coarse_y, ncell_dummy_fine, dummy_fine_data, &
&dummy_coarse_data, weights_high_low_data, weights_low_high_data, ndf_adspc1_dummy_fine, &
&undf_adspc1_dummy_fine, map_adspc1_dummy_fine, ndf_adspc2_dummy_coarse, undf_adspc2_dummy_coarse, &
&map_adspc2_dummy_coarse(:,cell), basis_adspc2_dummy_coarse_qr, undf_adspc3_weights_high_low, &
&map_adspc3_weights_high_low(:,cell), np_xy_qr, np_z_qr, weights_xy_qr, weights_z_qr)
    END DO
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    CALL weights_high_low_proxy%set_dirty()
    CALL weights_low_high_proxy%set_dirty()
    !
    !
    ! Deallocate basis arrays
    !
    DEALLOCATE (basis_adspc2_dummy_coarse_qr)
    !
  END SUBROUTINE invoke_weights_w2_inter_element_order_kernel_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_target of the target field of the intermesh kernel
  ! alongside that of that of the target field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_map_w2_fv_to_fe_kernel_type(target_field, source_field, weights, face_selector_ew, face_selector_ns)
      USE sci_map_w2_fv_to_fe_kernel_mod, ONLY: map_w2_fv_to_fe_code
      USE mesh_map_mod, ONLY: mesh_map_type
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(in) :: target_field, source_field, weights
      TYPE(integer_field_type), intent(in) :: face_selector_ew, face_selector_ns
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_target_field
      INTEGER(KIND=i_def), pointer, dimension(:) :: face_selector_ns_data => null()
      INTEGER(KIND=i_def), pointer, dimension(:) :: face_selector_ew_data => null()
      TYPE(integer_field_proxy_type) face_selector_ew_proxy, face_selector_ns_proxy
      REAL(KIND=r_def), pointer, dimension(:) :: weights_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: source_field_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: target_field_data => null()
      TYPE(field_proxy_type) target_field_proxy, source_field_proxy, weights_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null(), &
&map_adspc3_weights(:,:) => null(), map_adspc4_face_selector_ew(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field, &
&ndf_adspc3_weights, undf_adspc3_weights, ndf_adspc4_face_selector_ew, undf_adspc4_face_selector_ew
      INTEGER(KIND=i_def) ncell_source_field, ncpc_source_field_target_field_x, ncpc_source_field_target_field_y
      INTEGER(KIND=i_def), pointer :: cell_map_target_field(:,:,:) => null()
      TYPE(mesh_map_type), pointer :: mmap_source_field_target_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_target_field
      TYPE(mesh_type), pointer :: mesh_target_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_source_field
      TYPE(mesh_type), pointer :: mesh_source_field => null()
      !
      ! Initialise field and/or operator proxies
      !
      target_field_proxy = target_field%get_proxy()
      target_field_data => target_field_proxy%data
      source_field_proxy = source_field%get_proxy()
      source_field_data => source_field_proxy%data
      weights_proxy = weights%get_proxy()
      weights_data => weights_proxy%data
      face_selector_ew_proxy = face_selector_ew%get_proxy()
      face_selector_ew_data => face_selector_ew_proxy%data
      face_selector_ns_proxy = face_selector_ns%get_proxy()
      face_selector_ns_data => face_selector_ns_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_target_field = target_field_proxy%vspace%get_nlayers()
      !
      ! Look-up mesh objects and loop limits for inter-grid kernels
      !
      mesh_source_field => source_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_source_field = mesh_source_field%get_halo_depth()
      mesh_target_field => target_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_target_field = mesh_target_field%get_halo_depth()
      mmap_source_field_target_field => mesh_target_field%get_mesh_map(mesh_source_field)
      cell_map_target_field => mmap_source_field_target_field%get_whole_cell_map()
      ncell_source_field = mesh_source_field%get_last_halo_cell(depth=2)
      ncpc_source_field_target_field_x = mmap_source_field_target_field%get_ntarget_cells_per_source_x()
      ncpc_source_field_target_field_y = mmap_source_field_target_field%get_ntarget_cells_per_source_y()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
      map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
      map_adspc3_weights => weights_proxy%vspace%get_whole_dofmap()
      map_adspc4_face_selector_ew => face_selector_ew_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_target_field
      !
      ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
      undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_source_field
      !
      ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
      undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_weights
      !
      ndf_adspc3_weights = weights_proxy%vspace%get_ndf()
      undf_adspc3_weights = weights_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc4_face_selector_ew
      !
      ndf_adspc4_face_selector_ew = face_selector_ew_proxy%vspace%get_ndf()
      undf_adspc4_face_selector_ew = face_selector_ew_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh_target_field%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL map_w2_fv_to_fe_code(nlayers_target_field, cell_map_target_field(:,:,cell), ncpc_source_field_target_field_x, &
&ncpc_source_field_target_field_y, ncell_source_field, target_field_data, source_field_data, weights_data, face_selector_ew_data, &
&face_selector_ns_data, ndf_adspc1_target_field, undf_adspc1_target_field, map_adspc1_target_field(:,cell), &
&ndf_adspc2_source_field, undf_adspc2_source_field, map_adspc2_source_field, undf_adspc3_weights, map_adspc3_weights(:,cell), &
&undf_adspc4_face_selector_ew, map_adspc4_face_selector_ew(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL target_field_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_map_w2_fv_to_fe_kernel_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Passes in ndf_source of the source field of the intermesh kernel
  ! alongside that of that of the target field. Inter-grid kernels are not
  ! currently allowed GH_SCALAR arguments or this would not be necessary
  ! See PSyclone issues #2504 and #868
  SUBROUTINE invoke_map_w2_fe_to_fv_kernel_type(target_field, source_field, weights, face_selector_ew, face_selector_ns)
      USE sci_map_w2_fe_to_fv_kernel_mod, ONLY: map_w2_fe_to_fv_code
      USE mesh_map_mod, ONLY: mesh_map_type
      USE mesh_mod, ONLY: mesh_type
      TYPE(field_type), intent(in) :: target_field, source_field, weights
      TYPE(integer_field_type), intent(in) :: face_selector_ew, face_selector_ns
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_target_field
      INTEGER(KIND=i_def), pointer, dimension(:) :: face_selector_ns_data => null()
      INTEGER(KIND=i_def), pointer, dimension(:) :: face_selector_ew_data => null()
      TYPE(integer_field_proxy_type) face_selector_ew_proxy, face_selector_ns_proxy
      REAL(KIND=r_def), pointer, dimension(:) :: weights_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: source_field_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: target_field_data => null()
      TYPE(field_proxy_type) target_field_proxy, source_field_proxy, weights_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_target_field(:,:) => null(), map_adspc2_source_field(:,:) => null(), &
&map_adspc3_weights(:,:) => null(), map_adspc4_face_selector_ew(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_target_field, undf_adspc1_target_field, ndf_adspc2_source_field, undf_adspc2_source_field, &
&ndf_adspc3_weights, undf_adspc3_weights, ndf_adspc4_face_selector_ew, undf_adspc4_face_selector_ew
      INTEGER(KIND=i_def) ncell_target_field, ncpc_target_field_source_field_x, ncpc_target_field_source_field_y
      INTEGER(KIND=i_def), pointer :: cell_map_source_field(:,:,:) => null()
      TYPE(mesh_map_type), pointer :: mmap_target_field_source_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_target_field
      TYPE(mesh_type), pointer :: mesh_target_field => null()
      INTEGER(KIND=i_def) max_halo_depth_mesh_source_field
      TYPE(mesh_type), pointer :: mesh_source_field => null()
      !
      ! Initialise field and/or operator proxies
      !
      target_field_proxy = target_field%get_proxy()
      target_field_data => target_field_proxy%data
      source_field_proxy = source_field%get_proxy()
      source_field_data => source_field_proxy%data
      weights_proxy = weights%get_proxy()
      weights_data => weights_proxy%data
      face_selector_ew_proxy = face_selector_ew%get_proxy()
      face_selector_ew_data => face_selector_ew_proxy%data
      face_selector_ns_proxy = face_selector_ns%get_proxy()
      face_selector_ns_data => face_selector_ns_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_target_field = target_field_proxy%vspace%get_nlayers()
      !
      ! Look-up mesh objects and loop limits for inter-grid kernels
      !
      mesh_target_field => target_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_target_field = mesh_target_field%get_halo_depth()
      mesh_source_field => source_field_proxy%vspace%get_mesh()
      max_halo_depth_mesh_source_field = mesh_source_field%get_halo_depth()
      mmap_target_field_source_field => mesh_source_field%get_mesh_map(mesh_target_field)
      cell_map_source_field => mmap_target_field_source_field%get_whole_cell_map()
      ncell_target_field = mesh_target_field%get_last_halo_cell(depth=2)
      ncpc_target_field_source_field_x = mmap_target_field_source_field%get_ntarget_cells_per_source_x()
      ncpc_target_field_source_field_y = mmap_target_field_source_field%get_ntarget_cells_per_source_y()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_target_field => target_field_proxy%vspace%get_whole_dofmap()
      map_adspc2_source_field => source_field_proxy%vspace%get_whole_dofmap()
      map_adspc3_weights => weights_proxy%vspace%get_whole_dofmap()
      map_adspc4_face_selector_ew => face_selector_ew_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_target_field
      !
      ndf_adspc1_target_field = target_field_proxy%vspace%get_ndf()
      undf_adspc1_target_field = target_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_source_field
      !
      ndf_adspc2_source_field = source_field_proxy%vspace%get_ndf()
      undf_adspc2_source_field = source_field_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_weights
      !
      ndf_adspc3_weights = weights_proxy%vspace%get_ndf()
      undf_adspc3_weights = weights_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc4_face_selector_ew
      !
      ndf_adspc4_face_selector_ew = face_selector_ew_proxy%vspace%get_ndf()
      undf_adspc4_face_selector_ew = face_selector_ew_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh_source_field%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL map_w2_fe_to_fv_code(nlayers_target_field, cell_map_source_field(:,:,cell), ncpc_target_field_source_field_x, &
&ncpc_target_field_source_field_y, ncell_target_field, target_field_data, source_field_data, weights_data, face_selector_ew_data, &
&face_selector_ns_data, ndf_adspc1_target_field, undf_adspc1_target_field, map_adspc1_target_field, ndf_adspc2_source_field, &
&undf_adspc2_source_field, map_adspc2_source_field(:,cell), undf_adspc3_weights, map_adspc3_weights(:,cell), &
&undf_adspc4_face_selector_ew, map_adspc4_face_selector_ew(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL target_field_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_map_w2_fe_to_fv_kernel_type

  end module sci_psykal_light_mod
