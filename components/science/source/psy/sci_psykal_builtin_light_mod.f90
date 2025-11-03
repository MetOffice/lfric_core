!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module sci_psykal_builtin_light_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, int32
  use constants_mod, only : i_def, i_long, r_def
  use field_mod,     only : field_type, field_proxy_type

  implicit none

  public

contains

  !--------------------------------------------------------------------
  subroutine invoke_pointwise_convert_xyz2llr( coords)
    use coord_transform_mod, only: xyz2llr
    implicit none
    type(field_type), intent(inout) :: coords(3)

    type(field_proxy_type) :: x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: llr(3)

    do i = 1,3
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = x_p(1)%vspace%get_last_dof_annexed()

    do df = 1, undf
      call xyz2llr(x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), &
                   llr(1), llr(2), llr(3))
      x_p(1)%data(df) = llr(1)
      x_p(2)%data(df) = llr(2)
      x_p(3)%data(df) = llr(3)
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

  end subroutine invoke_pointwise_convert_xyz2llr

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #489. See that issue for further
  ! details.
  subroutine invoke_real32_field_min_max(field_min_norm, &
                                         field_max_norm, &
                                         real32_field)

    use scalar_real32_mod,  only: scalar_real32_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use field_real32_mod,   only: field_real32_type, &
                                  field_real32_proxy_type

    implicit none

    real(kind=real32),              intent(out)  :: field_min_norm
    real(kind=real32),              intent(out)  :: field_max_norm
    type(field_real32_type),         intent(in)  :: real32_field
    type(scalar_real32_type)                     :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real32), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real32), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real32)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start
    integer(kind=i_def)                          :: loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_real32_proxy_type)                :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = real32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_real32_field_min_max

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #489. See that issue for further
  !  details.
  subroutine invoke_real64_field_min_max(field_min_norm, &
                                         field_max_norm, &
                                         real64_field)

    use scalar_real64_mod,  only: scalar_real64_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use field_real64_mod,   only: field_real64_type, &
                                  field_real64_proxy_type

    implicit none

    real(kind=real64),               intent(out) :: field_min_norm
    real(kind=real64),               intent(out) :: field_max_norm
    type(field_real64_type),          intent(in) :: real64_field
    type(scalar_real64_type)                     :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real64), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real64), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real64)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start
    integer(kind=i_def)                          :: loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_real64_proxy_type)                :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = real64_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_real64_field_min_max

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #489. See that issue for further details.
  subroutine invoke_int32_field_min_max(field_min_norm, &
                                        field_max_norm, &
                                        int32_field)

    use scalar_int32_mod,   only: scalar_int32_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use field_int32_mod,    only: field_int32_type, &
                                  field_int32_proxy_type

    implicit none

    integer(kind=int32),               intent(out) :: field_min_norm
    integer(kind=int32),               intent(out) :: field_max_norm
    type(field_int32_type),            intent(in)  :: int32_field
    type(scalar_int32_type)                        :: global_min, &
                                                      global_max
    integer(kind=i_def)                            :: df
    integer(kind=int32), allocatable, dimension(:) :: l_field_min_norm
    integer(kind=int32), allocatable, dimension(:) :: l_field_max_norm
    integer(kind=int32)                            :: minv, maxv
    integer(kind=i_def)                            :: th_idx
    integer(kind=i_def)                            :: loop0_start
    integer(kind=i_def)                            :: loop0_stop
    integer(kind=i_def)                            :: nthreads
    type(field_int32_proxy_type)                   :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = int32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)

    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_int32_field_min_max

  !---------------------------------------------------------------------
  subroutine invoke_real32_local_field_min_max(field_min_norm, &
                                               field_max_norm, &
                                               real32_field)

    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_real32_mod,   only: field_real32_type, &
                                  field_real32_proxy_type

    implicit none

    real(kind=real32),              intent(out)  :: field_min_norm
    real(kind=real32),              intent(out)  :: field_max_norm
    type(field_real32_type),         intent(in)  :: real32_field
    integer(kind=i_def)                          :: df
    real(kind=real32), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real32), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real32)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start
    integer(kind=i_def)                          :: loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_real32_proxy_type)                :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = real32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                     field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                     field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx = 2, nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)
    !
  end subroutine invoke_real32_local_field_min_max

  !---------------------------------------------------------------------
  subroutine invoke_real64_local_field_min_max(field_min_norm, &
                                               field_max_norm, &
                                               real64_field)

    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_real64_mod,   only: field_real64_type, &
                                  field_real64_proxy_type

    implicit none

    real(kind=real64),               intent(out) :: field_min_norm
    real(kind=real64),               intent(out) :: field_max_norm
    type(field_real64_type),          intent(in) :: real64_field
    integer(kind=i_def)                          :: df
    real(kind=real64), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real64), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real64)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start
    integer(kind=i_def)                          :: loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_real64_proxy_type)                :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = real64_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                     field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                     field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx = 2, nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)
    !
  end subroutine invoke_real64_local_field_min_max

  !---------------------------------------------------------------------
  subroutine invoke_int32_local_field_min_max(field_min_norm, &
                                              field_max_norm, &
                                              int32_field)

    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_int32_mod,    only: field_int32_type, field_int32_proxy_type

    implicit none

    integer(kind=int32),               intent(out) :: field_min_norm
    integer(kind=int32),               intent(out) :: field_max_norm
    type(field_int32_type),            intent(in)  :: int32_field
    integer(kind=i_def)                            :: df
    integer(kind=int32), allocatable, dimension(:) :: l_field_min_norm
    integer(kind=int32), allocatable, dimension(:) :: l_field_max_norm
    integer(kind=int32)                            :: minv, maxv
    integer(kind=i_def)                            :: th_idx
    integer(kind=i_def)                            :: loop0_start
    integer(kind=i_def)                            :: loop0_stop
    integer(kind=i_def)                            :: nthreads
    type(field_int32_proxy_type)                   :: field_proxy
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = int32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    allocate (l_field_min_norm(nthreads))
    allocate (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                     field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                     field_proxy%data(df))
    end do
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    do th_idx = 2, nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    end do
    deallocate (l_field_min_norm, l_field_max_norm)
    !
  end subroutine invoke_int32_local_field_min_max

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #2674. See that issue for further
  ! details.
  subroutine invoke_copy_field_32_64(fsrce_32, fdest_64)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real32_mod,   only: field_real32_type, &
                                   field_real32_proxy_type
     use field_real64_mod,   only: field_real64_type, &
                                   field_real64_proxy_type

     implicit none

     type(field_real32_type), intent(in)     :: fsrce_32
     type(field_real64_type), intent(inout)  :: fdest_64

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real32_proxy_type)   :: fsrce_32_proxy
     type(field_real64_proxy_type)   :: fdest_64_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_32_proxy = fsrce_32%get_proxy()
     fdest_64_proxy = fdest_64%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_64_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     IF (fsrce_32_proxy%is_dirty(depth=1)) THEN
       ! only copy the owned dofs
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_annexed()
     ELSE
       ! copy the 1st halo row as well
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_halo(1)
     END IF
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     DO df=loop0_start,loop0_stop
       fdest_64_proxy%data(df) = real(fsrce_32_proxy%data(df), real64)
     END DO
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     CALL fdest_64_proxy%set_dirty()
     IF (.not. fsrce_32_proxy%is_dirty(depth=1)) THEN
       CALL fdest_64_proxy%set_clean(1)
     END IF
     !
  end subroutine invoke_copy_field_32_64

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #2674. See that issue for further
  ! details.
  subroutine invoke_copy_field_64_32(fsrce_64, fdest_32)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real32_mod,   only: field_real32_type, &
                                   field_real32_proxy_type
     use field_real64_mod,   only: field_real64_type, &
                                   field_real64_proxy_type

     implicit none

     type(field_real64_type), intent(in)     :: fsrce_64
     type(field_real32_type), intent(inout)  :: fdest_32

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real64_proxy_type)   :: fsrce_64_proxy
     type(field_real32_proxy_type)   :: fdest_32_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_64_proxy = fsrce_64%get_proxy()
     fdest_32_proxy = fdest_32%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_32_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     IF (fsrce_64_proxy%is_dirty(depth=1)) THEN
       ! only copy the owned dofs
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_annexed()
     ELSE
       ! copy the 1st halo row as well
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_halo(1)
     END IF
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     DO df=loop0_start,loop0_stop
       fdest_32_proxy%data(df) = real(fsrce_64_proxy%data(df), real32)
     END DO
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     CALL fdest_32_proxy%set_dirty()
     IF (.not. fsrce_64_proxy%is_dirty(depth=1)) THEN
       CALL fdest_32_proxy%set_clean(1)
     END IF
     !
  end subroutine invoke_copy_field_64_32

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #2674. See that issue for further
  ! details.
  subroutine invoke_copy_field_32_32(fsrce_32, fdest_32)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real32_mod,   only: field_real32_type, &
                                   field_real32_proxy_type

     implicit none

     type(field_real32_type), intent(in)     :: fsrce_32
     type(field_real32_type), intent(inout)  :: fdest_32

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real32_proxy_type)   :: fsrce_32_proxy
     type(field_real32_proxy_type)   :: fdest_32_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_32_proxy = fsrce_32%get_proxy()
     fdest_32_proxy = fdest_32%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_32_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     IF (fsrce_32_proxy%is_dirty(depth=1)) THEN
       ! only copy the owned dofs
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_annexed()
     ELSE
       ! copy the 1st halo row as well
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_halo(1)
     END IF
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     DO df=loop0_start,loop0_stop
       fdest_32_proxy%data(df) = real(fsrce_32_proxy%data(df), real32)
     END DO
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     CALL fdest_32_proxy%set_dirty()
     IF (.not. fsrce_32_proxy%is_dirty(depth=1)) THEN
       CALL fdest_32_proxy%set_clean(1)
     END IF
     !
  end subroutine invoke_copy_field_32_32

  !---------------------------------------------------------------------
  ! This is a PSyKAl-lite implementation of a built-in that will be
  ! implemented under PSYclone issue #2674. See that issue for further
  ! details.
  subroutine invoke_copy_field_64_64(fsrce_64, fdest_64)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real64_mod,   only: field_real64_type, &
                                   field_real64_proxy_type

     implicit none

     type(field_real64_type), intent(in)     :: fsrce_64
     type(field_real64_type), intent(inout)  :: fdest_64

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real64_proxy_type)   :: fsrce_64_proxy
     type(field_real64_proxy_type)   :: fdest_64_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_64_proxy = fsrce_64%get_proxy()
     fdest_64_proxy = fdest_64%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_64_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     IF (fsrce_64_proxy%is_dirty(depth=1)) THEN
       ! only copy the owned dofs
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_annexed()
     ELSE
       ! copy the 1st halo row as well
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_halo(1)
     END IF
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     DO df=loop0_start,loop0_stop
       fdest_64_proxy%data(df) = real(fsrce_64_proxy%data(df), real64)
     END DO
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     CALL fdest_64_proxy%set_dirty()
     IF (.not. fsrce_64_proxy%is_dirty(depth=1)) THEN
       CALL fdest_64_proxy%set_clean(1)
     END IF
     !
  end subroutine invoke_copy_field_64_64

end module sci_psykal_builtin_light_mod
