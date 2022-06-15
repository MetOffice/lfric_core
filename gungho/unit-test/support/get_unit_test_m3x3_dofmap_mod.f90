!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_m3x3_dofmap_mod
! A module containing a collection of helper routines that provide canned
! dofmaps (and stencil dofmaps) on a simple 3x3 mesh for use when writing
! unit tests of kernels

  use constants_mod, only : i_def

  implicit none

  private

  public :: get_w0_m3x3_dofmap,            &
            get_w1_m3x3_dofmap,            &
            get_w2_m3x3_dofmap,            &
            get_w2v_m3x3_dofmap,           &
            get_w2h_m3x3_dofmap,           &
            get_w2broken_m3x3_dofmap,      &
            get_w2trace_m3x3_dofmap,       &
            get_w3_m3x3_dofmap,            &
            get_wtheta_m3x3_dofmap,        &
            get_wchi_m3x3_dofmap,          &
            get_m3x3_stencil_dofmap_point, &
            get_m3x3_stencil_dofmap_cross, &
            get_m3x3_stencil_dofmap_region,&
            get_m3x3_stencil_dofmap_1dx,   &
            get_m3x3_stencil_dofmap_1dy
  contains

!---------------------------------------------------------------------

  subroutine get_w0_m3x3_dofmap(map_w0, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs on vertices
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w0(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: i,j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w0(8,9))
    map_w0 = reshape( [0,1,2,3, 0,0,0,0, &
                       1,4,5,2, 0,0,0,0, &
                       4,0,3,5, 0,0,0,0, &
                       3,2,6,7, 0,0,0,0, &
                       2,5,8,6, 0,0,0,0, &
                       5,3,7,8, 0,0,0,0, &
                       7,6,1,0, 0,0,0,0, &
                       6,8,4,1, 0,0,0,0, &
                       8,7,0,4, 0,0,0,0], [8,9] )
    do j=1,9
      do i=1,4
        map_w0(i+4,j) = map_w0(i,j)*(nlay+1)+2
        map_w0(i,j)   = map_w0(i,j)*(nlay+1)+1
      end do
    end do

  end subroutine get_w0_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w1_m3x3_dofmap(map_w1, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs on edges
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w1(:,:)
    integer(i_def), optional,    intent(in ) :: nlayers

    integer(i_def), parameter ::  ncol_edge(27) = &
               [1,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,0,1,1,1,1]
    integer(i_def) :: nlay
    integer(i_def) :: i,j,k,sum

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w1(12,9))
    map_w1 = reshape( [ 0, 1, 2, 3,   4, 5, 6, 7,  0,0,0,0, &
                        2, 8, 9,10,   5,11,12, 6,  0,0,0,0, &
                        9,13, 0,14,  11, 4, 7,12,  0,0,0,0, &
                       15, 3,16,17,   7, 6,18,19,  0,0,0,0, &
                       16,10,20,21,   6,12,22,18,  0,0,0,0, &
                       20,14,15,23,  12, 7,19,22,  0,0,0,0, &
                       24,17,25, 1,  19,18, 5, 4,  0,0,0,0, &
                       25,21,26, 8,  18,22,11, 5,  0,0,0,0, &
                       26,23,24,13,  22,19, 4,11,  0,0,0,0], [12,9] )
    do j=1,9
      do i=1,4
        sum=1
        do k=1,map_w1(i,j)
          sum=sum+nlay+ncol_edge(k)
        end do
        map_w1(i+8,j) = sum+1
        map_w1(i,j)   = sum
        sum=1
        do k=1,map_w1(i+4,j)
          sum=sum+nlay+ncol_edge(k)
        end do
        map_w1(i+4,j) = sum
      end do
    end do

  end subroutine get_w1_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w2_m3x3_dofmap(map_w2, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs on faces
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w2(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def), parameter ::  ncol_face(27) = &
               [0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,1,1]
    integer(i_def) :: nlay
    integer(i_def) :: i,j,k,sum

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w2(6,9))
    map_w2 = reshape( [ 0, 1, 2, 3,   4,  0, &
                        2, 5, 6, 7,   8,  0, &
                        6, 9, 0,10,  11,  0, &
                       12, 3,13,14,  15,  0, &
                       13, 7,16,17,  18,  0, &
                       16,10,12,19,  20,  0, &
                       21,14,22, 1,  23,  0, &
                       22,17,24, 5,  25,  0, &
                       24,19,21, 9,  26,  0], [6,9] )
    do j=1,9
      do i=1,4
        sum=1
        do k=1,map_w2(i,j)
          sum=sum+nlay+ncol_face(k)
        end do
        map_w2(i,j) = sum
      end do
      sum=1
      do k=1,map_w2(5,j)
        sum=sum+nlay+ncol_face(k)
      end do
      map_w2(5,j) = sum
      map_w2(6,j) = sum+1
    end do

  end subroutine get_w2_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w2v_m3x3_dofmap(map_w2v, nlayers)
    ! Generate a per-cell dofmap for a 3x3 mesh with dofs on top/bottom faces
    ! The dofmap for lowest order w2v is the same a for wtheta - so use that
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w2v(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    call get_wtheta_m3x3_dofmap(map_w2v, nlayers)

  end subroutine get_w2v_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w2h_m3x3_dofmap(map_w2h, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs on the side faces
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w2h(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: i,j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w2h(4,9))
    map_w2h = reshape( [ 0, 1, 2, 3,  &
                         2, 4, 5, 6,  &
                         5, 7, 0, 8,  &
                         9, 3,10,11,  &
                        10, 6,12,13,  &
                        12, 8, 9,14,  &
                        15,11,16, 1,  &
                        16,13,17, 4,  &
                        17,14,15, 7], [4,9] )
    do j=1,9
      do i=1,4
        map_w2h(i,j)   = map_w2h(i,j)*nlay+1
      end do
    end do

  end subroutine get_w2h_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w2broken_m3x3_dofmap(map_w2broken, nlayers)
    ! Create a per-cell dofmap for a 3x3 mesh with dofs in the cell volume
    ! but at the locations of the face dofs
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w2broken(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: i,j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w2broken(6,9))
    do j=1,9
      do i=1,6
        map_w2broken(i,j) = (j-1)*6*nlay+(i-1)*nlay+1
      end do
    end do

  end subroutine get_w2broken_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w2trace_m3x3_dofmap(map_w2trace, nlayers)
    ! Generate a per-cell dofmap for a 3x3 mesh with dofs for scalars, but
    ! on the faces. This dofmap at lowest order is the same as for w2
    ! - so use that
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w2trace(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    call get_w2_m3x3_dofmap(map_w2trace, nlayers)

  end subroutine get_w2trace_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_w3_m3x3_dofmap(map_w3, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with a dof in the cell volume
    implicit none
    integer(i_def), allocatable, intent(out) :: map_w3(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_w3(1,9))
    do j=1,9
      map_w3(1,j) = (j-1)*nlay+1
    end do

  end subroutine get_w3_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_wtheta_m3x3_dofmap(map_wtheta, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs on top/bottom faces
    implicit none
    integer(i_def), allocatable, intent(out) :: map_wtheta(:,:)
    integer(i_def), optional,    intent(in)  :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_wtheta(2,9))
    do j=1,9
      map_wtheta(1,j) = (j-1)*(nlay+1)+1
      map_wtheta(2,j) = (j-1)*(nlay+1)+2
    end do

  end subroutine get_wtheta_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_wchi_m3x3_dofmap(map_wchi, nlayers)
    ! Calculate a per-cell dofmap for a 3x3 mesh with dofs in the cell volume,
    ! but located at the cell vertices
    implicit none
    integer(i_def), allocatable, intent(out) :: map_wchi(:,:)
    integer(i_def), optional,   intent(in)   :: nlayers

    integer(i_def) :: nlay
    integer(i_def) :: i,j

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    allocate(map_wchi(8,9))
    do j=1,9
      do i=1,8
        map_wchi(i,j) = (j-1)*8*nlay+(i-1)*nlay+1
      end do
    end do

  end subroutine get_wchi_m3x3_dofmap

!---------------------------------------------------------------------

  subroutine get_m3x3_stencil_dofmap_point(stencil_map, dofmap)
    ! Generate a point stencil dofmap for any given input 3x3 dofmap
    implicit none
    integer(i_def), allocatable, intent(out) :: stencil_map(:,:,:)
    integer(i_def),              intent(in)  :: dofmap(:,:)

    integer(i_def) :: j

    allocate(stencil_map(size(dofmap,1),1,9))
    do j=1,9
      stencil_map(:,1,j)=dofmap(:,j)
    end do

  end subroutine get_m3x3_stencil_dofmap_point

!---------------------------------------------------------------------

  subroutine get_m3x3_stencil_dofmap_cross(stencil_map, dofmap)
    ! Generate a 5-point cross stencil dofmap for any given input 3x3 dofmap
    implicit none
    integer(i_def), allocatable, intent(out) :: stencil_map(:,:,:)
    integer(i_def),              intent(in)  :: dofmap(:,:)

    integer(i_def), parameter :: stencil(5,9)= reshape( [ 1,3,7,2,4, &
                                                          2,1,8,3,5, &
                                                          3,2,9,1,6, &
                                                          4,6,1,5,7, &
                                                          5,4,2,6,8, &
                                                          6,5,3,4,9, &
                                                          7,9,4,8,1, &
                                                          8,7,5,9,2, &
                                                          9,8,6,7,3], [5,9] )
    integer(i_def) :: i,j

    allocate(stencil_map(size(dofmap,1),5,9))
    do j=1,9
      do i=1,5
        stencil_map(:,i,j)=dofmap(:,stencil(i,j))
      end do
    end do

  end subroutine get_m3x3_stencil_dofmap_cross

!---------------------------------------------------------------------

  subroutine get_m3x3_stencil_dofmap_region(stencil_map, dofmap)
    ! Generate a 9-point region stencil dofmap for any given input 3x3 dofmap
    implicit none
    integer(i_def), allocatable, intent(out) :: stencil_map(:,:,:)
    integer(i_def),              intent(in)  :: dofmap(:,:)

    integer(i_def), parameter :: stencil(9,9)= reshape( [ 1,3,9,7,8,2,5,4,6, &
                                                          2,1,7,8,9,3,6,5,4, &
                                                          3,2,8,9,7,1,4,6,5, &
                                                          4,6,3,1,2,5,8,7,9, &
                                                          5,4,1,2,3,6,9,8,7, &
                                                          6,5,2,3,1,4,7,9,8, &
                                                          7,9,6,4,5,8,2,1,3, &
                                                          8,7,4,5,6,9,3,2,1, &
                                                          9,8,5,6,4,7,1,3,2], [9,9] )
    integer(i_def) :: i,j

    allocate(stencil_map(size(dofmap,1),9,9))
    do j=1,9
      do i=1,9
        stencil_map(:,i,j)=dofmap(:,stencil(i,j))
      end do
    end do

  end subroutine get_m3x3_stencil_dofmap_region

!---------------------------------------------------------------------

  subroutine get_m3x3_stencil_dofmap_1dx(stencil_map, dofmap)
    ! Generate a 5-point line (in the x-direction) stencil dofmap for
    ! any given input 3x3 dofmap
    implicit none
    integer(i_def), allocatable, intent(out) :: stencil_map(:,:,:)
    integer(i_def),              intent(in)  :: dofmap(:,:)

    integer(i_def), parameter :: stencil(5,9)= reshape( [ 1,3,2,2,3, &
                                                          2,1,3,3,1, &
                                                          3,2,1,1,2, &
                                                          4,6,5,5,6, &
                                                          5,4,6,6,4, &
                                                          6,5,4,4,6, &
                                                          7,9,8,8,9, &
                                                          8,7,9,9,7, &
                                                          9,8,7,7,8], [5,9] )
    integer(i_def) :: i,j

    allocate(stencil_map(size(dofmap,1),5,9))
    do j=1,9
      do i=1,5
        stencil_map(:,i,j)=dofmap(:,stencil(i,j))
      end do
    end do

  end subroutine get_m3x3_stencil_dofmap_1dx

!---------------------------------------------------------------------

  subroutine get_m3x3_stencil_dofmap_1dy(stencil_map, dofmap)
    ! Generate a 5-point line (in the y-direction) stencil dofmap for
    ! any given input 3x3 dofmap
    implicit none
    integer(i_def), allocatable, intent(out) :: stencil_map(:,:,:)
    integer(i_def),              intent(in)  :: dofmap(:,:)

    integer(i_def), parameter :: stencil(5,9)= reshape( [ 1,7,4,4,7, &
                                                          2,8,5,5,8, &
                                                          3,9,6,6,9, &
                                                          4,1,7,7,1, &
                                                          5,2,8,8,2, &
                                                          6,3,9,9,3, &
                                                          7,4,1,1,4, &
                                                          8,4,2,2,5, &
                                                          9,6,3,3,6], [5,9] )
    integer(i_def) :: i,j

    allocate(stencil_map(size(dofmap,1),5,9))
    do j=1,9
      do i=1,5
        stencil_map(:,i,j)=dofmap(:,stencil(i,j))
      end do
    end do

  end subroutine get_m3x3_stencil_dofmap_1dy

!---------------------------------------------------------------------

end module get_unit_test_m3x3_dofmap_mod
