!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module psykal_builtin_light_mod

  use constants_mod, only : i_def, i_long, r_def
  use field_mod,     only : field_type, field_proxy_type

  implicit none

  public

contains

  !----------------------------------------------------------------------------
  subroutine invoke_convert_cart2sphere_vector( field, coords)
    use coord_transform_mod, only: cart2sphere_vector
    implicit none
    type(field_type), intent(inout) :: field(3)
    type(field_type), intent(in)    :: coords(3)

    type(field_proxy_type) :: f_p(3), x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: vector_in(3), vector_out(3), xyz(3)

    do i = 1,3
      f_p(i) = field(i)%get_proxy()
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = f_p(1)%vspace%get_last_dof_annexed()

    do df = 1, undf
      vector_in(:)  = (/ f_p(1)%data(df), f_p(2)%data(df), f_p(3)%data(df) /)
      xyz(:)        = (/ x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df) /)
      vector_out(:) = cart2sphere_vector(xyz, vector_in)
      f_p(1)%data(df) = vector_out(1)
      f_p(2)%data(df) = vector_out(2)
      f_p(3)%data(df) = vector_out(3)
    end do

    call f_p(1)%set_dirty()
    call f_p(2)%set_dirty()
    call f_p(3)%set_dirty()

  end subroutine invoke_convert_cart2sphere_vector

  !----------------------------------------------------------------------------
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

end module psykal_builtin_light_mod
