!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!>@brief Contains wrapper routines for working on arrays of fields
!!       Effectively these are generally wrapper functions to pointwise kernels
module field_bundle_mod

  use constants_mod,                   only: i_def, r_def
  use field_mod,                       only: field_type
  use function_space_collection_mod,   only: function_space_collection
  use finite_element_config_mod,       only: element_order

  implicit none

contains

!> @brief Create a bundle y of fields on the same function space as bundle x
!> @param[in]  x Input bundle to clone
!> @param[out] y Output bundle to contain fields of the same type as x
!> @param[in]  bundle_size Number of fields in the bundle

  subroutine clone_bundle(x, y, bundle_size)

    use finite_element_config_mod, only: element_order

    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size)
    type(field_type), intent(inout) :: y(bundle_size)

    integer(i_def)            :: fs_handle, mesh_id
    integer                   :: i

    do i = 1,bundle_size   
      fs_handle = x(i)%which_function_space()
      mesh_id = x(i)%get_mesh_id()
      y(i) = field_type( vector_space = &
               function_space_collection%get_fs(mesh_id,element_order,fs_handle) )

    end do
  end subroutine clone_bundle
!=============================================================================!

!> @brief Set all fields in a bundle x to the scalar value a
!> @param [in] a Scalar
!> @param [inout] x Field bundle
!> @param [in] bundle_size Number of fields in the bundle x
  subroutine set_bundle_scalar(a, x, bundle_size)
    use psykal_lite_mod, only: invoke_set_field_scalar
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: x(bundle_size)
    real(kind=r_def), intent(in)    :: a
    integer ::i   
 
    do i = 1,bundle_size
      call invoke_set_field_scalar(a, x(i))
! A placeholder for PSyclone built-ins support (to be agreed on implementation)
!       call invoke( set_field_scalar(a, x(i)) )
    end do
  end subroutine set_bundle_scalar
!=============================================================================!

!> @brief Compute the inner product of two field bundles
!> @result a Inner product
!> @param [in] x First field bundle
!> @param [in] y Second field bundle
!> @param [in] bundle_size Number of fields in the bundle x
  function bundle_inner_product(x, y, bundle_size) result(a)
    use psykal_lite_mod, only: invoke_inner_prod
    implicit none
    integer,          intent(in) :: bundle_size
    type(field_type), intent(in) :: x(bundle_size), y(bundle_size)
    real(kind=r_def) :: a
    real(kind=r_def) :: b
    integer :: i

    a = 0.0_r_def
    do i = 1,bundle_size
      call invoke_inner_prod(x(i), y(i), b)
      a = a + b
    end do
  end function bundle_inner_product
!=============================================================================!

!> @brief Compute z = a*x + y for bundles x, y, z and scalar a
!> @param [in] a Scalar
!> @param [in] x First field bundle
!> @param [in] y Second field bundle
!> @param [inout] z Result field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine bundle_axpy(a, x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_axpy
    implicit none
    integer,          intent(in)    :: bundle_size
    real(kind=r_def), intent(in)    :: a
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_axpy(a, x(i), y(i), z(i))
    end do
  end subroutine bundle_axpy
!=============================================================================!

!> @brief Copy the data from bundle x to bundle y (y = x)
!> @param [in] x First field bundle
!> @param [inout] y Second field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine copy_bundle(x, y, bundle_size)
    use psykal_lite_mod, only: invoke_copy_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size)
    type(field_type), intent(inout) :: y(bundle_size)  
    integer :: i
    
    do i = 1,bundle_size
      call invoke_copy_field_data(x(i), y(i))
! A placeholder for PSyclone built-ins support (to be agreed on implementation)
!       call invoke( copy_field(x(i), y(i)) )
    end do
  end subroutine copy_bundle
!=============================================================================!
!> @brief Compute z = x - y for field bundles x, y and z
!> @param [in] x First field bundle
!> @param [in] y Second field bundle
!> @param [inout] z Result field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine minus_bundle(x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_minus_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_minus_field_data(x(i), y(i), z(i))
    end do
  end subroutine minus_bundle
!=============================================================================!
!> @brief Compute y = a*x for bundles x and y and scaler a
!> @param [in] a Scalar
!> @param [in] x First field bundle
!> @param [inout] y Second field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine bundle_ax(a, x, y, bundle_size)
    use psykal_lite_mod, only: invoke_copy_scaled_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: y(bundle_size)
    type(field_type), intent(in)    :: x(bundle_size)
    real(kind=r_def), intent(in)    :: a
    integer :: i
   
    do i = 1,bundle_size
      call invoke_copy_scaled_field_data( a, x(i), y(i) )
    end do
  end subroutine bundle_ax
!=============================================================================!
!> @brief Divide the dofs in bundle x by those in bundle y
!> @param [inout] x First field bundle
!> @param [inout] y Second field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine bundle_divide(x, y, bundle_size)
    use psykal_lite_mod, only: invoke_divide_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: x(bundle_size), y(bundle_size)
    integer :: i    

    do i = 1,bundle_size
      call invoke_divide_field_data(x(i),y(i))
    end do
  end subroutine bundle_divide
!=============================================================================!
!> @brief Write the min and max values of a field bundle to the log file
!> @param [in] x First field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine bundle_minmax(x, bundle_size)

    use psykal_lite_mod,           only: invoke_copy_field_data
    use log_mod,                   only: lOG_LEVEL_INFO   
    use finite_element_config_mod, only: element_order

    implicit none
    integer,          intent(in) :: bundle_size
    type(field_type), intent(in) :: x(bundle_size)
    type(field_type)             :: y
    integer(i_def)               :: fs_handle, mesh_id
    integer                      :: i

! This has strange syntax as psykal_lite_modclone doesnt like calls to
! type-bound procedures of field arrays
    do i = 1,bundle_size
      fs_handle = x(1)%which_function_space()
      mesh_id = x(i)%get_mesh_id()

      y = field_type( vector_space = &
               function_space_collection%get_fs(mesh_id,element_order,fs_handle) )

      call invoke_copy_field_data( x(1), y ) 
! A placeholder for PSyclone built-ins support (to be agreed on implementation)
!       call invoke( copy_field(x(1), y) )
      call y%log_minmax(LOG_LEVEL_INFO, 'field')
    end do
  end subroutine bundle_minmax

!=============================================================================!
!> @brief Compute z = a*x + b*y for bundles x, y, z and scalars a and b
!> @param [in] a Scalar
!> @param [in] b Scalar
!> @param [in] x First field bundle
!> @param [in] y Second field bundle
!> @param [inout] z Result field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine bundle_axpby(a, x, b, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_axpby
    implicit none
    integer,          intent(in)    :: bundle_size
    real(kind=r_def), intent(in)    :: a, b
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_axpby(a, x(i), b, y(i), z(i))
    end do
  end subroutine bundle_axpby
!=============================================================================!
!> @brief Compute z = x + y for field bundles x, y and z
!> @param [in] x First field bundle
!> @param [in] y Second field bundle
!> @param [inout] z Result field bundle
!> @param [in] bundle_size Number of fields in the bundle
  subroutine add_bundle(x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_plus_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_plus_field_data(x(i), y(i), z(i))
    end do
  end subroutine add_bundle
!=============================================================================!
 
end module field_bundle_mod

