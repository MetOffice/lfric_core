!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!>@brief Abstract vector type for fields to use the new solver API and 
!! extended vector types for particular solvers.

module field_vector_mod
  use constants_mod,                 only : i_def, r_def
  use vector_mod,                    only : abstract_vector_type
  use field_mod,                     only : field_type
  use function_space_collection_mod, only : function_space_collection
  use finite_element_config_mod,     only : element_order
  use log_mod,                       only : log_event, LOG_LEVEL_ERROR, &
                                            log_scratch_space
  use psykal_lite_mod,               only : invoke_set_field_scalar, &
                                            invoke_inc_axpy,         &
                                            invoke_inc_xpby,         &
                                            invoke_inner_prod,       &
                                            invoke_X_innerproduct_X, &
                                            invoke_scale_field_data

  implicit none
  private 

  type, public, extends(abstract_vector_type) :: field_vector_type
     private
     !> The array holding the fields
     type(field_type), allocatable :: vector(:)
     !> How many fields
     integer(kind=i_def)           :: nfields
     !> Whether a field has been set in this position
     logical,allocatable           :: field_set(:)     
   contains
!   public procedures of the type
     !> Import a field into the vector at the given position in the array
     !> @param[in] field  field_type, the field to copy
     !> @param[in] pos  the position in the array
     procedure, public  :: import_field
     !> Export a field from the vector at the given position in the array
     !> @param[out] field  field_type, the field to copy to
     !> @param[in] pos  the position in the array
     procedure          :: export_field

!   public procedures of the API overidden in the type
     !> set the vector to a scalar value
     !> @param[in] scalar  real the scalar value
     procedure, public  :: set_scalar => set_field_vector_scalar
     !> Compute y = y + alpha * x
     !> @param[in] alpha  real
     !> @param[inout] x  vector, an array of fields
     procedure, public  :: axpy => axpy_field_vector
     !> Compute norm of the field vector, returns a real scalar
     !! n = sqrt( sum_i( f_i*f_i )) where f_i is each field
     procedure, public  :: norm => norm_field_vector
     !> Compute the dot product of two field_vectors returns a real scalar
     !> @param[in] x field_vector to dot self with
     procedure, public  :: dot => dot_field_vector
     !> multiply a field vector by a scalar
     !> @param [in] scalar real
     procedure, public  :: scale => scale_field_vector
     !> Compute y = alpha*y + x
     !> @param[in] alpha  real
     !> @param[inout] x  vector, an array of fields
     procedure, public  :: aypx => aypx_field_vector


!   private procedure of the type which either implement the vector API
!   or allow encapsulation not to be broken
     procedure, private :: set_field_vector_scalar
     procedure, private :: axpy_field_vector
     procedure, private :: axpy_field
     procedure, private :: norm_field_vector
     procedure, private :: dot_field_vector
     procedure, private :: dot_field
     procedure, private :: aypx_field_vector
     procedure, private :: aypx_field
     procedure, private :: scale_field_vector
!   infractructure procedures 
     procedure :: field_vector_type_assign
     generic            :: assignment(=) => field_vector_type_assign
     final              :: field_vector_destroy
  end type field_vector_type
  
  interface field_vector_type
     module procedure field_vector_constructor
  end interface

contains

  ! compute the norm of a field vector  
  function norm_field_vector(self) result(normal)
    class(field_vector_type), intent(in) :: self
    real(kind=r_def)    :: normal
    integer(kind=i_def) :: fctr
    real(kind=r_def)    :: field_norm
    normal = 0.0_r_def
        
    do fctr = 1, self%nfields
       ! check we have a field set in each position
       if(.not.self%field_set(fctr)) then 
          write(log_scratch_space,'(A,I0,A)') & 
               "field_vector_mod:norm_field_vector: field at position", &
               fctr," has not been set"
          call log_event(log_scratch_space,LOG_LEVEL_ERROR)
       end if
       call invoke_X_innerproduct_X(field_norm, self%vector(fctr))
       normal=normal + field_norm
    end do
    normal = sqrt(normal)
  end function norm_field_vector

  ! compute the dot of inner product of a field vector
  function dot_field_vector(self, x) result(dot_prod)
    class(field_vector_type),    intent(in) :: self
    class(abstract_vector_type), intent(in) :: x
    real(kind=r_def)                        :: dot_prod
    real(kind=r_def)                        :: inner_prod_field
    integer(kind=i_def)                     :: fctr

    select type(x)
    type is(field_vector_type)
       dot_prod = 0.0_r_def
       do fctr = 1, self%nfields
         ! check we have a field set in each position
          if(.not.self%field_set(fctr)) then 
             write(log_scratch_space,'(A,I0,A)') & 
                  "field_vector_mod:dot_field_vector: field at position", &
                  fctr," has not been set"
             call log_event(log_scratch_space,LOG_LEVEL_ERROR)
          end if
          inner_prod_field = x%dot_field(self%vector(fctr), fctr)
          dot_prod = dot_prod + inner_prod_field
       end do
    class default
       write(log_scratch_space,'(A)') &
            "field_vector_mod:dot_field_vector: type of x is not field_vector_type"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end select
  end function dot_field_vector

  ! private function which computes the dot product of a field with a field at position
  ! pos in a field vector array
  function dot_field(self, y, pos) result(ip)
    class(field_vector_type), intent(in) :: self
    class(field_type),        intent(in) :: y
    integer(kind=i_def),      intent(in) :: pos
    real(kind=r_def)                     :: ip

    ! check we have a field set in this position
    if(.not.self%field_set(pos)) then 
       write(log_scratch_space,'(A,I0,A)') & 
            "field_vector_mod:axpy_field: field at position", &
            pos," has not been set"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    call invoke_inner_prod(y,self%vector(pos),ip)
    
  end function dot_field

  ! computes y = alpha * x + y on a field vector
  subroutine axpy_field_vector(self, alpha, x)
    class(field_vector_type),    intent(inout) :: self
    real(kind=r_def),            intent(in)    :: alpha
    class(abstract_vector_type), intent(inout) :: x
    integer(kind=i_def) :: fctr

    select type(x)
    type is(field_vector_type)
    
       do fctr = 1, self%nfields
          ! check we have a field set in each position
          if(.not.self%field_set(fctr)) then 
             write(log_scratch_space,'(A,I0,A)') & 
                  "field_vector_mod:axpy_field_vector: field at position", &
                  fctr," has not been set"
             call log_event(log_scratch_space,LOG_LEVEL_ERROR)
          end if
          ! pass the field to
          call x%axpy_field(alpha, self%vector(fctr), fctr)
       end do
    class default
       write(log_scratch_space,'(A)') &
            "field_vector_mod:axpy_field_vector: type of x is not field_vector_type"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end select
  end subroutine axpy_field_vector

  ! computes y = alpha *x(pos) + y
  ! where y is a field and x(self) is a field vector
  subroutine axpy_field(self, alpha, y, pos)
    class(field_vector_type), intent(inout) :: self
    real(kind=r_def),         intent(in)    :: alpha
    class(field_type),        intent(inout) :: y
    integer(kind=i_def),      intent(in)    :: pos

    ! check we have a field set in this position
       if(.not.self%field_set(pos)) then 
          write(log_scratch_space,'(A,I0,A)') & 
               "field_vector_mod:axpy_field: field at position", &
               pos," has not been set"
          call log_event(log_scratch_space,LOG_LEVEL_ERROR)
       end if
       call invoke_inc_xpby(y, alpha, self%vector(pos))

  end subroutine axpy_field

  ! sets a field to a scalar. If a field has not previously been set in
  ! the vector (array of fields) then there is no function space information
  ! so can't set the field to the scalar value. The procedure throws an error.
  subroutine set_field_vector_scalar(self, scalar)
    class(field_vector_type), intent(inout) :: self
    real(kind=r_def),         intent(in)    :: scalar
    integer(kind=i_def) :: fctr

    do fctr = 1, self%nfields
       ! check we have a field set in each position
       if(.not.self%field_set(fctr)) then 
          write(log_scratch_space,'(A,I0,A)') & 
               "field_vector_mod:set_field_vector_scalar: field at position", &
               fctr," has not been set"
          call log_event(log_scratch_space,LOG_LEVEL_ERROR)
       end if
       call invoke_set_field_scalar(scalar, self%vector(fctr))
! A placeholder for PSyclone built-ins support (to be agreed on implementation)
!        call invoke( set_field_scalar(scalar, self%vector(fctr)) )
       self%field_set(fctr) = .true.
    end do
    
  end subroutine set_field_vector_scalar

  ! computes y = alpha * y + x
  ! where self is y
  subroutine aypx_field_vector(self, alpha, x)
    class(field_vector_type),    intent(inout) :: self
    real(kind=r_def),            intent(in)    :: alpha
    class(abstract_vector_type), intent(inout) :: x
    integer(kind=i_def) :: fctr

    select type(x)
    type is(field_vector_type)
    
       do fctr = 1, self%nfields
          ! check we have a field set in each position
          if(.not.self%field_set(fctr)) then 
             write(log_scratch_space,'(A,I0,A)') & 
                  "field_vector_mod:aypx_field_vector: field at position", &
                  fctr," has not been set"
             call log_event(log_scratch_space,LOG_LEVEL_ERROR)
          end if
          ! pass the field to
          call x%aypx_field(alpha, self%vector(fctr), fctr)
       end do
    class default
       write(log_scratch_space,'(A)') &
            "field_vector_mod:aypx_field_vector: type of x is not field_vector_type"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end select
  end subroutine aypx_field_vector

  ! private procedure to compute y = alpha * y + x(pos)
  ! where y is a field and x is a vector of fields
  subroutine aypx_field(self, alpha, y, pos)
    class(field_vector_type), intent(inout) :: self
    real(kind=r_def),         intent(in)    :: alpha
    class(field_type),        intent(inout) :: y
    integer(kind=i_def),      intent(in)    :: pos

    ! check we have a field set in this position
    if(.not.self%field_set(pos)) then 
       write(log_scratch_space,'(A,I0,A)') & 
            "field_vector_mod:aypx_field: field at position", &
            pos," has not been set"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    call invoke_inc_axpy(alpha, y, self%vector(pos))
  end subroutine aypx_field

  ! multiply the field vector by a scalar
  subroutine scale_field_vector(self, scalar)
    class(field_vector_type), intent(inout) :: self
    real(kind=r_def),         intent(in)    :: scalar
    integer(kind=i_def) :: fctr

    do fctr = 1, self%nfields
       ! check we have a field set in each position
       if(.not.self%field_set(fctr)) then 
          write(log_scratch_space,'(A,I0,A)') & 
               "field_vector_mod:set_field_vector_scalar: field at position", &
               fctr," has not been set"
          call log_event(log_scratch_space,LOG_LEVEL_ERROR)
       end if
       call invoke_scale_field_data(scalar, self%vector(fctr))
    end do
    
  end subroutine scale_field_vector

  ! Copy a field into the field vector
  subroutine import_field(self, field, position)
    class(field_vector_type ), intent(inout) :: self
    type(field_type),             intent(in)    :: field
    integer(kind=i_def),          intent(in)    :: position

    if(position > self%nfields) then 
       write(log_scratch_space,'(A,2(":",I2))') & 
            "field_vector_mod:field position bigger than nfields", &
            self%nfields,position
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    
    self%vector(position) = field
    self%field_set(position) = .true.
  end subroutine import_field

   ! Copies a field from the nominated position. The procedure a field has been
   ! set in this postition
   ! @param[in] field type(field_type), the field to export
   ! @param[in] position integer postion in the field array where the field is 
   ! located.
  subroutine export_field(self, field, position) 
    class(field_vector_type),  intent(in)    :: self
    type(field_type),             intent(inout) :: field
    integer(kind=i_def),          intent(in)    :: position

    if(position > self%nfields) then 
       write(log_scratch_space,'(A,2(":",I2))') & 
            "field_vector_mod:export_field:field position bigger than nfields", &
            self%nfields,position
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if

    if(.not.self%field_set(position)) then 
       write(log_scratch_space,'(A,I0,A)') & 
            "field_vector_mod:export_field:field at position",position, &
            " has not been set"
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    
    field=self%vector(position)
  end subroutine export_field

  ! allocates the array of fields to the given size
  function field_vector_constructor(nfields) result(self)
    integer(kind=i_def), intent(in) :: nfields
    type(field_vector_type) :: self

    self%nfields = nfields
    allocate(self%vector(self%nfields))
    allocate(self%field_set(self%nfields))
    self%field_set = .false.
  end function field_vector_constructor

  ! The destructor/finalizer
  subroutine field_vector_destroy(self)
    type(field_vector_type), intent(inout) :: self
    if(allocated(self%vector)) then
       deallocate(self%vector)
    end if
    if(allocated(self%field_set)) then
       deallocate(self%field_set)
    end if
  end subroutine field_vector_destroy

  subroutine field_vector_type_assign(dest, source)
    class(field_vector_type), intent(out) :: dest
    class(field_vector_type), intent(in)  :: source
    integer :: pos
    
    ! make field_vector
    
    dest%nfields = source%nfields
    allocate(dest%vector(dest%nfields))
    allocate(dest%field_set(dest%nfields))
    ! if fields don't exist set false, otherwise copy and set true
    do pos = 1, source%nfields
       if(.not.source%field_set(pos)) then 
          dest%field_set(pos) = .false.
       else
          dest%vector(pos) = source%vector(pos)
          dest%field_set(pos) = .true.
       end if
    end do
    
  end subroutine field_vector_type_assign

end module field_vector_mod
