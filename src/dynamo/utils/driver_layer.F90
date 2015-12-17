!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> A container module for various routines that are part of the driver layer
!> and hence sit outside of the PsyKAl structure
module driver_layer

  implicit none

  private 
  public :: interpolated_output

contains
!> @brief An procedure to interpolate a field onto a defined set of points
!> @details An algorithm for computing interpolated output on a defined set of points
!>          and dumping to file.
!>          Effectively this routine evaluates a field using the known finite
!>          polynomial expansion of the field at a given point and
!>          dumps the value to file
!> @deprecated This is a tempoary implementation until a proper i/o + plotting
!> stategy is implemented
!> @param[in] n_out The number of output fields to generate from f
!> @param[in] f     A field to compute output data from
!> @param[in] mesh  The mesh object the model runs on
!> @param[in] chi   A 3D coordinate field
!> @param[in] fname The name of the field to be output
  subroutine interpolated_output(n_out, f, mesh, chi, fname) 

    use log_mod,                   only: log_event, log_scratch_space, LOG_LEVEL_INFO
    use constants_mod,             only: r_def, str_max_filename
    use field_mod,                 only: field_type
    use find_output_cell_mod,      only: find_output_cell
    use evaluate_output_field_mod, only: evaluate_output_field    
    use mesh_mod,                  only: mesh_type, domain_limits
    use coord_transform_mod,       only: xyz2llr

    implicit none
! Mesh
    type (mesh_type),   intent(in) :: mesh
! Dimension of input field
    integer,            intent(in) :: n_out
! Field to output
    type( field_type ), intent(in) :: f(n_out)
! Coodinate fields
    type( field_type ), intent(in) :: chi(3)  
! name of field
    character(str_max_filename), intent(in) :: fname

    integer                       :: nx(3), i, j, k, out_cell, dir 
    integer, parameter            :: OUTPUT_UNIT = 21
    real(kind=r_def), allocatable :: x_out(:,:,:,:), f_out(:,:,:,:)  
    real(kind=r_def)              :: dx(3)

    type (domain_limits) :: domain_size


! Create uniform grid for output (nx,ny,nz)
    nx(1) = 100
    nx(2) = 50
    nx(3) = 11 

    allocate( x_out(3,nx(3),nx(2),nx(1)), f_out(n_out,nx(3),nx(2),nx(1)) )

    domain_size = mesh%get_domain_size()

! Create regular domain
    dx(1) = (domain_size%maximum%x - domain_size%minimum%x)/real(nx(1)-1)
    dx(2) = (domain_size%maximum%y - domain_size%minimum%y)/real(nx(2)-1)
    dx(3) = (domain_size%maximum%z - domain_size%minimum%z)/real(nx(3)-1)

    do i = 1,nx(1)
      do j = 1,nx(2)
        do k = 1,nx(3)
          x_out(1,k,j,i) = real(i-1)*dx(1) + domain_size%minimum%x
          x_out(2,k,j,i) = real(j-1)*dx(2) + domain_size%minimum%y
          x_out(3,k,j,i) = real(k-1)*dx(3) + domain_size%minimum%z
        end do
      end do
    end do

! For each point on regular grid find the computational cell that it is located
! in and evaluate the field at that point -> This is very inefficient for large
! grids
    do i = 1,nx(1)
      do j = 1,nx(2)
! find grid cell each output point lives in
        out_cell = find_output_cell( chi, x_out(:,1,j,i) )
! evaluate field at output points
        do dir = 1,n_out
          call evaluate_output_field( mesh, f(dir), chi, x_out(:,:,j,i), out_cell, nx(3), f_out(dir,:,j,i) )
        end do
      end do
    end do


    write( log_scratch_space, '(A,A)' ) 'Writing interpolated output: ',trim(fname)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  
    open(OUTPUT_UNIT, file = trim(fname), status = "replace")
    write(OUTPUT_UNIT,*) 'II = ',nx(1),';'
    write(OUTPUT_UNIT,*) 'JJ = ',nx(2),';'
    write(OUTPUT_UNIT,*) 'KK = ',nx(3),';'
    write(OUTPUT_UNIT,*) 'LL = ',3+n_out,';'
    write(OUTPUT_UNIT,*) 'data = ['
    do i = 1,nx(1)
      do j = 1,nx(2)
        do k = 1,nx(3)    
          if ( n_out == 1) then
            write(OUTPUT_UNIT,'(4e18.8e3)') x_out(1,k,j,i), x_out(2,k,j,i), x_out(3,k,j,i), &
                                            f_out(1,k,j,i)
          else
            write(OUTPUT_UNIT,'(6e18.8e3)') x_out(1,k,j,i), x_out(2,k,j,i), x_out(3,k,j,i), &
                                            f_out(1,k,j,i), f_out(2,k,j,i), f_out(3,k,j,i)
          end if 
        end do
      end do
    end do
    write(OUTPUT_UNIT,*) '];'
    write(OUTPUT_UNIT,*) ' '
    write(OUTPUT_UNIT,*) 'x=zeros(LL,KK,JJ,II); id=1; '
    write(OUTPUT_UNIT,*) 'for i=1:II; for j=1:JJ; for k=1:KK;'
    write(OUTPUT_UNIT,*) 'for l=1:LL;'
    write(OUTPUT_UNIT,*) 'x(l,k,j,i) = data(id,l);'
    write(OUTPUT_UNIT,*) 'end;' 
    write(OUTPUT_UNIT,*) 'id = id + 1; end; end; end;'
    close(OUTPUT_UNIT)
  
  end subroutine interpolated_output 

end module driver_layer
