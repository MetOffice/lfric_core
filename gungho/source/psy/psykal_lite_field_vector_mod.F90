  module psykal_lite_field_vector_mod
    USE constants_mod, ONLY: r_def
    USE operator_mod, ONLY: operator_type, operator_proxy_type, columnwise_operator_type, columnwise_operator_proxy_type
    USE field_mod, ONLY: field_type, field_proxy_type
    IMPLICIT NONE
    CONTAINS
      subroutine invoke_0(field_norm, self_vector)
      USE scalar_mod, ONLY: scalar_type
      USE omp_lib, ONLY: omp_get_thread_num
      USE omp_lib, ONLY: omp_get_max_threads
      use mesh_mod, only: mesh_type
      implicit none
      REAL(KIND=r_def), intent(out) :: field_norm
      TYPE(field_type), intent(in) :: self_vector
      TYPE(scalar_type) global_sum
      INTEGER df
      REAL(KIND=r_def), allocatable, dimension(:,:) :: l_field_norm
      INTEGER th_idx
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      INTEGER nthreads
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      field_norm = 0.0_r_def
      ALLOCATE (l_field_norm(8,nthreads))
      l_field_norm = 0.0_r_def
      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        l_field_norm(1,th_idx) = l_field_norm(1,th_idx)+self_vector_proxy%data(df)*self_vector_proxy%data(df)
      END DO 
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        field_norm = field_norm+l_field_norm(1,th_idx)
      END DO 
      DEALLOCATE (l_field_norm)
      global_sum%value = field_norm
      field_norm = global_sum%get_sum()
      !
    END SUBROUTINE invoke_0
    subroutine invoke_1(inner_prod_field, self_vector, x_vector)
      USE scalar_mod, ONLY: scalar_type
      USE omp_lib, ONLY: omp_get_thread_num
      USE omp_lib, ONLY: omp_get_max_threads
      use mesh_mod, only: mesh_type
      implicit none      
      
      REAL(KIND=r_def), intent(out) :: inner_prod_field
      TYPE(field_type), intent(in) :: self_vector, x_vector
      TYPE(scalar_type) global_sum
      INTEGER df
      REAL(KIND=r_def), allocatable, dimension(:,:) :: l_inner_prod_field
      INTEGER th_idx
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      INTEGER nthreads
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy, x_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      x_vector_proxy = x_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      inner_prod_field = 0.0_r_def
      ALLOCATE (l_inner_prod_field(8,nthreads))
      l_inner_prod_field = 0.0_r_def

      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        l_inner_prod_field(1,th_idx) = l_inner_prod_field(1,th_idx)+self_vector_proxy%data(df)*x_vector_proxy%data(df)
      END DO 
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        inner_prod_field = inner_prod_field+l_inner_prod_field(1,th_idx)
      END DO 
      DEALLOCATE (l_inner_prod_field)
      global_sum%value = inner_prod_field
      inner_prod_field = global_sum%get_sum()
      !
    END SUBROUTINE invoke_1
    subroutine invoke_2(self_vector, alpha, y_vector)
      use mesh_mod, only: mesh_type
      implicit none      
      
      REAL(KIND=r_def), intent(in) :: alpha
      TYPE(field_type), intent(inout) :: self_vector
      TYPE(field_type), intent(in) :: y_vector
      INTEGER df
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy, y_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      y_vector_proxy = y_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        self_vector_proxy%data(df) = self_vector_proxy%data(df) + alpha*y_vector_proxy%data(df)
      END DO 
      !$omp end do
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      !$omp master
      CALL self_vector_proxy%set_dirty()
      !$omp end master
      !
      !$omp end parallel
      !
    END SUBROUTINE invoke_2
    subroutine invoke_3(self_vector, scalar)
      use mesh_mod, only: mesh_type
      implicit none      
      
      REAL(KIND=r_def), intent(in) :: scalar
      TYPE(field_type), intent(inout) :: self_vector
      INTEGER df
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        self_vector_proxy%data(df) = scalar
      END DO 
      !$omp end do
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      !$omp master
      CALL self_vector_proxy%set_dirty()
      !$omp end master
      !
      !$omp end parallel
      !
    END SUBROUTINE invoke_3
    subroutine invoke_4(alpha, self_vector, x_vector)
      use mesh_mod, only: mesh_type
      implicit none      
      
      REAL(KIND=r_def), intent(in) :: alpha
      TYPE(field_type), intent(inout) :: self_vector
      TYPE(field_type), intent(in) :: x_vector
      INTEGER df
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy, x_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      x_vector_proxy = x_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        self_vector_proxy%data(df) = alpha*self_vector_proxy%data(df) + x_vector_proxy%data(df)
      END DO 
      !$omp end do
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      !$omp master
      CALL self_vector_proxy%set_dirty()
      !$omp end master
      !
      !$omp end parallel
      !
    END SUBROUTINE invoke_4
    subroutine invoke_5(scale, self_vector)
      use mesh_mod, only: mesh_type
      implicit none      
      
      REAL(KIND=r_def), intent(in) :: scale
      TYPE(field_type), intent(inout) :: self_vector
      INTEGER df
      INTEGER ndf_any_space_1_self_vector, undf_any_space_1_self_vector
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER nlayers
      TYPE(field_proxy_type) self_vector_proxy
      !
      ! Initialise field and/or operator proxies
      !
      self_vector_proxy = self_vector%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = self_vector_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => self_vector%get_mesh()
      !
      ! Initialise number of DoFs for any_space_1_self_vector
      !
      ndf_any_space_1_self_vector = self_vector_proxy%vspace%get_ndf()
      undf_any_space_1_self_vector = self_vector_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=1,self_vector_proxy%vspace%get_last_dof_owned()
        self_vector_proxy%data(df) = scale*self_vector_proxy%data(df)
      END DO 
      !$omp end do
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      !$omp master
      CALL self_vector_proxy%set_dirty()
      !$omp end master
      !
      !$omp end parallel
      !
    END SUBROUTINE invoke_5
  end module psykal_lite_field_vector_mod
