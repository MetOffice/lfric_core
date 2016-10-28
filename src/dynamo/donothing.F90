program donothing
! dependency analyser doesn't look at the tests. These kernels only appear
! in tests so a fake "driver" is needed to spoof the dependency analyser
! so the kernels get compiled.
  use columnwise_op_app_kernel_mod, only : columnwise_op_app_kernel_type
  use columnwise_op_asm_kernel_mod, only : columnwise_op_asm_kernel_type
  use columnwise_op_appinv_kernel_mod, only : columnwise_op_appinv_kernel_type
  use columnwise_op_mul_kernel_mod, only : columnwise_op_mul_kernel_type
  use columnwise_op_scaledadd_kernel_mod, only : columnwise_op_scaledadd_kernel_type

  implicit none

  write(*,*) "doing nothing..."

end program donothing
