submodule (epsilon_mod) epsilon_two_submod
implicit none
contains
  module function new_cart_type(x,y)
    implicit none
    integer, intent(in) :: x
    integer, intent(in) :: y
    type(cart_type) :: new_cart_type
    new_cart_type%x = x
    new_cart_type%y = y
  end function new_cart_type
end submodule epsilon_two_submod
