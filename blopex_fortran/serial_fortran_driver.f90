program main
  use iso_c_binding
  use blopex_fortran_driver
  implicit none

  write(*,"(a)")"* blopex fortran driver"

  call blopex_lobpcg_solve()

contains

end program main
