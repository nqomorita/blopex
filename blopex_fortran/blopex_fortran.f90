module blopex_fortran_hold_vars
  !integer(4) :: N, NDOF, M
  !real(8) :: A
end module blopex_fortran_hold_vars

module blopex_fortran_driver
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none

contains

  subroutine blopex_lobpcg_solve()
    integer(c_int) :: n_eigs, maxit, mat_n
    real(c_double) :: tol
    external :: blopex_fortran_opA

    n_eigs = 1
    maxit = 100000
    mat_n = 10
    tol = 1.0d-6

    call blopex_lobpcg_solve_c(n_eigs, maxit, tol, mat_n, blopex_fortran_opA)
  end subroutine blopex_lobpcg_solve

end module blopex_fortran_driver

subroutine blopex_fortran_opA(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum
  real(c_double) :: a(*), b(*)
end subroutine blopex_fortran_opA
