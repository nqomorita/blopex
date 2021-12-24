module blopex_fortran_hold_vars
  use iso_c_binding
  integer(c_int) :: N_hold, NDOF_hold, M_hold
  integer(c_int), pointer :: index_hold(:), item_hold(:)
  real(c_double), pointer :: A_hold(:)
  real(c_double), allocatable :: Diag(:)
end module blopex_fortran_hold_vars

module blopex_fortran_driver
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none

contains

  subroutine blopex_lobpcg_solve(N, NZ, index, item, A, n_eigs, maxit, tol)
    implicit none
    integer(c_int) :: N, NZ
    integer(c_int), target :: index(0:N), item(NZ)
    integer(c_int) :: n_eigs, maxit, mat_n
    real(c_double) :: tol
    real(c_double), target :: A(NZ)
    external :: blopex_fortran_opA
    external :: blopex_fortran_opT

    if(N < n_eigs) n_eigs = N

    N_hold = N
    M_hold = n_eigs
    index_hold => index
    item_hold  => item
    A_hold => A

    call set_preconditioning(N, NZ, index, item, A, Diag)

    call blopex_lobpcg_solve_c(n_eigs, maxit, tol, N, blopex_fortran_opA, blopex_fortran_opT)
  end subroutine blopex_lobpcg_solve

  subroutine set_preconditioning(N, NZ, index, item, A, Diag)
    implicit none
    integer(c_int) :: N, NZ
    integer(c_int) :: index(0:N), item(NZ)
    integer(c_int) :: i, jS, jE, j, in
    real(c_double) :: A(NZ)
    real(c_double), allocatable :: Diag(:)

    allocate(Diag(N), source = 0.0d0)

    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        if(i == in)then
          Diag(i) = 1.0d0/A(j)
        endif
      enddo
    enddo
  end subroutine set_preconditioning

end module blopex_fortran_driver

subroutine blopex_fortran_opA(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, jS, jE, j, in, k, shift
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  b = 0.0d0
  do k = 1, M_hold
    shift = N_hold*(k-1)
    do i = 1, N_hold
      jS = index_hold(i-1) + 1
      jE = index_hold(i)
      do j = jS, jE
        in = item_hold(j) + shift
        b(i+shift) = b(i+shift) + A_hold(j)*a(in)
      enddo
    enddo
  enddo
end subroutine blopex_fortran_opA

subroutine blopex_fortran_opT(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, shift
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)
  !b = a
  do i = 1, M_hold
    shift = N_hold*(i-1)
    do j = 1, N_hold
      b(j+shift) = a(j+shift)*Diag(j)
    enddo
  enddo
end subroutine blopex_fortran_opT
