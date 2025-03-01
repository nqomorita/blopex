module blopex_fortran_hold_vars
  use iso_c_binding
  integer(c_int) :: N_hold, NDOF_hold, M_hold, prec_type_hold
  integer(c_int), pointer :: index_hold(:), item_hold(:)
  real(c_double), pointer :: A_hold(:)
  real(c_double), pointer :: B_hold(:)
  real(c_double), allocatable :: Diag(:)
end module blopex_fortran_hold_vars

module blopex_fortran_driver
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none

contains

  subroutine blopex_lobpcg_solve( &
      & N, NZ, index, item, A, n_eigs, maxit, tol, loglevel, &
      & eigen_val, eigen_vec, prec_type, B)
    implicit none
    integer(c_int) :: N, NZ, i, j
    integer(c_int) :: is_B, is_T, prec_type
    integer(c_int), target :: index(0:N), item(NZ)
    integer(c_int) :: n_eigs, maxit, mat_n, loglevel
    real(c_double) :: tol
    real(c_double), target :: A(NZ)
    real(c_double) :: eigen_val(n_eigs)
    real(c_double) :: eigen_vec(N,n_eigs)
    real(c_double) :: eigen_vec_temp(N*n_eigs)
    real(c_double), optional :: B(N)
    external :: blopex_fortran_opA
    external :: blopex_fortran_opB
    external :: blopex_fortran_opT

    if(N < n_eigs) n_eigs = N

    !call set_matvec_A(N, n_eigs, index, item, A)
    N_hold = N
    M_hold = n_eigs
    index_hold => index
    item_hold  => item
    A_hold => A

    is_B = 0
    if(present(B))then
      is_B = 1
      call set_matvec_B(N, B)
    endif

    is_T = 0
    if(prec_type > 0)then
      is_T = 1
      prec_type_hold = prec_type
      call set_preconditioning(N, NZ, index, item, A, Diag)
    endif

    call blopex_lobpcg_solve_c(n_eigs, maxit, tol, N, &
      & blopex_fortran_opA, blopex_fortran_opB, blopex_fortran_opT, &
      & is_B, is_T, &
      & loglevel, eigen_val, eigen_vec_temp)

    do i = 1, n_eigs
      do j = 1, N
        eigen_vec(j,i) = eigen_vec_temp(n_eigs*(i-1) + j)
      enddo
    enddo
  end subroutine blopex_lobpcg_solve

  subroutine set_matvec_B(N, B)
    implicit none
    integer(c_int) :: N
    real(c_double), target :: B(N)
    B_hold => B
  end subroutine set_matvec_B

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

subroutine blopex_fortran_opB(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, shift
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  do j = 1, M_hold
    shift = N_hold*(j-1)
    do i = 1, N_hold
      b(i+shift) = B_hold(i)*a(i+shift)
    enddo
  enddo
end subroutine blopex_fortran_opB

subroutine blopex_fortran_opT(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE, in, jn, k, shift
  real(c_double) :: s
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  if(prec_type_hold == 1)then
    !> diag
    do i = 1, M_hold
      shift = N_hold*(i-1)
      do j = 1, N_hold
        b(j+shift) = a(j+shift)*Diag(j)
      enddo
    enddo

  elseif(prec_type_hold == 2)then
    !> sor
    do k = 1, M_hold
      shift = N_hold*(k-1)
      do i = 1, N_hold
        s = a(i+shift)
        jS = index_hold(i-1) + 1
        jE = index_hold(i  )
        do j = jS, jE
          jn = item_hold(j)
          if(jn < i)then
            s = s - A_hold(j)*a(jn+shift)
          endif
        enddo
        a(i+shift) = s*Diag(i)
      enddo

      do i = 1, N_hold
        a(i+shift) = a(i+shift)/Diag(i)
      enddo

      do i = N_hold, 1, -1
        s = a(i+shift)
        jS = index_hold(i-1) + 1
        jE = index_hold(i  )
        do j = jE, jS, -1
          jn = item_hold(j)
          if(i < jn)then
            s = s - A_hold(j)*a(jn+shift)
          endif
        enddo
        a(i+shift) = s*Diag(i)
        b(i+shift) = a(i+shift)
      enddo
    enddo
  endif
end subroutine blopex_fortran_opT
