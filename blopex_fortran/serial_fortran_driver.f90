program main
  use iso_c_binding
  use blopex_fortran_driver
  implicit none
  integer(4) :: N, NZ
  integer(4) :: n_eigs, maxit, loglevel
  integer(4), allocatable :: index(:), item(:)
  real(8) :: tol
  real(8), allocatable :: A(:)
  character :: fin*100

  write(*,"(a)")"* blopex fortran driver"

  call get_input_arg(fin, n_eigs, maxit, loglevel, tol)

  call input_from_matrix_market_csr(fin, N, NZ, index, item, A)

  call blopex_lobpcg_solve(N, NZ, index, item, A, n_eigs, maxit, tol, loglevel)

contains

  subroutine get_input_arg(fin, n_eigs, maxit, loglevel, tol)
    implicit none
    integer(4) :: count
    integer(4) :: n_eigs, maxit, loglevel
    real(8) :: tol
    character :: fin*100

    count = iargc()
    if(count == 1)then
      call getarg(1, fin)
    else
      stop "please enter input file name"
    endif

    open(10, file = "setting.dat", status = "old")
      read(10,*)n_eigs
      read(10,*)maxit
      read(10,*)tol
      read(10,*)loglevel
    close(10)

    write(*,"(a,i12)")     "* n_eigs  :", n_eigs
    write(*,"(a,i12)")     "* maxiter :", maxit
    write(*,"(a,1pe12.5)") "* tol     :", tol
    write(*,"(a,i12)")     "* loglevel:", loglevel
  end subroutine get_input_arg

  subroutine input_from_matrix_market_mm(fin, N, A)
    implicit none
    integer(4) :: N, NZ, i, j, in, ierr
    real(8) :: val
    real(8), allocatable :: A(:,:)
    character :: fin*100, ctemp*10
    logical :: is_first

    open(20, file = trim(fin), status = "old")
    in = 0
    do
      read(20, "(a)", iostat = ierr) ctemp
      if(0 /= ierr) exit

      if(ctemp(1:1) == "%") cycle

      if(is_first)then
        backspace(20)
        read(20,*) N, i, NZ
        allocate(A(N,N), source = 0.0d0)
        is_first = .false.
      else
        in = in + 1
        backspace(20)
        read(20,*) i, j, val
        A(i,j) = val
        A(j,i) = val
      endif
      if(in == NZ) exit
    enddo
    close(20)
  end subroutine input_from_matrix_market_mm

  subroutine input_from_matrix_market_csr(fin, N, NZ, index, item, A)
    implicit none
    integer(4) :: N, NZ, i, j, in
    integer(4), allocatable :: index(:), item(:)
    real(8), allocatable :: A(:), tmp(:,:)
    character :: fin*100

    call input_from_matrix_market_mm(fin, N, tmp)

    NZ = 0
    do i = 1, N
      do j = 1, N
        if(dabs(tmp(j,i)) > 1.0d-10) NZ = NZ + 1
      enddo
    enddo

    allocate(index(0:N), source = 0)
    allocate(item(NZ), source = 0)
    allocate(A(NZ), source = 0.0d0)

    in = 0
    do i = 1, N
      do j = 1, N
        if(dabs(tmp(j,i)) > 1.0d-10)then
          in = in + 1
          index(i) = index(i) + 1
          item(in) = j
          A(in) = tmp(j,i)
        endif
      enddo
    enddo

    do i = 1, N
      index(i) = index(i) + index(i-1)
    enddo
  end subroutine input_from_matrix_market_csr
end program main
