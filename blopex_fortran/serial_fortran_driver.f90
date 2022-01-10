program main
  use iso_c_binding
  use blopex_fortran_driver
  implicit none
  integer(4) :: N, NZ
  integer(4) :: n_eigs, maxit, loglevel
  integer(4), allocatable :: index(:), item(:)
  real(8) :: tol
  real(8), allocatable :: A(:), B(:), eigen_val(:), eigen_vec(:,:)
  character :: finA*100, finB*100
  logical :: is_B
  logical :: is_prec = .false.

  write(*,"(a)")"* blopex fortran driver"

  call get_input_arg(finA, finB, n_eigs, maxit, loglevel, tol, is_B)

  call input_from_matrix_market_csr(finA, N, NZ, index, item, A)

  if(is_B) call input_from_matrix_B(finB, N, B)

  allocate(eigen_val(n_eigs), source = 0.0d0)
  allocate(eigen_vec(N,n_eigs), source = 0.0d0)

  call blopex_lobpcg_solve(N, NZ, index, item, A, n_eigs, maxit, tol, loglevel, &
    & eigen_val, eigen_vec, is_prec, B)

  call output(N, n_eigs, eigen_val, eigen_vec)

contains

  subroutine output(N, n_eigs, eigen_val, eigen_vec)
    implicit none
    integer(4) :: N, n_eigs, i, j
    real(8) :: eigen_val(:), eigen_vec(:,:)
    character :: cnum*4

    call system('if [ ! -d output ]; then (echo "* create output"; mkdir -p output); fi')

    write(*,"(a)")"* eigen value"
    write(*,"(1pe12.5)")eigen_val

    open(20, file = "output/eigen_value.txt", status = "replace")
    write(20,"(i4)")n_eigs
    do i = 1, n_eigs
      write(20, "(i4,1pe12.4)")i, eigen_val(i)
    enddo
    close(20)

    do i = 1, n_eigs
      write(cnum,"(i0)")i
      open(20, file = "output/eigen_vector."//trim(cnum)//".txt", status = "replace")
      write(20,"(i8)")N
      do j = 1, N
        write(20, "(1pe12.4)")eigen_vec(j,i)
      enddo
      close(20)
    enddo
  end subroutine output

  subroutine get_input_arg(finA, finB, n_eigs, maxit, loglevel, tol, is_B)
    implicit none
    integer(4) :: count
    integer(4) :: n_eigs, maxit, loglevel
    real(8) :: tol
    character :: finA*100, finB*100
    logical :: is_B

    count = iargc()
    if(count == 1)then
      call getarg(1, finA)
    elseif(count == 2)then
      call getarg(1, finA)
      call getarg(2, finB)
      is_B = .true.
    else
      stop "please enter input file name of matrix A and matrix B (optional)"
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
    write(*,"(a,l12)")     "* read B  :", is_B
    write(*,"(a,l12)")     "* apply M :", is_prec
  end subroutine get_input_arg

  subroutine input_from_matrix_market_mm(fin, N, A)
    implicit none
    integer(4) :: N, NZ, i, j, in, ierr
    real(8) :: val
    real(8), allocatable :: A(:,:)
    character :: fin*100, ctemp*10
    logical :: is_first

    is_first = .true.
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

  subroutine input_from_matrix_B(fin, N, B)
    implicit none
    integer(4) :: N, i
    real(8), allocatable :: B(:)
    character :: fin*100

    open(20, file = trim(fin), status = "old")
      read(20,*)i
      if(i /= N) stop "input_from_matrix_B: "
      allocate(B(N), source = 0.0d0)
      do i = 1, N
        read(20,*) B(i)
      enddo
    close(20)
  end subroutine input_from_matrix_B
end program main
