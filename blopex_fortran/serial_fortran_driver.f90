program main
  use iso_c_binding
  use blopex_fortran_driver
  implicit none
  integer(4) :: N, NZ
  integer(4) :: n_eigs, maxit, loglevel, n_bc
  integer(4), allocatable :: index(:), item(:), i_bc(:)
  real(8) :: tol
  real(8), allocatable :: A(:), B(:), eigen_val(:), eigen_vec(:,:)
  character :: finA*100, finB*100, finBC*100
  logical :: is_B = .false.
  logical :: is_BC = .false.
  logical :: is_prec = .false.

  write(*,"(a)")"* blopex fortran driver"

  call get_input_arg(finA, finB, finBC, n_eigs, maxit, loglevel, tol, is_B, is_BC)

  if(is_BC)then
    call input_BC(finBC, n_bc, i_bc)
    call input_from_matrix_market_csr_bc(finA, N, NZ, index, item, A, n_bc, i_bc)
    !if(is_B) call input_from_matrix_B_bc(finB, N, B, n_bc, i_bc)
  else
    call input_from_matrix_market_csr(finA, N, NZ, index, item, A)
    if(is_B) call input_from_matrix_B(finB, N, B)
  endif

  allocate(eigen_val(n_eigs), source = 0.0d0)
  allocate(eigen_vec(N,n_eigs), source = 0.0d0)

  if(is_B)then
    call blopex_lobpcg_solve(N, NZ, index, item, A, n_eigs, maxit, tol, loglevel, &
      & eigen_val, eigen_vec, is_prec, B)
  else
    call blopex_lobpcg_solve(N, NZ, index, item, A, n_eigs, maxit, tol, loglevel, &
      & eigen_val, eigen_vec, is_prec)
  endif

  if(is_BC)then
    !call convet_eigval(N, n_bc, i_bc, eigen_vec)
    call output(N, n_eigs, eigen_val, eigen_vec)
  else
    call output(N, n_eigs, eigen_val, eigen_vec)
  endif

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

  subroutine get_input_arg(finA, finB, finBC, n_eigs, maxit, loglevel, tol, is_B, is_BC)
    implicit none
    integer(4) :: count
    integer(4) :: n_eigs, maxit, loglevel, prec
    real(8) :: tol
    character :: finA*100, finB*100, finBC*100
    logical :: is_B, is_BC

    count = iargc()
    if(count == 1)then
      call getarg(1, finA)
    elseif(count == 2)then
      call getarg(1, finA)
      call getarg(2, finBC)
      is_BC = .true.
    elseif(count == 3)then
      call getarg(1, finA)
      call getarg(2, finBC)
      call getarg(3, finB)
      is_BC = .true.
      is_B = .true.
    else
      stop "please enter input file name of matrix A and matrix B (optional)"
    endif

    open(10, file = "setting.dat", status = "old")
      read(10,*)n_eigs
      read(10,*)maxit
      read(10,*)tol
      read(10,*)loglevel
      read(10,*)prec
    close(10)

    if(prec == 1) is_prec = .true.

    write(*,"(a,i12)")     "* n_eigs  :", n_eigs
    write(*,"(a,i12)")     "* maxiter :", maxit
    write(*,"(a,1pe12.5)") "* tol     :", tol
    write(*,"(a,i12)")     "* loglevel:", loglevel
    write(*,"(a,l12)")     "* read BC :", is_BC
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

  subroutine input_from_matrix_market_csr_bc(fin, N, NZ, index, item, A, n_bc, i_bc)
    implicit none
    integer(4) :: N, NZ, i, j, in, n_bc, ip, jp
    integer(4) :: N_new
    integer(4) :: i_bc(:)
    integer(4), allocatable :: index(:), item(:), perm(:)
    real(8), allocatable :: A(:), tmp(:,:)
    logical, allocatable :: is_use(:)
    character :: fin*100

    call input_from_matrix_market_mm(fin, N, tmp)

    allocate(is_use(N), source = .true.)
    do i = 1, n_bc
      in = i_bc(i)
      is_use(in) = .false.
    enddo

    N_new = N - n_bc
    allocate(perm(N), source = 0)
    in = 1
    do i = 1, N
      if(.not. is_use(i)) cycle
      perm(i) = in
      in = in + 1
    enddo

    NZ = 0
    do i = 1, N
      if(.not. is_use(i)) cycle
      do j = 1, N
        if(.not. is_use(j)) cycle
        if(dabs(tmp(j,i)) > 1.0d-10) NZ = NZ + 1
      enddo
    enddo

    allocate(index(0:N_new), source = 0)
    allocate(item(NZ), source = 0)
    allocate(A(NZ), source = 0.0d0)

    in = 0
    do i = 1, N
      if(.not. is_use(i)) cycle
      do j = 1, N
        if(.not. is_use(j)) cycle
        if(dabs(tmp(j,i)) > 1.0d-10)then
          in = in + 1
          ip = perm(i)
          index(ip) = index(ip) + 1
          jp = perm(j)
          item(in) = jp
          A(in) = tmp(j,i)
        endif
      enddo
    enddo

    do i = 1, N_new
      index(i) = index(i) + index(i-1)
    enddo

    N = N_new
  end subroutine input_from_matrix_market_csr_bc

  subroutine input_BC(fin, n_bc, i_bc)
    implicit none
    integer(4) :: n_bc, i, n, ndof, id, idof
    integer(4), allocatable :: i_bc(:)
    real(8) :: val
    character :: fin*100

    open(20, file = trim(fin), status = "old")
      read(20,*)n_bc, ndof
      allocate(i_bc(n_bc), source = 0)
      do i = 1, n_bc
        read(20,*) id, idof, val
        i_bc(i) = ndof*(id-1) + idof
      enddo
    close(20)
  end subroutine input_BC

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
