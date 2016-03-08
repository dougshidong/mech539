program main
    use grid
    use globals
    use solvers
    use output

    implicit none

    conv_tol = 1e-6
    call question1()
    write(*,*) 'DONE'
end program main

subroutine question1()
    use grid
    use globals
    use solvers
    use output

    implicit none

    double precision, allocatable :: phi(:,:)
    double precision, allocatable :: resi(:), times(:)

    integer :: isolv = 2
    integer :: maxit = nint(1e6)

    character(len=20) :: fname
    fname = 'solq1.dat'

    call gridin(1)

    allocate(A(imax, jmax), mu(imax, jmax))
    allocate(ac(imax, jmax), bc(imax, jmax), cc(imax, jmax), & 
             dc(imax, jmax), ec(imax, jmax), gc(imax, jmax))

    allocate(resi(maxit), times(maxit))

    allocate(phi(imax, jmax))
    phi = 0.00D0

    M_in = 0.83d0
    call eval_constants()

    call solve_murmcol(phi, resi, times, isolv, conv_tol)

    call solout(phi, fname)

    deallocate(x, y)
    deallocate(A, mu)
    deallocate(ac, bc, cc, dc, ec, gc)
    deallocate(resi, times)
    deallocate(phi)
end subroutine question1

subroutine question2()
    use grid
    use globals
    use solvers
    use output

    implicit none

    double precision, allocatable :: phi(:,:)
    double precision, allocatable :: resi(:), times(:)

    integer :: i
    integer :: isolv = 2
    integer :: maxit = nint(1e6)

    character(len=20) :: fname

    call gridin(1)

    allocate(A(imax, jmax), mu(imax, jmax))
    allocate(ac(imax, jmax), bc(imax, jmax), cc(imax, jmax), & 
             dc(imax, jmax), ec(imax, jmax), gc(imax, jmax))

    allocate(resi(maxit), times(maxit))

    allocate(phi(imax, jmax))

    do i = 0, 5
        phi = 0.00D0

        M_in = 0.80d0 + i * 0.02d0
        write(fname, '( "solq2_", I2, ".dat" )' )  80 + i * 2
        call eval_constants()

        call solve_murmcol(phi, resi, times, isolv, conv_tol)

        call solout(phi, fname)
    enddo

    deallocate(x, y)
    deallocate(A, mu)
    deallocate(ac, bc, cc, dc, ec, gc)
    deallocate(resi, times)
    deallocate(phi)
end subroutine question2
