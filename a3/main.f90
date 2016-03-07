program main
    use grid
    use globals
    use solvers

    implicit none

    double precision, allocatable :: phi(:,:)
    double precision, allocatable :: resi(:)
    double precision    :: del_min = 0.020D0
    double precision    :: conv_tol = 3e-4

    integer             :: isolv

    allocate(x(imax), y(jmax))
    allocate(phi(imax, jmax))
    allocate(A(imax, jmax), mu(imax, jmax))
    allocate(ac(imax, jmax), bc(imax, jmax), cc(imax, jmax), & 
             dc(imax, jmax), ec(imax, jmax), gc(imax, jmax))

    allocate(resi(nint(1e7)))

    phi = 0.5D0
    phi(1:2, :) = 0.0d0
    phi(imax, :) = 0.0d0
    phi(:, jmax) = 0.0d0


    call eval_constants()

    call generategrid(del_min)

    isolv = 2
    call solve_murmcol(phi, resi, isolv, conv_tol)

    call solout(x, y, phi)




    deallocate(x, y)
    deallocate(phi)

end program main
