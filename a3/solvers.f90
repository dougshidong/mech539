module solvers
 
    double precision, allocatable :: A(:, :)
    integer, allocatable :: mu(:, :)
    double precision, dimension(:,:), allocatable :: ac, bc, cc, dc, ec, gc

    contains


    subroutine solve_murmcol(phi, resi, times, iterations, isolv, tol)
    use grid
    use boundaryConditions

    implicit none

    double precision :: phi(:, :)
    double precision :: resi(:)
    integer          :: isolv, iterations
    integer          :: i, j, iter
    double precision :: tol, error, phitemp, dphi, res
    double precision :: times(:), tstart

    double precision, dimension(2:jmax-1) :: xvec, lower, upper, diag, rhs

    error = 9999
    iter = 0
    call cpu_time(tstart)
    do while(error > tol)

        call setLower(phi) 
        call evalAmu(phi)
        call evalCoeff()

        iter = iter + 1

        res = 0.0d0

        if(isolv == 2) then ! Gauss-Seidel
            do j = 2, jmax - 1
                do i = 3, imax - 1
                    phitemp = phi(i,j)
                    phi(i, j) =( -bc(i,j) * phi(i,j+1)   &
                                 -cc(i,j) * phi(i,j-1)   &
                                 -dc(i,j) * phi(i-1,j)   &
                                 -ec(i,j) * phi(i+1,j)   &
                                 -gc(i,j) * phi(i-2,j) ) &
                               /  ac(i,j)
                    dphi = abs(phitemp - phi(i,j))
                    res = max(res, dphi)
                end do ! I Loop
            end do ! J Loop
        endif

        if(isolv == 3) then ! Line Implicit Gauss-Seidel
            do i = 3, imax - 1

                diag  = ac(i, 2:jmax-1)
                lower = cc(i, 2:jmax-1)
                upper = bc(i, 2:jmax-1)
                do j = 2, jmax - 1
                    rhs(j) = -(gc(i, j) * phi(i-2, j) &
                             + dc(i, j) * phi(i-1, j) &
                             + ec(i, j) * phi(i+1, j))
                enddo
                rhs(2) = rhs(2) - cc(i, 2) * phi(i, 1)
                rhs(jmax-1) = rhs(jmax-1) - bc(i, jmax-1) * phi(i, jmax)

                call TDMA(lower, diag, upper, rhs, xvec)
                dphi = maxval(abs(phi(i, 2:jmax-1) - xvec(:)))
                res = max(res, dphi)
                phi(i, 2:jmax-1) = xvec(:)
            enddo
        endif
        
        error = res
        resi(iter) = res
        if(mod(iter, 10000) == 0) write(*,*) iter, res
        call cpu_time(times(iter))
    end do ! While Loop

    iterations = iter

    times(:) = times(:) - tstart

    endsubroutine

    subroutine evalAmu(phi)
        use globals
        use grid
        
        implicit none

        double precision, intent(in)  :: phi(:, :)
        integer                       :: i, j
        double precision :: Min2
        Min2 = M_in**2

        mu = 0
        do j = 2, jmax - 1
            do i = 2, imax - 1
                A(i, j) = (1.0d0 - Min2) &
                          - (gam + 1.0d0) * Min2 / U_in &
                          * (phi(i+1, j) - phi(i-1,j)) / (x(i+1) - x(i-1))
                if(A(i, j) < 0) mu(i,j) = 1
            end do
        end do
    end subroutine


    subroutine evalCoeff()
        use globals
        use grid
        implicit none

        double precision :: dxp, dxc, dxm, dxmm, dx12
        double precision :: dyp, dyc, dym
        double precision :: muij, Aij, muijm, Aijm
        integer          :: i, j

        do j = 2, jmax - 1
            dyp = y(j+1) - y(j)
            dyc = y(j+1) - y(j-1)
            dym = y(j)   - y(j-1)
            do i = 3, imax - 1
                dxp  = x(i+1) - x(i)
                dxc  = x(i+1) - x(i-1)
                dxm  = x(i)   - x(i-1)
                dxmm = x(i)   - x(i-2)
                dx12 = x(i-1) - x(i-2)
                muij = mu(i, j)
                muijm = mu(i-1, j)
                Aij = A(i, j)
                Aijm = A(i-1, j)
                ac(i, j) = - 2.0d0 * (1 - muij) * Aij / (dxc * dxp) &
                           - 2.0d0 * (1 - muij) * Aij / (dxc * dxm) &
                           - 2.0d0                            / (dyc * dyp) &
                           - 2.0d0                            / (dyc * dym) &
                           + 2.0d0 * muijm * Aijm / (dxm * dxmm)

                bc(i, j) =   2.0d0 / (dyc * dyp)

                cc(i, j) =   2.0d0 / (dyc * dym)

                dc(i, j) =   2.0d0 * (1 - muij) * Aij / (dxc * dxm) &
                           - 2.0d0 * muijm * Aijm   / (dxm * dxmm) &
                           - 2.0d0 * muijm * Aijm   / (dxmm * dx12)

                ec(i, j) =   2.0d0 * (1 - muij) * Aij / (dxc * dxp)

                gc(i, j) =   2.0d0 * muijm * Aijm    / (dx12 * dxmm) 
            end do
        end do
    end subroutine


    subroutine TDMA(a, b, c, d, x)
        implicit none

        double precision, dimension(:) :: a, b, c, d, x
        double precision :: m
        integer :: n, k

        ! a is lower diagonal with indices a[2:n]
        ! b is diagonal with indices b[1:n]
        ! c is upper diagonal with indices c[1:n-1]

        n = size(x)

        do k = 2, n
            m = a(k) / b(k-1)
            b(k) = b(k) - m * c(k-1)
            d(k) = d(k) - m * d(k-1)
        enddo

        x(n) = d(n) / b(n)
        do k = n-1, 1, -1
            x(k) = (d(k) - c(k) * x(k+1)) / b(k)
        enddo

    end subroutine

end module solvers
