module solvers
 
    double precision, allocatable :: A(:, :)
    integer, allocatable :: mu(:, :)
    double precision, dimension(:,:), allocatable :: ac, bc, cc, dc, ec, gc

    contains


    subroutine solve_murmcol(phi, resi, isolv, tol)
    use grid
    use boundaryConditions

    implicit none

    double precision :: phi(:, :)
    double precision, allocatable:: phiold(:, :)
    double precision :: resi(:)
    integer          :: isolv
    integer          :: i, j, iter
    double precision :: tol, error, phitemp, phigs, dphi, res
    double precision :: w = 1.5

    if(isolv == 1) allocate(phiold(size(phi,1), size(phi,2)))
    error = 9999
    iter = 0
    do while(error > tol)
        call evalCoeff()
        
        iter = iter + 1

        if(isolv == 1) phiold = phi
        res = 0.0d0

        do i = 3, imax - 1
            do j = 2, jmax - 1
                phitemp = phi(i,j)
                if(isolv == 1) then
                    phi(i, j) =( -bc(i,j) * phiold(i,j+1)   &
                                 -cc(i,j) * phiold(i,j-1)   &
                                 -dc(i,j) * phiold(i-1,j)   &
                                 -ec(i,j) * phiold(i+1,j)   &
                                 -gc(i,j) * phiold(i-2,j) ) &
                               /  ac(i,j)
                endif
                if(isolv == 2) then
                    phi(i, j) =( -bc(i,j) * phi(i,j+1)   &
                                 -cc(i,j) * phi(i,j-1)   &
                                 -dc(i,j) * phi(i-1,j)   &
                                 -ec(i,j) * phi(i+1,j)   &
                                 -gc(i,j) * phi(i-2,j) ) &
                               /  ac(i,j)
                endif
                if(isolv == 3) then
                    phigs     =( -bc(i,j) * phi(i,j+1)   &
                                 -cc(i,j) * phi(i,j-1)   &
                                 -dc(i,j) * phi(i-1,j)   &
                                 -ec(i,j) * phi(i+1,j)   &
                                 -gc(i,j) * phi(i-2,j) ) &
                               /  ac(i,j)
                    phi(i, j) = (1.0d0-w) * phi(i,j) + w * phigs
                endif
                dphi = abs(phitemp - phi(i,j))
                res = max(res, dphi)
            end do ! I Loop
        end do ! J Loop

        call setLower(phi) 
        
        error = res
        resi(iter) = res
        if(mod(iter, 1) == 0) write(*,*) iter, res

    end do ! While Loop

    endsubroutine

    subroutine evalAmu(phi)
        use globals
        use grid
        
        implicit none

        double precision, intent(in)  :: phi(:, :)
        integer                       :: i, j

        mu = 0
        do i = 2, imax - 1
            do j = 1, jmax
                A(i, j) = (1.0d0 - M_in * M_in) - (gam + 1.0d0) * M_in * M_in / U_in &
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
        integer          :: i, j

        do i = 3, imax - 1
            dxp = x(i+1) - x(i)
            dxc = x(i+1) - x(i-1)
            dxm = x(i) - x(i-1)
            dxmm = x(i) - x(i-2)
            dx12 = x(i-1) - x(i-2)
            do j = 2, jmax - 1
                dyp = y(j+1) - y(j)
                dyc = y(j+1) - y(j-1)
                dym = y(j) - y(j-1)
                ac(i, j) = - 2.0d0 * (1.0d0 - mu(i, j)) * A(i, j) / (dxc * dxp) &
                           - 2.0d0 * (1.0d0 - mu(i, j)) * A(i, j) / (dxc * dxm) &
                           - 2.0d0                                / (dyc * dyp) &
                           - 2.0d0                                / (dyc * dym) &
                           + 2.0d0 * mu(i-1, j) * A(i-1, j) / (dxm * dxmm)

                bc(i, j) =   2.0d0 / (dyc * dyp)

                cc(i, j) =   2.0d0 / (dyc * dym)

                dc(i, j) =   2.0d0 * (1.0d0 - mu(i, j)) * A(i, j) / (dxc * dxm) &
                           - 2.0d0 * mu(i-1, j) * A(i-1, j)       / (dxm * dxmm) &
                           - 2.0d0 * mu(i-1, j) * A(i-1, j)       / (dxmm * dx12)

                ec(i, j) =   2.0d0 * (1.0d0 - mu(i, j)) * A(i, j) / (dxc * dxp)

                gc(i, j) =   2.0d0 * (mu(i-1, j)) * A(i-1,j)      / (dx12 * dxmm) 
            end do
        end do
    end subroutine



end module solvers
