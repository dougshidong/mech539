module boundaryConditions
    contains

    function dydx(x)

        implicit none

        double precision :: x, dydx

        if(x < 20.0d0-1e-10 .OR. x > 21.0d0+1e-10) then
            dydx = 0.0d0
        else
            dydx = 0.08d0 * (-4.0d0 * x + 82.0d0)
        end if

    end function dydx

    subroutine setLower(phi)
        use grid
        use globals

        implicit none

        double precision :: phi(:, :)
        integer          :: i

        do i = 3, imax - 1
            phi(i, 1) = phi(i, 2) - U_in * dydx(x(i)) * (y(2) - y(1)) 
        end do
    end subroutine setLower


end module boundaryConditions
