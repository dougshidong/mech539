module boundaryConditions
    contains

    function dydx(x)

        implicit none

        double precision :: x, dydx

        if(x < 20.0D0 .OR. x > 21.0D0) then
            dydx = 0.0D0
        else
            dydx = 0.08D0 * (-4.0D0 * x + 82)
        end if

    end function dydx

    subroutine setLower(phi)
        use grid
        use globals

        implicit none

        double precision :: phi(:, :)
        integer          :: i

        do i = 3, imax - 1
            phi(i, 1) = U_in * dydx(x(i)) * (y(2) - y(1)) + phi(i, 2)
        end do
    end subroutine setLower


end module boundaryConditions
