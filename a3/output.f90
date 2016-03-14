module output
    contains
    subroutine solout(phi, fname, resi, times, iterations)
        use grid

        implicit none

        double precision :: phi(:, :)
        double precision :: resi(:), times(:)
        double precision, allocatable, dimension(:, :) :: u, v, Mach, p, cp
        integer :: igridf, i, j, iterations
        character(len=20) :: fname

        allocate(u(imax,jmax), v(imax,jmax), &
                 Mach(imax,jmax), p(imax,jmax), &
                 cp(imax,jmax))
        call get_M_u_v_p(phi, u, v, Mach, p, cp)

        igridf = 10
        ! Output for Tecplot
!       open(unit = igridf, file = trim(fname))
!       write(igridf, 301)
!       write(igridf, 302)
!       write(igridf, 303)
!       write(igridf, 304)
!       write(igridf,305) ' I =', size(x),', J =', size(y),', F=POINT'
!       do j = 1, size(y) 
!         do i = 1, size(x) 
!          write(igridf,306) x(i), y(j), phi(i,j), u(i,j), v(i,j), &
!                            Mach(i,j), p(i,j), cp(i,j)
!         enddo
!       enddo
!       close(igridf)
       
        ! Output for Python
        open(unit = igridf, file = trim(fname))

!!      write(igridf, *) size(x)
!!      write(igridf, *) size(y)

        do i = 1, size(x) 
          do j = 1, size(y) 
!           write(igridf, *) x(i), y(j), phi(i,j), u(i,j), v(i,j), &
!                            Mach(i,j), p(i,j), cp(i,j)
            write(igridf, *) x(i), y(j), cp(i,j)
          enddo
        enddo
        close(igridf)
        
        open(unit = igridf, file = trim('conv'//fname))

!       write(igridf, *) iterations
        do i = 1, iterations 
            write(igridf, *) resi(i), times(i)
        enddo
        close(igridf)


        deallocate(u,v,Mach,p)
        return

301 format(19H TITLE = "Solution")
302 format(16H FILETYPE = FULL)
303 format(12H VARIABLES =,38H "X" "Y" "PHI" "U" "V" "MACH" "P" "CP")
304 format(5H ZONE)
305 format(2(A,I3),A)
306 format(8(E16.10,2X))

    end subroutine solout

    subroutine get_M_u_v_p(phi, u, v, Mach, p, cp)
        use globals
        use grid

        implicit none

        integer :: i, j
        double precision :: phi(:, :)
        double precision, dimension(:, :) :: u, v, Mach, p, cp
        double precision :: sound

        do i = 1, imax
            do j = 1, jmax
                ! Evaluate U
                if(i == 1) then
                    u(i,j) = U_in + (phi(i+1,j) - phi(i,j)) / (x(i+1) - x(i))
                else if(i == imax) then
                    u(i,j) = U_in + (phi(i,j) - phi(i-1,j)) / (x(i) - x(i-1))
                else
                    u(i,j) = U_in + (phi(i+1,j) - phi(i-1,j)) / (x(i+1) - x(i-1))
                endif

                ! Evaluate V
                if(j == 1) then
                    v(i,j) = (phi(i,j+1) - phi(i,j)) / (y(j+1) - y(j))
                else if(j == jmax) then
                    v(i,j) = (phi(i,j) - phi(i,j-1)) / (y(j) - y(j-1))
                else
                    v(i,j) = (phi(i,j-1) - phi(i,j-1)) / (y(j+1) - y(j-1))
                endif
 
                ! Evaluate speed of sound
                sound = sqrt(a02 - 0.50d0 * (gam - 1.0d0) * u(i,j)**2 + v(i,j)**2.0)
!               sound = sqrt(a02 - (gam - 1.0d0) * U_in * (u(i,j) - U_in))
                
                ! Evaluate Mach Number
                Mach(i,j) = sqrt(u(i,j)**2 + v(i,j)**2) / sound

                ! Evaluate Pressure
                p(i,j) = p_in * (1.0d0 - gam * M_in * M_in &
                         * (u(i,j) / U_in - 1.0d0))
                         
                ! Evaluate Pressure Coefficient
                cp(i,j) = (p(i,j) / p_in - 1.0d0) &
                          / (0.50d0 * gam / a_in**2.0d0 * U_in**2.0d0)
                
                if(i == 1) then
                    cp(i,j) = -2.0d0 &
                        * (phi(i+1,j) - phi(i,j)) / (x(i+1) - x(i)) &
                        / U_in
                else if(i == imax) then
                    cp(i,j) = -2.0d0 &
                        * (phi(i,j) - phi(i-1,j)) / (x(i) - x(i-1)) &
                        / U_in
                else
                    cp(i,j) = -2.0d0 &
                        * (phi(i+1,j) - phi(i-1,j)) / (x(i+1) - x(i-1)) &
                        / U_in
                endif

            enddo
        enddo

    end subroutine get_M_u_v_p

end module output
