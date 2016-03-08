module grid
    double precision, allocatable :: x(:), y(:)
    integer             :: imax = 80, jmax = 40

    contains
    subroutine generategrid(del_min)

    ! ************************************************
    !
    ! 	MECH 539 Project 3 
    !       Generate Grid to Solve TSD Equations
    !
    !       Siva Nadarajah
    !       3 March 2016
    !
    ! ************************************************

    implicit none
    integer            :: i, j, ichord, ifarf
    double precision   :: del_min ! Minimum Grid Spacing
    double precision, parameter :: ggrowth=1.1 ! Exponential Growth
    double precision, parameter :: le=20.0d0, te=21.0d0 ! Leading Edge, Trailing Edge
    double precision, parameter :: xs=0.0d0, xe=50.0d0 ! Leading Edge, Trailing Edge
    double precision, parameter :: ys=0.0d0, ye=50.0d0 ! Leading Edge, Trailing Edge
    double precision   :: chord
    double precision, allocatable :: xspacing(:), yspacing(:)

    !  Generate distribution of x-coordinates
    x(imax/2)=0.0d0

    chord = te - le
    ichord = nint((chord/2)/del_min)
    if(abs((chord/2)/del_min - ichord).gt.1e-4) then
        write(*,*) 'MINIMUM SPACING IS NOT A MULTIPLE OF CHORD LENGTH'
        write(*,*) chord/2/del_min, '   ',ichord
        stop
    end if
    
!   ifarf = imax/2 - ichord
!   allocate(xspacing(ifarf))


!   do i = 1, ifarf
!       xspacing(i) = del_min
!   enddo

    !  Loop over from imax/2 to imax.
    do i = (imax/2 + 1), imax
        if (x(i-1) <= chord/2.0d0) then
          ! If x is on the airfoil surface then use constant delta
          x(i) = x(i-1) + del_min 
        else
          ! If x is downstream then stretch the grid until the right farfield boundary
          x(i) = x(i-1) + (x(i-1) - x(i-2)) * 1.32
        end if

        if (i.lt.imax) x(imax-i) = -x(i) ! Mirror the grid about imax/2.
    end do

    x=x+20.5D0 ! Move the midpoint of the grid to 20.5

    !  Generate distribution of y-coordinates
    !- Y-Grid

    ! Initialize first and second y coordinates. 
    y(1) = 0.0d0
    y(2) = del_min
    do j = 3, jmax
      ! Stretch the grid until the upper farfield boundary
      y(j) = y(j-1) + (y(j-1) - y(j-2)) * 1.160
    end do

    call gridout() 
    end subroutine generategrid

    subroutine gridout()

    integer :: igridf, i, j

    igridf = 10
    open(unit = igridf,file = "grid.dat")
    
    write(igridf, 201)
    write(igridf, 202)
    write(igridf, 203)
    write(igridf, 204)
    write(igridf,205) ' I =', size(x),', J =', size(y),', F=POINT'
    do j = 1, size(y) 
      do i = 1, size(x) 
        write(igridf,206) x(i), y(j)
      enddo
    enddo

    close(igridf)

    return

201 format(15H TITLE = "GRID")
202 format(16H FILETYPE = FULL)
203 format(12H VARIABLES =,8H "X" "Y")
204 format(5H ZONE)
205 format(2(A,I3),A)
206 format(2(E16.10,2X))

    end subroutine gridout

    subroutine solout(x, y, phi)

    double precision, dimension (:) :: x, y
    double precision :: phi(:, :)
    double precision, allocatable, dimension(:, :) :: u, v, Mach, p
    integer :: igridf, i, j

    allocate(u(imax,jmax), v(imax,jmax), Mach(imax,jmax), p(imax,jmax))
    call get_M_u_v_p(phi, u, v, Mach, p)

    igridf = 10
    open(unit = igridf,file = "sol.dat")
    
    write(igridf, 301)
    write(igridf, 302)
    write(igridf, 303)
    write(igridf, 304)
    write(igridf,305) ' I =', size(x),', J =', size(y),', F=POINT'
    do j = 1, size(y) 
      do i = 1, size(x) 
       write(igridf,306) x(i), y(j), phi(i,j), u(i,j), v(i,j), Mach(i,j), p(i,j)
!       write(igridf,306) x(i), y(j), phi(i,j)
      enddo
    enddo

    close(igridf)

    deallocate(u,v,Mach,p)
    return

301 format(19H TITLE = "Solution")
302 format(16H FILETYPE = FULL)
303 format(12H VARIABLES =,33H "X" "Y" "PHI" "U" "V" "MACH" "P")
!303 format(12H VARIABLES =,14H "X" "Y" "PHI")
304 format(5H ZONE)
305 format(2(A,I3),A)
306 format(7(E16.10,2X))
!306 format(3(E16.10,2X))

    end subroutine solout


    subroutine get_M_u_v_p(phi, u, v, Mach, p)
        use globals

        implicit none

        integer :: i, j
        double precision :: phi(:, :)
        double precision, dimension(:, :) :: u, v, Mach, p
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
                if((a02 - (gam - 1.0d0) * u(i,j)**2 + v(i,j)**2.0) < 0) then
                write(*,*) i, j, (a02 - (gam - 1.0d0) * u(i,j)**2 + v(i,j)**2.0)
                endif

                sound = sqrt(a02 - (gam - 1.0d0) * u(i,j)**2 + v(i,j)**2.0)
                sound = sqrt(a02 - (gam - 1.0d0) * U_in * (u(i,j) - U_in))
                ! Evaluate Mach Number
                Mach(i,j) = sqrt(u(i,j)**2 + v(i,j)**2) / sound

                ! Evaluate Pressure
                p(i,j) = p_in * (1.0d0 - gam * M_in * M_in &
                         * (u(i,j) / U_in - 1.0d0))


            enddo
        enddo

    end subroutine get_M_u_v_p

end module grid
