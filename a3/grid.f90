module grid
    double precision, allocatable :: x(:), y(:)
    integer             :: imax, jmax

    contains

    subroutine generategrid()
        use globals

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
        integer            :: i, j, ichord
        double precision   :: chord

        imax = 80
        jmax = 40

        allocate(x(imax), y(jmax))

        !  Generate distribution of x-coordinates
        x(imax/2)=0.0d0

        chord = 1.0d0
        ichord = nint((chord/2)/del_min)
        if(abs((chord/2)/del_min - ichord).gt.1e-4) then
            write(*,*) 'MINIMUM SPACING IS NOT A MULTIPLE OF CHORD LENGTH'
            write(*,*) chord/2/del_min, '   ',ichord
            stop
        end if
        
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

        x = x + 20.5D0 ! Move the midpoint of the grid to 20.5

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
        implicit none
        integer :: igridf, i, j

        igridf = 10
        open(unit = igridf,file = "grid.dat")
        
        write(igridf, 201)
        write(igridf, 202)
        write(igridf, 203)
        write(igridf, 204)
        write(igridf,205) ' I =', imax,', J =', jmax,', F=POINT'
        do j = 1, jmax 
          do i = 1, imax
            write(igridf,206) x(i), y(j)
          enddo
        enddo

        close(igridf)

!       open(unit = igridf,file = "coarsegrid.dat")
!       write(igridf,*) imax
!       write(igridf,*) jmax
!       do i = 1, imax
!           write(igridf,*) x(i)
!       enddo
!       do j = 1, jmax
!           write(igridf,*) y(j)
!       enddo
!       close(igridf)

        return

201 format(15H TITLE = "GRID")
202 format(16H FILETYPE = FULL)
203 format(12H VARIABLES =,8H "X" "Y")
204 format(5H ZONE)
205 format(2(A,I3),A)
206 format(2(E16.10,2X))

    end subroutine gridout

    subroutine gridin(gridmult)
        implicit none
        integer :: gridmult
        integer :: i, j, k, igridf
        integer :: imaxtemp, jmaxtemp
        double precision, allocatable :: xtemp(:), ytemp(:)
        double precision :: dxi

        igridf = 10

        open(unit = igridf,file = "coarsegrid.dat")

        read(igridf,*) imaxtemp
        read(igridf,*) jmaxtemp

        allocate(xtemp(imaxtemp), ytemp(jmaxtemp))

        do i = 1, imaxtemp
            read(igridf,*) xtemp(i)
        enddo

        do j = 1, jmaxtemp
            read(igridf,*) ytemp(j)
        enddo

        close(igridf)

        imax = (imaxtemp - 1) * gridmult + 1
        jmax = (jmaxtemp - 1) * gridmult + 1
        
        allocate(x(imax), y(jmax))

        do i = 1, imaxtemp - 1
            dxi = (xtemp(i+1) - xtemp(i)) / dble(gridmult)
            do k = 1, gridmult
                x((i-1)*gridmult + k) = xtemp(i) + (k - 1) * dxi
            enddo
        enddo
        x(imax) = xtemp(imaxtemp)

        do j = 1, jmaxtemp - 1
            dxi = (ytemp(j+1) - ytemp(j)) / dble(gridmult)
            do k = 1, gridmult
                y((j-1)*gridmult + k) = ytemp(j) + (k - 1) * dxi
            enddo
        enddo
        y(jmax) = ytemp(jmaxtemp)

        deallocate(xtemp, ytemp)

        call gridout()
        return
    end subroutine gridin

end module grid
