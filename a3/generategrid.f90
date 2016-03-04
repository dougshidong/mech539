module grid
    contains

    subroutine generategrid(imax, jmax, del_min, x, y)

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
    integer            :: i, j
    integer            :: imax, jmax ! imax must be divisible by 2
    double precision   :: del_min ! Minimum Grid Spacing
    double precision   :: x(1:imax), y(1:jmax) 

    !  Generate distribution of x-coordinates
    x(imax/2)=0.0D0

    !  Loop over from imax/2 to imax.
    do i = (imax/2 + 1), imax
      if (x(i-1) <= 0.50D0) then
        ! If x is on the airfoil surface then use constant delta
        x(i) = x(i-1) + del_min 
      else
        ! If x is downstream then stretch the grid until the right farfield boundary
        x(i) = x(i-1) + (x(i-1) - x(i-2)) * 1.1**1.25D0 
      end if

      if (i.lt.imax) x(imax-i) = -x(i) ! Mirror the grid about imax/2.
    end do

    x=x+20.5D0 ! Move the midpoint of the grid to 20.5

    !  Generate distribution of y-coordinates
    !- Y-Grid

    ! Initialize first and second y coordinates. 
    y(1) = 0.0D0
    y(2) = del_min
    do j = 3, jmax
      ! Stretch the grid until the upper farfield boundary
      y(j) = y(j-1) + (y(j-1) - y(j-2)) * 1.1**1.25D0
    end do

    call gridout(x, y) 
    end subroutine generategrid

    subroutine gridout(x, y)

    double precision, dimension (:) :: x, y
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

    return

201 format(15H TITLE = "GRID")
202 format(16H FILETYPE = FULL)
203 format(12H VARIABLES =,8H "X" "Y")
204 format(5H ZONE)
205 format(2(A,I3),A)
206 format(2(E16.10,2X))

    end subroutine gridout

end module grid
