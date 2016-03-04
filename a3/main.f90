program main
    use grid

    implicit none

    double precision, allocatable :: x(:), y(:)
    double precision    :: del_min = 0.020D0
    integer             :: imax = 150, jmax = 66

    allocate(x(imax), y(jmax))

    call generategrid(imax, jmax, del_min, x, y)

    deallocate(x, y)

end program main
