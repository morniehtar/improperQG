module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: a = 0.d0, b = 1.d0
    !Integration step
    real(8), parameter :: step = 1d-7
    !Threads number
    integer, parameter :: thr = 4
end module config

module com
    implicit none
    public
    real(8) :: xx
end module com

program output
    use config
    use omp_lib
    use com
    implicit none

    abstract interface
        real(8) function func(x)
            implicit none
            real(8), intent(in) :: x
        end function func
    end interface

    procedure(func), pointer :: fptr1, fptr2

    fptr1 => func1
    fptr2 => func2
    print *, 2*QG10(fptr1, 0.d0, 2.d0, step) + QG10(fptr2, 2.d0, 3.d0, step)
contains

    pure real(8) function func1(x)
        implicit none
        real(8), intent(in) :: x
        func1 = sin(x)/x
    end function func1

    pure real(8) function func2(x)
        implicit none
        real(8), intent(in) :: x
        func2 = (sin(x)+cos(x))/x
    end function func2

    real(8) function QG10(fptr, min, max, step)
        implicit none

        procedure(func), pointer, intent(in) :: fptr
        real(8), intent(in) :: min, max, step

        real(8) :: x, half, mid
        real(8) :: tmp, er, next
        real(8), dimension(10) :: a, w
        integer(kind=8) :: i, j, n

        data w /0.066671344308688d0, 0.149451349150581d0, 0.219086362515982d0, &
                0.269266719309996d0, 0.295524224714753d0, 0.295524224714753d0, &
                0.269266719309996d0, 0.219086362515982d0, 0.149451349150581d0, &
                0.066671344308688d0/
        data a /0.97390652851717174d0, 0.86506336668898454d0, 0.67940956829902444d0, &
                0.43339539412924721d0, 0.14887433898163122d0, -0.14887433898163119d0, &
                -0.43339539412924721d0, -0.67940956829902444d0, -0.86506336668898454d0, &
                -0.97390652851717174d0/

        n = (max-min)/step
        half = step/2
        QG10=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) shared(half, n, min, step, a, w, QG10, fptr) num_threads(thr)
            er = 0
            !$omp do reduction(+: QG10)
            do i = 1, n
                x = min + (i-1)*step
                mid = x + half
                do j = 1, 10
                    next = half*w(j)*fptr(a(j)*half+mid) - er
                    tmp = QG10 + next
                    er = (tmp - QG10) - next
                    QG10 = tmp
                end do
            end do
            !$omp end do
        !$omp end parallel

    end function QG10

end program output
