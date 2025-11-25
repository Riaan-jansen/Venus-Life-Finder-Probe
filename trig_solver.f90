! to solve eqn 5.13 in Murray thesis
!
! tan(ax) = b*cot(cx)
! tan(ax)*tan(cx) - b = 0
! f(x) = 0

! bisector method
! - find interval such that f(a) * f(b) < 0
! - compute midpoint, m
! - eval f(m)
! - if f(a) * f(m) < 0, root is between a and m
! - else it is between m and b
! - repeat until residual is small
! - return that m as the root


! ======================================================================
! double precision function f(x, a, b, c)
!     double precision, intent(in) :: a, b, c
!     double precision, intent(in) :: x

! ! operates on a scalar value of x and returns one real*8
!     f = dtan(x*a) * dtan(x*c) - b
! end function f


double precision function f(x, c)
    double precision, intent(in) :: c ! to look at more parameters - just add arguments!
    double precision, intent(in) :: x        ! variable (angle)

    f = c*dcos(x) + dsin(2*x)

end function f


double precision function f_a6(x, c)
    double precision, intent(in) :: c  ! constants (variables)
    double precision, intent(in) :: x        ! variable (angle)

    ! up to 6th order - can increase
    ! f = sin(2x) + c*cos(x)
    f_a6 = 2*x - 8*(x**3)/6 + 32*(x**5)/120 + &
         & c * (1 - (x**2)/2 + (x**4)/4 - (x**6)/720)

end function f_a6


double precision function f_a2(x, c)
    double precision, intent(in) :: c  ! constants (variables)
    double precision, intent(in) :: x        ! variable (angle)

    ! up to 2nd order
    f_a2 = 2*x + c * (1 - (x**2)/2)

end function f_a2


double precision function f_perp(x, c)
    double precision, intent(in) :: c ! to look at more parameters - just add arguments!
    double precision, intent(in) :: x        ! variable (angle)

    double precision :: r1, r2, R, g, a, b, k

    r1 = 1.0
    r2 = 3.0
    R = 5.0
    g = 9.8
    omega = 12.0

    a = r/(4*g)
    b = ((r1**2)/(r1**2 + r2**2) - 0.5)
    k = ((r2**2)/(r1**2 + r2**2) + 0.5)

    f_perp = (dcos(x)**2)/dsin(2*x) * b/a  + dcos(x)/dsin(2*x) * (k/a) - omega**2

end function f_perp

! ======================================================================


subroutine interval(x, N, c, root)
!
    integer, intent(in) :: N
    double precision, intent(in) :: x(N)
    double precision, intent(in) :: c
    double precision, intent(out) :: root

! apparently you need to redeclare type of func in subroutine
    double precision f  ! without ::

! only calls bisect if the interval covers a root (sign change)
    do i = 1, N - 1
        if (f(x(i), c) * f(x(i+1), c) .lt. 0.0) then
            call bisect(f, x(i), x(i+1), c, root)
        end if
    end do

end subroutine interval


subroutine bisect(f, x1, x2, c, root)
!
    double precision, intent(in) :: x1, x2, c
    double precision, intent(out) :: root  ! the answer
! values of function evaluated at x1, etc, and midpoint
    double precision :: f1, f2, fm, m
    integer :: i
    
    double precision, parameter :: residual = 1.0d-8 
    integer, parameter :: max_iter = 10000
    
    double precision f

    f1 = f(x1, c)
    f2 = f(x2, c)

    if (f1 * f2 .gt. 0.0) then
        print *, "no root in this interval"
        root = x1
        return
    end if

    do i = 1, max_iter
        m = (x1 + x2) / 2.0
        fm = f(m, c)
    ! function is close enough to zero m can be approx the root
    ! not sure about the .or. keep testing
        if (abs(fm) .lt. residual .or. abs(x2 - x1) .lt. residual) then
            root = m
            return
        end if
        if (f1 * fm .lt. 0.0) then
            root = m
            f2 = fm
        else
            root = m
            f1 = fm
        end if
    end do

    root = m

end subroutine bisect


program main


    implicit none

    ! for writing to file and plotting
    character(len=20), parameter :: filename = 'datadump.txt'
    character(len=20), parameter :: python_exe = 'python'

    ! variables
    double precision, parameter :: omega = 30.0  ! rad/s
    double precision, parameter :: r = 1.0       ! half-length of rod
    double precision, parameter :: R1 = 0.1      ! size of masses
    double precision, parameter :: R2 = 0.05

    double precision            :: c  ! , v_z, a, b, d, 

    ! physical constants
    double precision, parameter :: pi = 3.1415926
    double precision, parameter :: g = 9.8

    ! need global assignment
    double precision, parameter :: x0 = -pi/2, xN = pi/2
    double precision, parameter :: dx = 0.05

    ! assigning variable values
    c = ((R2**2 - R1**2)/(R2**2 + R1**2)) * 2*g / (r*omega**2)

    !write(*,'(A,G10.4)') "const. = ", c

    ! plot
    call analytical_comparison(filename, 1)

contains

    subroutine secant(root_avg)
        ! takes root argument to pass value out for plotting
        double precision, intent(inout) :: root_avg

        ! defines the steps and range to search over
        integer, parameter :: N = 10000
        integer, parameter :: iter = 500

        double precision :: x(N)
        double precision :: root(iter)
        
        double precision :: start, end
        integer :: i

        ! initialise the array x to search over
        do i = 1, N
            x(i) = x0 + (xN - x0) * (i - 1) / (N - 1)
        end do
    
        call cpu_time(start)
    
        ! calls bisect if the interval covers a sign change
        root_avg = 0.0
        do i = 1, iter
            call interval(x, N, c, root(i))
            root_avg = root_avg + root(i)
        end do
    ! this still only returns the last found root
    
        root_avg = root_avg / iter
    
        call cpu_time(end)
    
        ! print *, "=========================================="
        ! print *, "=============== root finder =============="
        ! print *, "=========================================="
        ! print *, "root =", root_avg
        ! print *, "=========================================="
        ! print *, "time =", (end-start)

    end subroutine secant


    subroutine analytical_comparison(filename, option)
        ! write the numerical and analytical solutions

        character(len=20), intent(in) :: filename

        integer, intent(in)           :: option

        integer, parameter            :: fileunit = 10

        integer          :: info, start, end, i

        double precision :: f, f_a6, f_a2, f_perp

        double precision :: x, root

        101 FORMAT(4(ES13.4,2X))

        open(fileunit, file=filename, iostat=info)
        if (info.NE.0) stop

        ! finding and writing the root for plotting
        call secant(root)

        write(fileunit, '(4(G12.4))') root, 0, 0, 0

        ! plotting f(theta) over the search interval
        start = int(x0/dx)
        end = int(xN/dx)

        do i=start, end
            x = dx*i
            write(fileunit,101) x, f(x, c), f_a6(x, c), f_a2(x, c)
            ! write(fileunit,101) x, f_perp(x, c)
        end do

        close(fileunit)
    
        call plot(filename, option)

    end subroutine


    subroutine plot(filename, option)
        ! calls the python script for plotting <filename>

        character(len=20), intent(in) :: filename

        integer, intent(in)           :: option

        character(len=128)            :: command

        character(len=4)              :: option_char

        ! integer to string
        write(option_char,'(I0)') option

        ! do plot
        write(command,'(A)') trim(python_exe)//' plot_probe.py "'//  &
        &trim(filename) //'" '//trim(option_char)

        call system(command)

    end subroutine plot

end program main
