! to integrate the 3 coupled ODEs for rigid body motion


program main

    implicit none

    ! initial conditions
    double precision, dimension(3) :: inertia, torque, r

    ! timestep, time-end
    double precision :: h, t_end, t

    integer :: i, n

    ! simulation controls
    h = 0.01

    t_end = 3.0

    n = int(t_end / h)

    ! initial conditions
    r(1) = 0.01
    r(2) = 2.0
    r(3) = 12.0

    inertia(1) = 1.0
    inertia(2) = 1.0
    inertia(3) = 0.2

    torque(1) = 0.0
    torque(2) = 0.0
    torque(3) = 4.0

    ! integrator
    ! write(*,'(A33)') "omega_x    omega_y    omega_z"
    do i=1, N

        t = i * h

        write(*,'(4(G12.3,2X))') t, r(1), r(2), r(3)

        r = r + rk4(t, r, h)

    end do

contains
    
    pure function f(t, r)
        ! ODE to solve where f is an array of differentials
        double precision, intent(in) :: r(3), t
        double precision :: f(3)

        ! three coupled ODEs
        f(1) = ( torque(1) - (inertia(3) - inertia(2)) * r(2)*r(3) ) / inertia(1)
        f(2) = ( torque(2) - (inertia(1) - inertia(3)) * r(3)*r(1) ) / inertia(2)
        f(3) = ( torque(3) - (inertia(2) - inertia(1)) * r(1)*r(2) ) / inertia(3)

    end function f


    pure function rk4(t, r, h)

        double precision, intent(in) :: r(3), t, h
        double precision :: rk4(3)
        double precision, dimension(3) :: k1, k2, k3, k4

        k1 = h * f(t, r)
        k2 = h * f(t + 0.5 * h, r + 0.5 * k1)
        k3 = h * f(t + 0.5 * h, r + 0.5 * k2)
        k4 = h * f(t + h, r + k3)

        rk4 = (k1 + (2 * k2) + (2 * k3) + k4) / 6

    end function rk4

end program main
