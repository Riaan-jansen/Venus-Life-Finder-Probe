! program to plot trajectories of probes dropped from balloon
! we want to see the distance the probes fall from the balloon
! from 70 km to 50 km

! ######################################################################
!                              TO RUN
!
! gfortran -g -c trajectories.f90 && gfortran -g -o main trajectories.o
! ./main
!
! python trajectory_plot.py trajectory.dat
!
! ######################################################################

module PARAMETERS

    implicit none

    ! ##################################################################
    !                            USER PARAMETERS 
    ! ##################################################################
        
    ! number of probes
    integer, parameter :: n_probes = 20

    ! steady state terminal velocity
    double precision, parameter :: v_z = 3.0

    ! forces from results table (4 480 -1 3)
    double precision, parameter :: fx = -0.001236166
    double precision, parameter :: fy = 0.000181528

    ! steady-state orientation
    double precision, parameter :: cone = 4.0
    double precision, parameter :: pitch = -1.0

    ! probe mass
    double precision, parameter :: mass = 0.03

    ! probe length and width
    double precision, parameter :: length = 0.02
    double precision, parameter :: width = 0.01
    double precision, parameter :: thickness = 0.005

    ! grav acceleration
    double precision, parameter :: g = 9.81

    ! wind strength
    double precision, parameter :: wind_factor = 100.0

    double precision, parameter :: rho = 1.225  ! kg/m3

    double precision, parameter :: c_d = 0.5

    ! ##################################################################

    double precision, parameter :: pi = ACOS(-1d0)

end module PARAMETERS


program main

    use PARAMETERS

    implicit none

    double precision, dimension(3, n_probes) :: position, velocity, force

    double precision :: step, t_end, t_start

    integer :: i, N, tumble_n

    character(len=64), parameter :: filename = "trajectory.dat"

    double precision :: fz, tumble_time, bigstep, z_start, z_end

    fz = mass*g

    ! drop test
    z_start = 70000.0
    z_end = 50000.0

    ! time settings
    ! tumble - defined as transition to steady-state
    tumble_time = terminal_time(v_z)
    step = 0.05
    tumble_n = int(tumble_time / step)

    ! once transitioned to steady-state
    bigstep = 1.0
    t_start = 0.0
    t_end = (z_start - z_end) / v_z

    n = int((t_end - tumble_time) / bigstep) + tumble_n

    ! initialise conditions
    call initialise(position, velocity, force)

    ! integrates over time and writes to file
    call integrate(n)

contains

    function terminal_time(v_z)
        ! input: vz
        ! output: returns REAL time to achieve terminal velocity

        double precision terminal_time

        double precision, intent(in) :: v_z

        terminal_time = v_z / g

    end function terminal_time


    function area(cone, pitch)
        ! input: cone and pitch angle (degrees)
        ! output: relative area for x, y, z

        double precision, intent(in) :: pitch, cone

        double precision, dimension(3) :: area

        integer i

        ! A_x
        area(1) = ABS(length * width * SIN(cone*pi/180))

        area(2) = ABS(length * width * SIN(pitch*pi/180))

        area(3) = ABS(length * width * COS(cone*pi/180) * &
                    & COS(pitch*pi/180))
        
    end function area


    function drag_force(velocity)
        ! input: velocity array (3, number of probes)
        ! output: drag force array (same size)

        double precision, dimension(3, n_probes) :: drag_force

        ! integer, intent(in) :: dim

        double precision, dimension(3, n_probes), intent(inout) :: velocity

        integer i, j

        double precision, dimension(3) :: rel_area

        double precision, dimension(n_probes) :: v_Sq

        double precision :: v_norm, v_maximum

        rel_Area = area(cone, pitch)

        v_maximum = 100.0

        do j = 1, n_probes

            ! for runaway velocity but fixed drag shouldnt be needed
            ! velocity(:,j) = MIN(velocity(:,j), v_maximum)

            v_sq(j) = velocity(1,j)**2 + velocity(2,j)**2 + velocity(3,j)**2

            v_norm = sqrt(v_sq(j))

            if (v_norm .GT. 1.0E-10) then
                drag_force(:,j) = -0.5 * c_d * rho * rel_area(:) * &
                                    & v_sq(j) * velocity(:,j) / v_norm
            else
                drag_force(:,j) = 0.0
            end if

        end do

    end function drag_force


    subroutine initialise(position, velocity, force)
        ! set up probe at origin. introduce initial velocity

        double precision, dimension(3, n_probes), intent(inout) :: position,     &
            & velocity, force

        integer :: seed_num
        integer, allocatable :: seed(:)
            
        call random_seed(size=seed_num)   ! n is processor-dependent
        allocate(seed(seed_num))
        seed = 12                  ! putting arbitrary seed to all elements
        call random_seed(put=seed) ! effectively initializing the RNG with the seed

        position = 0d0
        position(3,:) = z_start

        force(1,:) = 0d0
        force(2,:) = 0d0
        force(3,:) = -mass*g

        call random_number(velocity)

        ! +/- direction with scaled up magnitude
        velocity = velocity - 0.5
        velocity = velocity * wind_factor
        velocity(3,:) = 0d0

    end subroutine initialise


    subroutine verlet(position, velocity, force, step)
        ! verlet integrator for position and velocity of each probe

        double precision, dimension(3, n_probes), intent(inout) :: position,     &
            & velocity, force

        double precision, intent(in) :: step

        double precision, dimension(3, n_probes) :: v_half

        integer :: i
    
        ! eqn of motion
        v_half = velocity + step/2 * force/mass

        position = position + step * v_half

        velocity = v_half + step/2 * force/mass

    end subroutine verlet


    subroutine step_forces(time, force, velocity)
        ! pushing forces towards the steady-state solution

        double precision, dimension(3, n_probes), intent(inout) :: force

        double precision, dimension(3, n_probes), intent(inout) :: velocity

        double precision, intent(in) :: time

        double precision, dimension(3, n_probes) :: drag

        if (time .LE. tumble_time) then

            force(1,:) = fx * time/tumble_time
            force(2,:) = fy * time/tumble_time
            force(3,:) = (fz * time/tumble_time) - mass*g

        else

            drag = drag_force(velocity)
    
            ! write(*,*) "# drag", drag
    
            ! drag force
            ! f = Cd/2*rho*A*v^2
            force(1,:) = drag(1,:)
            
            force(2,:) = drag(2,:)

        end if


    end subroutine step_forces


    subroutine integrate(n)
        ! run integrator for n steps
        ! passing by reference

        integer, intent(in) :: n

        integer :: i, io

        double precision :: time

        character(len=64) :: format

        double precision :: kinetic_energy, potential_energy, total_energy

        ! float for 12 spaces at 3dp, 3 whitespaces, repeated N times
        write(format,'("( ",I0,"(e14.6,2x))")') n_probes*3+1

        open(newunit=io, file=filename, status="replace", action="write")

        write(io,*) "# n_probes: ", n_probes

        write(io,*) "# ", "time    ", "x    ", "y    ", "z    ..."

        ! main integration loop - need some way of changing forces here
        do i = 1, N

            time = i * step

            write(io,format) time, position

            ! write(io,*) "# velocity", velocity
            ! write(io,*) "# forces", force

            call step_forces(time, force, velocity)

            ! pushing forces towards the steady-state solution
            if (i .GE. tumble_n) then

                ! increase step size once everything is settled
                step = bigstep
                force(3,:) = 0.0
                velocity(3,:) = -v_z
            end if
            
            call verlet(position, velocity, force, step)

            kinetic_energy = 0.5 * mass * SUM(velocity**2)
            potential_energy = mass * g * SUM(position(3,:))
            total_energy = kinetic_energy + potential_energy

            !write(io,*) "# Time: ", time, "Total Energy: ", total_energy


        end do

        close(io)

    end subroutine integrate

end program main