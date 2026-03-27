program gaussian_generator
    implicit none

    ! Precision and constants
    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: n = 1000

    ! Variables
    real(kind=dp) :: u1, u2, z0, z1, pi
    integer :: i

    ! ---- Executable section ----
    pi = acos(-1.0d0)

    ! Initialize the intrinsic random number generator
    call random_seed()

    print *, "Generating", n, "Gaussian random numbers..."

    open(unit=10, file="gaussian_data.txt", status="replace")

    do i = 1, n
        ! Generate two uniform random numbers in (0, 1]
        call random_number(u1)
        call random_number(u2)

        ! Ensure u1 is not exactly zero to avoid log(0)
        if (u1 < 1.0d-12) u1 = 1.0d-12

        ! Box-Muller transform
        z0 = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
        z1 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * pi * u2)

        write(10,*) z0, z1

        ! Example usage (optional)
        ! print *, "Z0:", z0, " Z1:", z1
    end do

    print *, "Generation complete."

close(10)

end program gaussian_generator
