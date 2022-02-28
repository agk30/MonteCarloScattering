subroutine gaussian_speed(outgoing_speed)
implicit none

double precision :: average_speed, sigma, z1, z2
double precision, intent(out) :: outgoing_speed

average_speed = 1900D0
sigma = 100D0

call gaussian_distribution(average_speed, sigma, z1, z2)
outgoing_speed = z1

end subroutine

subroutine gaussian_distribution(mean, sigma, z1, z2)
    implicit none

    double precision :: rand1, rand2, v1, v2, rSquared, z1, z2, mean, sigma

    do
        call random_number(rand1)
        call random_number(rand2)

        v1 = (2.0*rand1) - 1.0
        v2 = (2.0*rand2) - 1.0

        rSquared = (v1**2.0) + (v2**2.0)

        if (rSquared .lt. 1) then
            z1 = v1*SQRT((-2.0*log(rSquared))/rSquared)
            z2 = v2*SQRT((-2.0*log(rSquared))/rSquared)

            z1 = mean + sigma*z1
            z2 = mean + sigma*z2

            EXIT
        end if
    end do
end subroutine

program test

    double precision :: outgoing_speed

    call random_seed

    call gaussian_speed(outgoing_speed)

    print *, outgoing_speed
end program