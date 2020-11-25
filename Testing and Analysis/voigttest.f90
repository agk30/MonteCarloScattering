subroutine lorentzian_distribution(speed)
    implicit none

    double precision :: probability, speed, gamma, rand1

    call random_number(rand1)

    gamma = 40.0D0

    speed = gamma*tan(3.1415*(rand1-0.5D0))

end subroutine lorentzian_distribution

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

program voigttest
    implicit none

    double precision :: mean, sigma, speed, gaussSpeed, gaussSpeed2, rand
    integer, dimension(-200:200) :: array
    integer :: i, j

    call random_seed()

    array = 0
    mean = 0
    sigma = 40

    do j = 1, 10000000

        call random_number(rand)

        if (rand .gt. 0) then

            call lorentzian_distribution(speed)

        else

            call gaussian_distribution(mean, sigma, speed, gaussSpeed2)

        end if
        !print *, speed

        !speed = gaussSpeed
        !print *, speed



        do i = -200, 200
            
            if ((speed .gt. i) .and. (speed .lt. (i+1))) then

                array(i) = array(i) + 1

                EXIT

            end if

        end do

    end do

    array = array / 2.0

    do i = -200, 200
        print *, array(i)
    end do

end program voigttest