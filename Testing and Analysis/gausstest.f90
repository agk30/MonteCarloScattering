subroutine gaussian_distribution(mean, variance, z1, z2)
    implicit none

    double precision :: rand1, rand2, v1, v2, rSquared, z1, z2, mean, variance

    do
        call random_number(rand1)
        call random_number(rand2)

        v1 = (2.0*rand1) - 1.0
        v2 = (2.0*rand2) - 1.0

        rSquared = (v1**2.0) + (v2**2.0)

        if (rSquared .lt. 1) then
            z1 = v1*SQRT((-2.0*log(rSquared))/rSquared)
            z2 = v2*SQRT((-2.0*log(rSquared))/rSquared)

            z1 = mean + variance*z1
            z2 = mean + variance*z2

            EXIT
        end if
    end do

end subroutine

program gausstest
    implicit none

    double precision :: mean, variance, z1, z2, start, finish, time
    integer :: i, j
    double precision, dimension(61) :: array

    mean = 0
    variance = 1
    array = 0

    call cpu_time(start)

    do i = 1, 1000000
        call gaussian_distribution(mean, variance, z1, z2)

        do j = -30, 30
            if ((z1 .gt. ((j-1)*0.1)) .and. ((z1 .lt. (j*0.1)))) then
                array(j+31) = array(j+31) + 1
                EXIT
            end if
        end do

        do j = -30, 30
            if ((z2 .gt. ((j-1)*0.1)) .and. ((z2 .lt. (j*0.1)))) then
                array(j+31) = array(j+31) + 1
                EXIT
            end if
        end do
    end do

    call cpu_time(finish)

    time = finish-start

    do i = 1, 61
        print *, array(i)
    end do

    !print *, z1, z2

end program gausstest