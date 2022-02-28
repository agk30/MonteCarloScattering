program lss
    implicit none

    integer :: i
    double precision, dimension(10) :: data_array, background
    double precision :: factor, factor2, factor3, up, down, start1, stop1, start2, stop2

    up = 1.5
    down = 0.5

    data_array = [1.0,2.7,3.0,4.0,4.1,2.7,3.2,2.1,3.6,4.4]
    !background = [0.0,0.7,0.0,0.4,1.1,0.7,0.3,1.2,1.3,0.4]
    background = 1.0*data_array

    call cpu_time(start1)
    do i = 1, 3000000
        call least_squares_fit3(data_array, background, down, up, 3D-8, factor3)
    end do
    call cpu_time(stop1)

    call cpu_time(start2)
    do i = 1, 3000000
        call least_squares_fit(data_array, background, down, up, 3D-8, factor)
    end do
    call cpu_time(stop2)

    print *, factor, factor2, factor3
    print *, (stop1-start1), (stop2-start2)

    contains

    subroutine least_squares_fit3 (data_array, background_array, low_limit, high_limit, tol, factor)
        integer :: i, array_size
        double precision, intent(in) :: low_limit, high_limit, tol
        double precision, intent(out) :: factor
        double precision, intent(in), dimension(:) :: data_array, background_array
        double precision :: x0, x1, x2, x3
        double precision :: sum_sqr_x1, sum_sqr_x2
        double precision, parameter :: golden = 0.38196601125

        array_size = size(data_array)

        x0 = low_limit
        x3 = high_limit

        x1 = x0 + (abs(x3 - x0) * golden)
        x2 = x3 - (abs(x3 - x0) * golden)

        do while (abs(x3-x0) .gt. tol*(abs(x1)+abs(x2)))
            sum_sqr_x1 = 0
            sum_sqr_x2 = 0

            do i = 1, array_size
                sum_sqr_x1 = sum_sqr_x1 + ((data_array(i) - (x1*background_array(i)))**2.0)
                sum_sqr_x2 = sum_sqr_x2 + ((data_array(i) - (x2*background_array(i)))**2.0)
            end do

            if (sum_sqr_x1 .lt. sum_sqr_x2) then
                x3 = x2
                x2 = x1
                x1 = x0 + (abs(x3 - x0) * golden)
            else
                x0 = x1
                x1 = x2
                x2 = x3 - (abs(x3 - x0) * golden)
            end if
        end do

        if (sum_sqr_x1 .lt. sum_sqr_x2) then
            factor = x1
        else
            factor = x2
        end if
    end subroutine

    subroutine least_squares_fit2 (data_array, background_array, low_limit, high_limit, factor)
            
        integer :: i, array_size
        double precision, intent(in) :: low_limit, high_limit
        double precision, intent(out) :: factor
        double precision, intent(in), dimension(:) :: data_array, background_array
        double precision, dimension(3) :: bracket
        double precision :: x_bracket, sum_sqr_x, sum_sqr_2, left_dif, right_dif, gap
        double precision, parameter :: golden = 0.38196601125
        logical :: accepted

        array_size = size(data_array)

        bracket(1) = low_limit
        bracket(3) = high_limit

        bracket(2) = bracket(1) + ((bracket(3) - bracket(1)) * golden)
        x_bracket = bracket(3) - ((bracket(3) - bracket(1)) * golden)

        accepted = .FALSE.

        do while (.not. accepted)
            sum_sqr_2 = 0
            sum_sqr_x = 0
            
            do i = 1, array_size
                sum_sqr_2 = sum_sqr_2 + ((data_array(i) - (bracket(2)*background_array(i)))**2.0)
                sum_sqr_x = sum_sqr_x + ((data_array(i) - (x_bracket*background_array(i)))**2.0)
            end do

            gap = bracket(3) - bracket(1)

            if (x_bracket .gt. bracket(2)) then
                if (sum_sqr_x .lt. sum_sqr_2) then
                    bracket(1) = bracket(2)
                    bracket(2) = x_bracket
                else
                    bracket(3) = x_bracket
                end if
            else
                if (sum_sqr_x .lt. sum_sqr_2) then
                    bracket(2) = x_bracket
                    bracket(3) = bracket(2)
                else
                    bracket(1) = x_bracket
                end if
            end if

            left_dif = bracket(2) - bracket(1)
            right_dif = bracket(3) - bracket(2)

            if (left_dif .gt. right_dif) then
                x_bracket = bracket(1) + ((bracket(3) - bracket(1)) * golden)
            else
                x_bracket = bracket(3) - ((bracket(3) - bracket(1)) * golden)
            end if

            !if (gap .lt. (3E-10*abs(x_bracket - bracket(2)))) then
            if (gap .lt. 1D-14) then
                accepted = .TRUE.
                factor = x_bracket
                !print *, factor, sum_sqr_x, gap
            end if
        end do
    end subroutine

    subroutine least_squares_fit (data_array, background_array, low_limit, high_limit, tolerance, factor)
        implicit none
    
        integer :: i, array_size
        double precision, intent(in) :: low_limit, high_limit, tolerance
        double precision, intent(out) :: factor
        double precision, intent(in), dimension(:) :: data_array, background_array
        double precision, dimension(3) :: bracket
        double precision :: sum_square_lower, sum_square_upper, upper, lower
    
        array_size = SIZE(data_array)
    
        bracket(3) = high_limit
        bracket(1) = low_limit
    
        bracket(2) = bracket(3) - ((bracket(3) - bracket(1)) / 2.0)
    
        sum_square_lower = 0
        sum_square_upper = 0
    
        lower = bracket(2) - ((bracket(2) - bracket(1)) / 2.0)
        upper = bracket(3) - ((bracket(3) - bracket(2)) / 2.0)
    
        do
            sum_square_lower = 0
            sum_square_upper = 0
    
            do i = 1, array_size
                sum_square_lower = sum_square_lower + ((data_array(i) - (lower*background_array(i)))**2.0)
                sum_square_upper = sum_square_upper + ((data_array(i) - (upper*background_array(i)))**2.0)
            end do
    
            if (sum_square_lower .lt. sum_square_upper) then
                if (abs(bracket(3) - bracket(1)) .lt. tolerance) then
                    factor = lower
                   ! print *, factor, sum_square_lower
                    EXIT
                    
                else
                    bracket(3) = bracket(2)
                    bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)
    
                    upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                    lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
                end if
            else
                if (abs(bracket(3) - bracket(1)) .lt. tolerance) then
                    factor = upper
                    !print *, factor, sum_square_upper
                    EXIT
                    
                else
                    bracket(1) = bracket(2)
                    bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)
    
                    upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                    lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
                end if
            end if
        end do
    end subroutine least_squares_fit
end program lss