program find_factor
    implicit none

    double precision, dimension(5) :: test_data, test_background
    double precision :: factor, tolerance, lower, upper

    tolerance = 1E-6

    upper = 1.5
    lower = 0.5

    test_data = [1.0, 2.0, 3.0, 4.0, 5.0]
    test_background = [1.0, 2.0, 3.0, 4.0, 5.0]

    call least_squares_fit(test_data(:), test_background(:), lower, upper, tolerance, factor)

    print *, factor

    contains

        subroutine least_squares_fit (data_array, background_array, low_limit, high_limit, tolerance, factor)
            implicit none
        
            integer :: i, array_size, bracket_index
            double precision, intent(in) :: low_limit, high_limit
            double precision, intent(out) :: factor, tolerance
            double precision, intent(in), dimension(:) :: data_array, background_array
            double precision, dimension(3) :: bracket, sum_square_initial
            double precision :: current_min, prev_min, sum_square_lower, sum_square_upper, upper, lower
        
            array_size = SIZE(data_array)
        
            bracket(3) = high_limit
            bracket(1) = low_limit
        
            bracket(2) = bracket(3) - ((bracket(3) - bracket(1)) / 2.0)
        
            sum_square_lower = 0
            sum_square_upper = 0
        
            prev_min = 50.0
        
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
                    prev_min = current_min
                    current_min = lower
        
                    if (abs(prev_min - current_min) .lt. tolerance) then
                        factor = lower
                        EXIT
                    else
                        bracket(3) = bracket(2)
                        bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)
        
                        upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                        lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
                    end if
                else
                    prev_min = current_min
                    current_min = upper
                    if (abs(prev_min - current_min) .lt. tolerance) then
                        factor = upper
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
end program find_factor