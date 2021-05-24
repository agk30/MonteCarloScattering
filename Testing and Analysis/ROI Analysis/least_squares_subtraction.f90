program lss
    implicit none

    integer :: i, array_size, bracket_index
    double precision, allocatable, dimension(:) :: data_array, background_array
    double precision, dimension(3) :: bracket, sum_square_initial
    double precision :: tolerance, residual, upper, lower
    double precision :: current_min, prev_min, sum_square_lower, sum_square_upper

    allocate(data_array(5))
    allocate(background_array(5))

    data_array = [1.0, 2.0, 3.0, 4.0, 5.0]
    background_array = [1.1, 2.2, 3.0, 2.5, 4.0]

    array_size = SIZE(data_array)

    bracket(3) = 1.5
    bracket(1) = 0.5

    bracket(2) = bracket(3) - ((bracket(3) - bracket(1)) / 2.0)

    tolerance = 1E-4

    sum_square_lower = 0
    sum_square_upper = 0

    prev_min = 50.0

    lower = bracket(2) - ((bracket(2) - bracket(1)) / 2.0)
    upper = bracket(3) - ((bracket(3) - bracket(2)) / 2.0)

    do bracket_index = 1, 3
        do i = 1, array_size
            residual = data_array(i) - (background_array(i)*bracket(bracket_index))

            sum_square_initial(bracket_index) = sum_square_initial(bracket_index) + residual**2.0
        end do
    end do

    print *, sum_square_initial

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
                print *, lower
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

                print *, upper
                EXIT

            else
                bracket(1) = bracket(2)
                bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)

                upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
            end if
        end if
    end do
end program lss