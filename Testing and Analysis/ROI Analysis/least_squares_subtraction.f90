program find_factor
    implicit none

    integer :: i, counter, iostat, num_entries, start_delay, stop_delay, header_lines
    double precision, dimension(2,9) :: test_data, test_background
    double precision :: factor, tolerance, lower, upper
    character(100) :: str

    tolerance = 1E-12

    upper = 1.5
    lower = 0.5
    counter = 1
    header_lines = 4
    iostat = 0

    start_delay = 80
    stop_delay = 96

    num_entries = 0 - header_lines

    open(unit=50,file="/home/adam/code/MonteCarloScattering/Testing and Analysis/ROI Analysis/ROI 03 Final angle 45.00 SurfaceIn.csv")
    open(unit=60,file="/home/adam/code/MonteCarloScattering/Testing and Analysis/ROI Analysis/ROI 03 Final angle 45.00 SurfaceOut.csv")

    do
        read(50,*,IOSTAT=iostat)

        if (iostat /= 0) then
            exit
        else
            num_entries = num_entries + 1
        end if
    end do

    rewind(50)

    do i = 1, header_lines

        read(50,*)
        read(50,*)
        read(50,*)
        read(50,*)

    end do

    do i = 1, header_lines

        read(60,*)
        read(60,*)
        read(60,*)
        read(60,*)

    end do


    do i = 1, num_entries
        read(50,*) test_data(:,counter)
        if (nint(test_data(1,counter)) .ge. start_delay) then
            if (nint(test_data(1,counter)) .le. stop_delay) then
                counter = counter + 1
            else
                EXIT
            end if
        end if
    end do
    
    counter = 1

    do i = 1, num_entries
        read(60,*) test_background(:,counter)
        if (nint(test_background(1,counter)) .ge. start_delay) then
            if (nint(test_background(1,counter)) .le. stop_delay) then
                counter = counter + 1
            else
                EXIT
            end if
        end if
    end do

    call least_squares_fit(test_data(2,:), test_background(2,:), lower, upper, tolerance, factor)

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
        
            prev_min = 500000.0
        
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