program SurfaceOutCorrection
    !This program reads in a surface-in measurement a surface-out measurement that's supposed to be subtracted from it, then
    !finds a factor by which the surface-out measurement has to be corrected to minimise the sum of squares of the
    !difference between the two, then, if specified, it corrects the surface-out measurement and outputs a relevant csv
    implicit none

    integer :: i, counter, iostat, num_entries, start_delay, stop_delay, header_lines, delay, num_delays
    double precision, allocatable, dimension(:,:) :: corr_data, SurfaceIn, SurfaceOut
    double precision :: factor, tolerance, lower, upper, signal, ROId, final_angle
    double precision, dimension(2) :: ROI
    character(200) :: header, str, ROI_dist, angle
    logical :: savecorrected

    !These parameters HAVE to be specified by the user
    header_lines = 4                            !Number of header lines in the input files
    start_delay = 82                            !Delay from which to start calculating the correction
    stop_delay = 96                             !Delay at which to finish calculating the correction
    num_delays = ((96-80)/2)+1                  !Number of delays (for 2 us timesteps)
    num_entries = 0 - header_lines              !Number of timpeoints, i.e. images, counted below
    tolerance = 1E-12                           !Tolerance of the least squares fit
    upper = 1.5                                 !Upper limit of the brackets for least squares fit 
    lower = 0.5                                 !Lower limit of the brackets for least squares fit
    
    allocate(SurfaceIn(2,num_delays))
    allocate(SurfaceOut(2,num_delays))
    SurfaceIn = 0
    SurfaceOut = 0
    
    savecorrected = .TRUE.    !TRUE if you want to correct the input SurfaceOut csv and creat a corrected output csv
    !savecorrected = .FALSE.                 

    iostat = 0
    counter = 1

    !Surface-in file
    open(unit=50,file="C:\Users\Maks\Documents\MC OH scattering code\MonteCarloScattering\Testing and Analysis\ROI Analysis\Data Files\ROI 03 Final angle  45.00 SurfaceIn.csv")
    !Surface-out file
    open(unit=60,file="C:\Users\Maks\Documents\MC OH scattering code\MonteCarloScattering\Testing and Analysis\ROI Analysis\Data Files\ROI 03 Final angle  45.00 SurfaceOut.csv")

    !Counts the number of timepoints/images in the input csv file
    do                              
        read(50,*,IOSTAT=iostat)

        if (iostat /= 0) then
            exit
        else
            num_entries = num_entries + 1
        end if
    end do    

    allocate(corr_data(2,num_entries))
    corr_data = 0
   
    rewind(50)

    !Skips the headers of the input files
    do i = 1, header_lines
        read(50,*)
    end do

    !do i = 1, header_lines
     !   read(60,*)
    !end do

    read(60,*)
    read(60,*) ROId, final_angle
    read(60,*)
    read(60,*)

    !Creates the arrays used to find the correction factor through sum of least squares
    do i = 1, num_entries
        read(50,*) SurfaceIn(:,counter)
            if (nint(SurfaceIn(1,counter)) .ge. start_delay) then
                if (nint(SurfaceIn(1,counter)) .le. stop_delay) then
                    counter = counter + 1
                else
                    EXIT
                end if
            end if
    end do

    counter = 1

    do i = 1, num_entries
        read(60,*) SurfaceOut(:,counter)
            if (nint(SurfaceOut(1,counter)) .ge. start_delay) then
                if (nint(SurfaceOut(1,counter)) .le. stop_delay) then
                    counter = counter + 1
                else
                    EXIT
                end if
            end if
    end do

    !Calls the subroutine that does the least squares analysis
    call least_squares_fit(SurfaceIn(2,:), SurfaceOut(2,:), lower, upper, tolerance, factor)

    print *, factor

    !Multiplies the signals in the input surface-out file and creates a new, corrected surface-out csv
    if (savecorrected) then
        write(ROI_dist,'(F6.2)') ROId
        write(angle,'(F6.2)') final_angle
        open(70,file="Data Files/Corr ROI "//trim(ROI_dist)//" Final angle "//trim(angle)//'.csv')
        
        rewind(60)

        do i = 1, header_lines
            read(60,'(a)') header
            write(70,'(a)') trim(header)
        end do

        do i = 1, num_entries 
            read(60,*) corr_data(:,i)
            print *, corr_data(1,i), corr_data(2,i)
            delay = int(corr_data(1,i))
            signal = corr_data(2,i)
            signal = signal*factor
                write(str,'(I3,a,ES12.5)') delay,", ",signal
                write(70,'(a)') trim(str)
        end do
    end if
    
    contains

        subroutine least_squares_fit (data_array, background_array, low_limit, high_limit, tolerance, factor)
            implicit none
        
            integer :: i, array_size
            double precision, intent(in) :: low_limit, high_limit
            double precision, intent(out) :: factor, tolerance
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
end program SurfaceOutCorrection