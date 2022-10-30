include "hwlib/mathConstants.f90"
include "inputs.f90"
include "hwlib/speeds.f90"
include "hwlib/directions.f90"
include "hwlib/sheet_intersection.f90"
include "hwlib/imaging.f90"
include "hwlib/SGArray.f90"
include "hwlib/testingMods.f90"
include "hwlib/m_config.f90"

program MCScattering
    use inputs
    use speeds
    use directions
    use mathConstants
    use sheetIntersection
    use imaging
    use sgconv
    use mod_tests
    use m_config
    use OMP_LIB
 
    implicit none

    integer :: i, j, k, l, vectorsPerParticle, NumberOfTimePoints, startTimePoint, endTimePoint
    integer :: startVector, runTimeMin, tStepInt, max_int
    double precision :: tWheel, rand1, deflectionAngle
    double precision :: mostLikelyProbability, startTime, endTime, runTime, &
     entryTime, exitTime, totalTraj, runTimeSec
    double precision, dimension(3) :: sheetDimensions, sheetCentre
    ! particle vectors, speeds and start times are given in these arrays with (1,:) for ingoing and (2,:)
    ! for scattered for use in do loop
    double precision, dimension(2,3) :: particleVector, particleStartPos
    double precision, dimension(2) :: particleSpeed, particleTime
    ! intersection of planes for top (:,1,:) bottom (:,2,:) front (:,3,:) and back (:,4,:)
    double precision, dimension(2,4,3) :: intersection
    double precision, dimension(:,:,:,:), allocatable :: image
    double precision, dimension(:,:), allocatable :: ifinput, ifoutput
    integer, dimension(:,:), allocatable :: angleSpeedDist
    logical, dimension(4) :: hitsSheet
    double precision, dimension(100,100) :: surface_heatmap
    double precision :: heatmapx, heatmapy

    character(:), allocatable :: time, date, timeOutput, output_image_path, proper_path, raw_path, blur_path, if_path, input_string, parent_path
    character(300) :: cwd, input_string_temp
    character(10) :: runTime_string, runTimeSec_string, runTimeMin_string
    character(17) :: date_time

    double precision :: speed_total, void_dump

    !New gaussian values
    !parameters for guassians used in fit
    !double precision, dimension(:), allocatable :: m_s, w_s, std_s
    double precision, dimension(:), allocatable :: m_t_scatter, w_t_scatter, std_t_scatter
    integer :: n_t_scatter
    !number of guassians to be used for time and speed calculations

    double precision :: t, x, w_low, w_upper, w_sum, rand
    double precision :: arrivalTime
    double precision :: avg_speed_counter

    integer, dimension(8) :: values

    type(CFG_t) :: input_param

    ! Bins
    double precision :: bin_size
    integer :: n_bins
    double precision, dimension(2) :: bin_range
    double precision, dimension(:), allocatable :: bin_list

    integer, dimension(:,:), allocatable :: surface_grid
    double precision :: surface_grid_interval
    integer :: surface_grid_index_x, surface_grid_index_y, half_grid_size
    logical :: x_found, y_found, do_surf_grid

    double precision :: lor_pos_modifier, lor_rand, lor_speed, lor_speed2

    ! cheat beam is where we ignore most molecular beam generating geometry. Basically, this gives us a good looking beam profile at the expense of cheating the system a bit.
    logical :: cheat_beam

    cheat_beam = .TRUE.

    do_surf_grid = .TRUE.
    surface_grid_interval = 1E-4
    half_grid_size = 200

    if (do_surf_grid) then
        allocate(surface_grid(-half_grid_size:half_grid_size,-half_grid_size:half_grid_size))
    end if

    bin_size = 10.0
    bin_range(1) = 0
    bin_range(2) = 3000
    n_bins = int(( bin_range(2) - bin_range(1) ) / bin_size)
    allocate(bin_list(n_bins))

    allocate(m_t_scatter(1))
    allocate(w_t_scatter(1))
    allocate(std_t_scatter(1))

    m_t_scatter(1) = 0.0D-6
    w_t_scatter(1) = 1.0
    std_t_scatter(1) = 1.5D-6
    n_t_scatter = 1

    ! TODO put in licensing statement.

    normalRun = .TRUE.

    ! Without calling random seed, random number sequences can often be repeated
    call random_seed

    call cpu_time(startTime)

    call load_inputs(input_param)

    call get_command_argument(1, input_string_temp)
    input_string = trim(input_string_temp)

    print *, input_string

    ! Assigns all parameters from input files into main program variables

    !*****************************************************************************************************
    ! This section prepares a start message then allocate arrays as needed and other necessary parameters
    !*****************************************************************************************************

    call date_time_string(date_time)

    call directory_setup(imagePath, date_time, linux, output_image_path, raw_path, blur_path, if_path, input_string, parent_path)

    call CFG_write(input_param, trim(parent_path)//"/input_values.cfg", .FALSE., .FALSE.)

    ! caution when using getcwd:
    ! ifort handles this function just fine if you use a dynamic and unallocated string
    ! but gfortran throws a fit if you try this, must use a properly allocated string
    call getcwd(cwd)
    cwd = trim(cwd)
    print "(a)", "Writing to "//trim(parent_path)

    !allocate(proper_path(len(output_image_path)+2))

    !proper_path = trim('"'//output_image_path//'"')

    open (5, file='outputpath.txt')
    write (5, "(a)") output_image_path

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1

    if (.not. hush) then
        if (.not. show_beam) then
            print "(a)", "Scattering only"
        end if

        if (.not. writeImages) then
            print "(a)", "Image writing disabled"
        end if
    end if

    ! allocates the image array, which is shared from the imaging class
    allocate(image(zPx,xPx,NumberOfTimePoints,3))
    allocate(ifinput(zPx,xPx))
    allocate(ifoutput(zPx,xPx))

    image = 0
    
    ! Establishes array for sheet position and dimensions
    sheetCentre = 0D0
    sheetCentre(3) = sheetCentreZ

    sheetDimensions(1) = 0D0
    sheetDimensions(2) = halfSheetHeight*2D0
    sheetDimensions(3) = sheetWidth

    if (scattering) then
        ! Calculates probability of most probable speed at given temperature for use in thermal desorption subroutines
        mostLikelyProbability = MB_most_likely(temp, mass)
        ! Sets the number of loops in later do loop depending on the number of vectors per particle
        vectorsPerParticle = 2
    else
        vectorsPerParticle = 1
    end if

    if (show_beam) then
        startVector = 1
    else
        startVector = 2
    end if

    if (.not. hush) then
        print "(a)", "Starting compute"
    end if

    ! for prevention of int overflows when using very large times of origin
    max_int = huge(startTimePoint)

    !*****************************************************************************************************
    ! Scattering calculations begin here
    !*****************************************************************************************************

    !$OMP PARALLEL SHARED(image) PRIVATE(particleTime, particleSpeed, particleStartPos, particleVector)

    !$OMP DO
    do i = 1, ncyc
        if (normalRun .eqv. .TRUE.) then
            ! sets the ingoing speed and start time
            !call ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))
            call ingoing_speed_from_Gauss&
            (w_s, m_s, std_s, w_t, m_t, std_t, n_s, n_t, gauss_time, gauss_dist, pulseLength, particleSpeed(1), particleTime(1), time_offset)
            
            !call random_gauss_speed(2040D0, 150D0, particleSpeed(1))
            !call gaussian_distribution(0D0, 4D-6, particleTime(1), particleTime(2))
            !particleTime(1) = particleTime(1) + 22D-6
            !call random_number(rand)
            !particleTime(1) = 0

            ! Generates the ingoing direction unit vector of the molecule, along with its start point in space.
            if (.not. cheat_beam) then
                call ingoing_direction(valveRad, valvePos, skimRad, skimPos, colRad, colPos, particleVector(1,:), particleStartPos(1,:))

                if (trans_speed_modify) then
                    ! adds a transverse speed to the molcule as it exits the final apperture.
                    call transverse_speed(trans_gauss_mean, trans_gauss_sigma, trans_lor_gamma, l_g_fraction, colPos, (valvePos - colPos), particleTime(1), particleSpeed(1), particleStartPos(1,:), particleVector(1,:))
                end if

            else

                particleStartPos(1,:) = 0
                particleStartPos(1,3) = valvePos

                particleVector(1,:) = 0
                particleVector(1,3) = -1

                call transverse_speed(trans_gauss_mean, trans_gauss_sigma, trans_lor_gamma, l_g_fraction, colPos, (valvePos - colPos), particleTime(1), particleSpeed(1), particleStartPos(1,:), particleVector(1,:))
                call random_number(rand1)

                call rotation_z(particleVector(1,:), (90*rand1), particleVector(1,:))

            end if

            if (incidenceAngle .ne. 0) then
                ! changes the angle of incidence and starting point of the particle using a rotation matrix
                call rotation(particleVector(1,:), incidenceAngle, particleVector(1,:))
                call rotation(particleStartPos(1,:), incidenceAngle, particleStartPos(1,:))
            end if

            speed_total = speed_total + particleSpeed(1)

            ! time taken to travel to the wheel (NOT time of origin for scattered particle)
            tWheel = abs(particleStartPos(1,3) / (particleSpeed(1)*particleVector(1,3)))
            
            ! Establishes scattered particle parameters based on ingoing beam particle
            particleTime(2) = particleTime(1) + tWheel
            particleStartPos(2,1) = particleStartPos(1,1) + (particleVector(1,1)*tWheel*particleSpeed(1))
            particleStartPos(2,2) = particleStartPos(1,2) + (particleVector(1,2)*tWheel*particleSpeed(1))
            particleStartPos(2,3) = 0

            if (do_surf_grid) then

                x_found = .FALSE.
                y_found = .FALSE.

                inner1: do l = -half_grid_size, half_grid_size
                    if ((particleStartPos(2,1) .gt. l*surface_grid_interval) .and. (particleStartPos(2,1) .lt. (l+1)*surface_grid_interval)) then
                        surface_grid_index_x = l
                        x_found = .TRUE.
                        EXIT inner1
                    end if
                end do inner1

                inner2: do l = -half_grid_size, half_grid_size
                    if ((particleStartPos(2,2) .gt. l*surface_grid_interval) .and. (particleStartPos(2,2) .lt. (l+1)*surface_grid_interval)) then
                        surface_grid_index_y = l
                        y_found = .TRUE.
                        EXIT inner2
                    end if
                end do inner2

                if (x_found .and. y_found) then
                    surface_grid(surface_grid_index_x, surface_grid_index_y) = surface_grid(surface_grid_index_x, surface_grid_index_y) + 1
                end if
            end if
            ! Decides whicih scattering regime to simulate
            if (scattering) then
                call random_number(rand1)
                
                ! first case: TD scattering
                if (rand1 .gt. scatterFraction) then
                    ! Obtains Maxwell Boltzmann speed as well as scattered direction

                    if (MB_scatter_speed) then
                        call MB_speed(maxSpeed, temp, mass, mostLikelyProbability, particleSpeed(2))
                    else
                        call ingoing_speed_from_Gauss&
                        (w_s_scatter, m_s_scatter, std_s_scatter, w_t_scatter, m_t_scatter, std_t_scatter, n_s_scatter, n_t_scatter, gauss_time_scatter, gauss_dist_scatter, pulseLength, particleSpeed(2), void_dump, time_offset_scatter)
                    end if

                    !bin_list(find_bin_index(particleSpeed(2), bin_range, bin_size, n_bins)) = bin_list(find_bin_index(particleSpeed(2), bin_range, bin_size, n_bins)) + 1
                    !if (particleSpeed(2) .lt. 0) then
                    !    print *, particleSpeed(2)
                    !end if

                    call cosine_distribution(cosinePowerTD, particleVector(2,:))
                    avg_speed_counter = avg_speed_counter + particleSpeed(2)
                ! second case: IS scattering
                else 
                    correctDirection = .false.

                    ! rejects parrticles not scattering in positive z-direction from surface
                    do while (correctDirection .eqv. .false.)
                        ! sets impulsive scattering direction based on some cosine distribution in IS subroutine
                        call cosine_distribution(cosinePowerIS, particleVector(2,:))
                        ! rotates scattered vector about the y-axis (this may not respresent scattered distribution properly)
                        call rotation(particleVector(2,:), exitAngle, particleVector(2,:))

                        if (particleVector(2,3) .gt. 0) then
                            correctDirection = .true.
                        end if
                    end do

                    ! sets IS speed based on scattered direction using soft sphere model
                    call deflection_angle(particleVector(1,:), particleVector(2,:), deflectionAngle)
                    call soft_sphere_speed(massMol, energyTrans, surfaceMass, particleSpeed(1), deflectionAngle, particleSpeed(2))
                end if
            end if
        ! Fixed run parameters
        else
            if (fixedIngoingSpeed) then
                particleSpeed(1) = speedIn
            end if

            if (fixedOutgoingSpeed) then
                particleSpeed(2) = speedOut
            end if

            if (fixedStartPos) then
                particleStartPos(1,:) = (/startx,starty,startz/)
            end if

            if (fixedScatterPos) then
                particleStartPos(2,:) = (/scatterx,scattery,scatterz/)
            end if

            if (fixedCreationTime) then
                particleTime(1) = creationTime
            end if

            if (fixedScatterTime) then
                particleTime(2) = scatterTime
            end if
        end if

        ! Loops through ingoing trajectories (j=1) then scattered trajectories (j=2)
        do j = startVector, vectorsPerParticle
            ! Finds coordinates of intersection with sheet planes and whether or not it lies within the sheet
            call sheet_intersection(particleVector(j,:), particleStartPos(j,:), sheetCentre,&
                sheetDimensions, intersection(j,:,:))
            call within_sheet(intersection(j,:,:), sheetCentre, sheetDimensions, hitsSheet)
            if (ANY(hitsSheet)) then
                ! If any sheet faces are hit, then intersection times are calculated
                call intersection_time(hitsSheet, intersection(j,:,:), particleStartPos(j,:), particleVector(j,:), &
                particleSpeed(j), particleTime(j), entryTime, exitTime)
                ! Finds corresponding image timepoints for entry and exit times
                call start_end_timepoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, tStep, &
                startTimePoint, endTimePoint)

                ! if t0 is very large (in order of seconds) it can end up with a integer overflow for the start and end timepoints
                ! to combat this, just have to find the max value you can get out of the int (see huge() function at start of program)
                ! make sure startTimePoint is larger than the lowest possible int, and endTimePoint is less than the maximum possible int
                if ((startTimePoint .gt. -max_int) .and. (endTimePoint .lt. max_int)) then
                    ! Finds where in the sheet the particle is located and writes position to image array
                    call position_in_probe(image(:,:,:,1), startTimePoint, &
                    endTimePoint, xPx, zPx, particleTime(j), &
                    probeStart, tStep, particleSpeed(j), pxMmRatio, particleVector(j,:), particleStartPos(j,:),&
                    sheetDimensions, testMods, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose, surface_z)
                end if
            end if
        end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

    !$OMP BARRIER

    call cpu_time(endTime)
    runTime = endTime - startTime

    if (runTime .lt. 60D0) then
        write (runTime_string, "(F7.2)") runTime
        !write (runTime_string, "(a,a,F7.2,a,a)") "Compute finished in"," ", runTime," ", "seconds"
        timeOutput = "Compute finished in "//trim(adjustl(runTime_string))//" seconds"
    else
        runTimeMin = floor(runTime/60D0)
        runTimeSec = mod(runTime,60D0)
        if (runTimeMin == 1) then
            write (runTimeMin_string, "(I1)") runTimeMin
            write (runTimeSec_string, "(F6.2)") runTimeSec
            !write (runTime_string, "(a,I1,a,F6.2,a)") "Compute finished in ", runTimeMin, " minute and ", runTimeSec, " seconds."
            timeOutput = "Compute finished in "//trim(adjustl(runTimeMin_string))//" minute and "//trim(adjustl(runTimeSec_string))//" seconds."
        else
            write (runTimeMin_string, "(I4)") runTimeMin
            write (runTimeSec_string, "(F6.2)") runTimeSec
            !write (runTime_string, "(a,I4,a,F6.2,a)") "Compute finished in ", runTimeMin, " minutes and ", runTimeSec, " seconds."
            timeOutput = "Compute finished in "//trim(adjustl(runTimeMin_string))//" minutes and "//trim(adjustl(runTimeSec_string))//" seconds."
        end if
    end if

    if (.not. hush) then
        print "(a)", (trim(timeOutput))
    end if

    !*****************************************************************************************************
    ! Image processing begins, followed by writing out image files
    !*****************************************************************************************************

    ! convolutes image with a Gaussian blur
    do k = 1, NumberOfTimePoints
        call convim(image(:,:,k,1), xPx, zPx, gaussDev, image(:,:,k,2))
    end do

     !prepares smoothed IF image
    open(2001,file=trim(ifPath))

    do i = 1, xPx      
        read(2001,*) (ifinput(i,j),j=1,420)       
    end do

    call sg_array(xPx, zPx, ksize, matrixPath, ifinput, ifoutput)
    call sg_convolve(xPx, zPx, NumberOfTimePoints, image, ifoutput)

    ! writes image arrays out into files if writeimages is set to .true.
    if (writeImages) then
        tStepInt = int(tStep*1D6)
        call write_image(image, xPx, zPx, probeStart, probeEnd, tstep, NumberOfTimePoints, date_time, raw_path, blur_path, if_path)
    end if

    totalTraj = real(ncyc)*real(vectorsPerParticle)

    if (.not. hush) then
        print "(ES8.1E2,a,a)", totalTraj, " ", "Total trajectories"
    end if

    !print *, avg_speed_counter/ncyc

    open(1000110010, file="bins.txt")
    write(1000110010,'(ES12.5)') bin_list

    open(unit = 1102020, file = trim(parent_path)//"/surface_grid.txt")

    if (do_surf_grid) then
        do i = -half_grid_size,half_grid_size
            do j = -half_grid_size,half_grid_size
                write(1102020,'(I6,a)',advance='no') surface_grid(i,j)," "
            end do

            write(1102020,*)
        end do

        print *, SUM(surface_grid)
    end if

end program MCScattering