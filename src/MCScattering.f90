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
 
    implicit none

    integer :: i, j, k, vectorsPerParticle, NumberOfTimePoints, startTimePoint, endTimePoint
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

    character(:), allocatable :: time, date, timeOutput, output_image_path, cwd, proper_path
    character(10) :: runTime_string, runTimeSec_string, runTimeMin_string
    character(17) :: date_time

    !New gaussian values
    !parameters for guassians used in fit
    !double precision, dimension(:), allocatable :: m_s, w_s, std_s
    double precision, dimension(:), allocatable :: m_t, w_t, std_t
    !number of guassians to be used for time and speed calculations
    integer :: n_t

    double precision :: t, x, w_low, w_upper, w_sum
    double precision :: arrivalTime

    integer, dimension(8) :: values

    type(CFG_t) :: input_param

    allocate(m_t(1))
    allocate(w_t(1))
    allocate(std_t(1))

    m_t(1) = 0.0D-6
    w_t(1) = 1.0
    std_t(1) = 1.5D-6
    n_t = 1

    ! TODO put in licensing statement.

    normalRun = .TRUE.

    ! Without calling random seed, random number sequences can often be repeated
    call random_seed

    call cpu_time(startTime)

    ! Assigns all parameters from input files into main program variables

    call load_inputs

    !*****************************************************************************************************
    ! This section prepares a start message then allocate arrays as needed and other necessary parameters
    !*****************************************************************************************************

    call date_time_string(date_time)
    call directory_setup(imagePath, date_time, input_param, linux, output_image_path)

    call getcwd(cwd)
    print "(a)", "Writing to "//cwd//output_image_path

    !allocate(proper_path(len(output_image_path)+2))

    !proper_path = trim('"'//output_image_path//'"')

    open (5, file='outputpath.txt')
    write (5, "(a)") output_image_path

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1

    if (.not. hush) then
        if (.not. fullSim) then
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

    if (fullSim) then
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

    do i = 1, ncyc
        ! For fixing parameters, hopefully modern science can find a better way of doing this
        if (normalRun .eqv. .TRUE.) then
            ! sets the ingoing speed and start time
            !call ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))
            call ingoing_speed_from_Gauss&
            (w_s, m_s, std_s, w_t, m_t, std_t, n_s, n_t, gauss_time, gauss_dist, pulseLength, particleSpeed(1), particleTime(1))

            ! Generates the ingoing direction unit vector of the molecule, along with its start point in space.
            call ingoing_direction(valveRad, valvePos, skimRad, skimPos, colRad, colPos, particleVector(1,:), particleStartPos(1,:))

            ! adds a transverse speed to the molcule as it exits the final apperture.
            call transverse_temp(0D0, 40D0, 40D0, 0.5D0, colPos, (valvePos - colPos), particleTime(1), particleSpeed(1), &
            particleStartPos(1,:), particleVector(1,:))

            ! changes the angle of incidence and starting point of the particle using a rotation matrix
            call rotation(particleVector(1,:), incidenceAngle, particleVector(1,:))
            call rotation(particleStartPos(1,:), incidenceAngle, particleStartPos(1,:))

            ! time taken to travel to the wheel (NOT time of origin for scattered particle)
            tWheel = abs(particleStartPos(1,3) / (particleSpeed(1)*particleVector(1,3)))
            
            ! Establishes scattered particle parameters based on ingoing beam particle
            particleTime(2) = particleTime(1) + tWheel
            particleStartPos(2,1) = particleStartPos(1,1) + (particleVector(1,1)*tWheel*particleSpeed(1))
            particleStartPos(2,2) = particleStartPos(1,2) + (particleVector(1,2)*tWheel*particleSpeed(1))
            particleStartPos(2,3) = 0

            ! Decides whicih scattering regime to simulate
            if (scattering) then
                call random_number(rand1)
                
                ! first case: TD scattering
                if (rand1 .gt. scatterFraction) then
                    ! Obtains Maxwell Boltzmann speed as well as scattered direction
                    call MB_speed(maxSpeed, temp, mass, mostLikelyProbability, particleSpeed(2))
                    call cosine_distribution(cosinePowerTD, particleVector(2,:))
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
    ! prepares smoothed IF image
    open(2000,file=trim(ifPath))

    do i = 1, xPx      
        read(2000,*) (ifinput(i,j),j=1,420)       
    end do

    call sg_array(xPx, zPx, ksize, matrixPath, ifinput, ifoutput)
    call sg_convolve(xPx, zPx, NumberOfTimePoints, image, ifoutput)

    ! writes image arrays out into files if writeimages is set to .true.
    if (writeImages) then
        tStepInt = int(tStep*1D6)
        call write_image(image, xPx, zPx, probeStart, probeEnd, tstep, NumberOfTimePoints, date_time, imagePath)
    end if

    totalTraj = real(ncyc)*real(vectorsPerParticle)

    if (.not. hush) then
        print "(ES8.1E2,a,a)", totalTraj, " ", "Total trajectories"
    end if


end program MCScattering