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

    ! Variables concerning input parameters
    integer :: ncyc, ksize, polyOrder, cosinePowerTD, cosinePowerIS, runNumber, xPx, zPx
    double precision :: incidenceAngle, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
    double precision :: colPos, skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, &
     pxMmRatio, maxSpeed, gaussDev, massMol, energyTrans, surfaceMass, exitAngle, scatterFraction, scatterIntensity, fLifeTime, &
      captureGateOpen, captureGateClose
    logical :: scattering, testMods, writeImages, fullSim
    character(200) :: imagePath, matrixPath, ifPath

    integer :: i, j, k, vectorsPerParticle, NumberOfTimePoints, startTimePoint, endTimePoint
    integer :: startVector, runTimeMin, tStepInt
    double precision :: tWheel, rand1, deflectionAngle, perpSpeed, colTime, modifStartTime
    double precision :: mostLikelyProbability, mostLikelyProbabilityPerp, startTime, endTime, runTime, acceptanceRatio, &
     entryTime, exitTime, totalTraj, runTimeSec
    double precision, dimension(3) :: sheetDimensions, sheetCentre, perpVector
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
    logical :: correctDirection, normalRun, linux
    character(200) :: time, date, runNumber_string

    integer, dimension(8) :: values

    type(CFG_t) :: input_param

    ! TODO put in licensing statement.

    normalRun = .TRUE.

    ! Without calling random seed, random number sequences can often be repeated
    call random_seed

    call cpu_time(startTime)

    ! Assigns all parameters from input files into main program variables

    call load_inputs(input_param)

    !Inputs for bash script running
    call CFG_get(input_param, "runNumber", runNumber)

    !Experimental inputs
    call CFG_get(input_param, "skimPos", skimPos)
    call CFG_get(input_param, "valvePos", valvePos)
    call CFG_get(input_param, "colPos", colPos)
    call CFG_get(input_param, "skimRad", skimRad)
    call CFG_get(input_param, "valveRad", valveRad)
    call CFG_get(input_param, "colRad", colRad)
    call CFG_get(input_param, "sheetCentre", sheetCentreZ)
    call CFG_get(input_param, "halfSheetHeight", halfSheetHeight)
    call CFG_get(input_param, "sheetWidth", sheetWidth)
    call CFG_get(input_param, "pulseLength", pulseLength)

    ! Imaging inputs
    call CFG_get(input_param, "pxMmRatio", pxMmRatio)
    call CFG_get(input_param, "probeStart", probeStart)
    call CFG_get(input_param, "probeEnd", probeEnd)
    call CFG_get(input_param, "tStep", tStep)
    call CFG_get(input_param, "gaussDev", gaussDev)
    call CFG_get(input_param, "ksize", ksize)
    call CFG_get(input_param, "polyOrder", polyOrder)
    call CFG_get(input_param, "scattering", scattering)
    call CFG_get(input_param, "fullSim", fullSim)
    call CFG_get(input_param, "testMods", testMods)
    call CFG_get(input_param, "writeImages", writeImages)
    call CFG_get(input_param, "scatterIntensity", scatterIntensity)
    call CFG_get(input_param, "fLifeTime", fLifeTime)
    call CFG_get(input_param, "captureGateOpen", captureGateOpen)
    call CFG_get(input_param, "captureGateClose", captureGateClose)

    ! Mathematical Inputs
    call CFG_get(input_param, "xPx", xPx)
    call CFG_get(input_param, "zPx", zPx)
    call CFG_get(input_param, "incidenceAngle", incidenceAngle)
    call CFG_get(input_param, "cosinePowerTD", cosinePowerTD)
    call CFG_get(input_param, "cosinePowerIS", cosinePowerIS)
    call CFG_get(input_param, "x0", x0)
    call CFG_get(input_param, "aMax", aMax)
    call CFG_get(input_param, "aMin", aMin)
    call CFG_get(input_param, "h", h)
    call CFG_get(input_param, "s", s)
    call CFG_get(input_param, "dist", dist)
    call CFG_get(input_param, "mass", mass)
    call CFG_get(input_param, "massMol", massMol)
    call CFG_get(input_param, "energyTrans", energyTrans)
    call CFG_get(input_param, "surfaceMass", surfaceMass)
    call CFG_get(input_param, "exitAngle", exitAngle)
    call CFG_get(input_param, "temp", temp)
    call CFG_get(input_param, "ncyc", ncyc)
    call CFG_get(input_param, "maxSpeed", maxSpeed)
    call CFG_get(input_param, "scatterFraction", scatterFraction)

    ! File paths go here
    call CFG_get(input_param, "linux", linux)
    call CFG_get(input_param, "imagePath", imagePath)
    call CFG_get(input_param, "matrixPath", matrixPath)
    call CFG_get(input_param, "ifPath", ifPath)

    write(runNumber_string, '(i0)') runNumber

    call directory_setup(imagePath, runNumber_string, input_param, linux)

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1

    if (.not. fullSim) then
        print "(a)", "Scattering only"
    end if
    
    if (.not. writeImages) then
        print "(a)", "Image writing disabled"
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

    print "(a)", "Starting compute"

    do i = 1, ncyc
        if (normalRun .eqv. .TRUE.) then
            ! sets the ingoing speed and start time
            call ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))

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
        else
            call point_source(particleVector(2,:), particleStartPos(2,:), particleSpeed(2), particleTime(2))
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
                ! Finds where in the sheet the particle is located and writes position to image array
                call position_in_probe(image(:,:,:,1), startTimePoint, &
                endTimePoint, xPx, zPx, particleTime(j), &
                probeStart, tStep, particleSpeed(j), pxMmRatio, particleVector(j,:), particleStartPos(j,:),&
                sheetDimensions, testMods, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose)
            end if
        end do
    end do

    call cpu_time(endTime)
    runTime = endTime - startTime

    if (runTime .lt. 60D0) then
        print "(a,a,F7.2,a,a)", "Compute finished in"," ", runTime," ", "seconds"
    else
        runTimeMin = floor(runTime/60D0)
        runTimeSec = mod(runTime,60D0)
        if (runTimeMin == 1) then
            print "(a,I1,a,F6.2,a)", "Compute finished in ", runTimeMin, " minute and ", runTimeSec, " seconds."
        else
            print "(a,I4,a,F6.2,a)", "Compute finished in ", runTimeMin, " minutes and ", runTimeSec, " seconds."
        end if
    end if

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
        call write_image(image, xPx, zPx, probeStart, probeEnd, tstep, NumberOfTimePoints, runNumber, imagePath)
    end if

    totalTraj = real(ncyc)*real(vectorsPerParticle)

    print "(ES8.1E2,a,a)", totalTraj, " ", "Total trajectories"

end program MCScattering