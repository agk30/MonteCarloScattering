include "mathConstants.f90"
include "getInputs.f90"
include "getSpeeds.f90"
include "getDirections.f90"
include "sheetIntersection.f90"
include "imaging.f90"
include "SGArray.f90"
include "testingMods.f90"
include "m_config.f90"

program MCScattering
    use getInputs
    use getSpeeds
    use getDirections
    use mathConstants
    use sheetIntersection
    use imaging
    use sgconv
    use mod_tests
 
    implicit none

    ! Variables concerning input parameters
    integer :: ncyc, ksize, polyOrder, cosinePowerTD, cosinePowerIS
    double precision :: incidenceAngle, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
    double precision :: colPos, skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, &
     pxMmRatio, maxSpeed, gaussDev, massMol, energyTrans, surfaceMass, exitAngle, scatterFraction, scatterIntensity, fLifeTime, &
      captureGateOpen, captureGateClose
    logical :: scattering, testMods, writeImages, fullSim

    integer :: i, j, k, vectorsPerParticle, NumberOfTimePoints,&
     xPx, zPx, startTimePoint, endTimePoint
    integer :: startVector
    double precision :: tWheel, rand1, deflectionAngle
    double precision :: mostLikelyProbability, startTime, endTime, runTime, acceptanceRatio, &
     entryTime, exitTime, totalTraj
    double precision, dimension(3) :: sheetDimensions, sheetCentre
    ! particle vectors, speeds and start times are given in these arrays with (1,:) for ingoing and (2,:)
    ! for scattered for use in do loop
    double precision, dimension(2,3) :: particleVector, particleStartPos
    double precision, dimension(2) :: particleSpeed, particleTime
    ! intersection of planes for top (:,1,:) bottom (:,2,:) front (:,3,:) and back (:,4,:)
    double precision, dimension(2,4,3) :: intersection
    double precision, dimension(:,:,:,:), allocatable :: image
    double precision, dimension(:,:), allocatable :: ifoutput
    logical, dimension(4) :: hitsSheet
    logical :: correctDirection

    ! TODO put in licensing statement.

    ! Without calling random seed, random number sequences can often be repeated
    call random_seed

    call cpu_time(startTime)

    ! Loads parameters from input file into main body of code for use in other functions
    call loadInputs (xPx, zPx, incidenceAngle, ncyc, x0, aMax, aMin, &
     h, s, dist, pulseLength, mass, massMol, energyTrans, surfaceMass, exitAngle, temp, skimPos, valvePos, colPos, &
      skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth,&
       probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering, gaussDev, ksize, polyOrder, testMods,&
        writeImages, fullSim, scatterFraction, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose, &
         cosinePowerTD, cosinePowerIS)

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1

    if (.not. fullSim) then
        print "(a)", "Scattering only"
    end if
    
    if (.not. writeImages) then
        print "(a)", "Image writing disabled"
    end if

    ! allocates the image array, which is shared from the imaging class
    allocate(image(zPx,xPx,NumberOfTimePoints,3))
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
        mostLikelyProbability = MBMostLikely(temp, mass)
        
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
        ! sets the ingoing speed, directional unit vector and start time and point 
        call ingoingSpeed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))
        call ingoingDirection(valveRad, valvePos, skimRad, skimPos, colRad, colPos, particleVector(1,:), particleStartPos(1,:))

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

        if (scattering) then
            call random_number(rand1)

            ! alter this if statement to have a higher or lower fraction of TD vs IS scattering
            ! TODO replace as input variable
            
            ! first case: TD scattering
            if (rand1 .gt. scatterFraction) then
                ! Obtains Maxwell Boltzmann speed as well as scattered direction
                call MBSpeed(maxSpeed, temp, mass, mostLikelyProbability, particleSpeed(2))
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
                call getDeflectionAngle(particleVector(1,:), particleVector(2,:), deflectionAngle)
                call softSphereSpeed(massMol, energyTrans, surfaceMass, particleSpeed(1), deflectionAngle, particleSpeed(2))
            end if
        end if

        ! Loops through ingoing trajectories (j=1) then scattered trajectories (j=2)
        do j = startVector, vectorsPerParticle
            ! Finds coordinates of intersection with sheet planes and whether or not it lies within the sheet
            call getSheetIntersection(particleVector(j,:), particleStartPos(j,:), sheetCentre, sheetDimensions, intersection(j,:,:))
            call withinSheet(intersection(j,:,:), sheetCentre, sheetDimensions, hitsSheet)

            if (ANY(hitsSheet)) then
                ! If any sheet faces are hit, then intersection times are calculated
                call getIntersectionTime(hitsSheet, intersection(j,:,:), particleStartPos(j,:), particleVector(j,:), &
                 particleSpeed(j), particleTime(j), entryTime, exitTime)
                ! Finds corresponding image timepoints for entry and exit times
                call startEndTimePoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, tStep, &
                 startTimePoint, endTimePoint)
                ! Finds where in the sheet the particle is located and writes position to image array

                call getPosInProbe(image(:,:,:,1), NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, particleTime(j), &
                 probeStart, tStep, particleSpeed(j), pxMmRatio, particleVector(j,:), particleStartPos(j,:),&
                  sheetDimensions, testMods, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose)
            end if
        end do
    end do

    call cpu_time(endTime)

    runTime = endTime - startTime

    print "(a,a,F5.2,a,a)", "Compute finished in"," ", runTime," ", "seconds"

    ! convolutes image with a Gaussian blur
    do k = 1, NumberOfTimePoints
        call convim(image(:,:,k,1), xPx, zPx, gaussDev, image(:,:,k,2))
    end do

    ! prepares smoothed IF image
    call sgarray(xPx, zPx, ksize, polyOrder, ifoutput)

    ! TODO put in its own subroutine somewhere, keep out of main body
    ! convolutes image with smoothed IF image
    do i = 1, 420
        do j = 1, 420
            do k = 1, NumberOfTimePoints
                if (image(i,j,k,2) .gt. 0) then
                    image(i,j,k,3) = image(i,j,k,2) * ifoutput(i-50,j)
                end if
            end do
        end do
    end do

    ! TODO is this really needed? makes sure there are no negative values in image. Don't be lazy here
    do i = 1, 420
        do j = 1, 420
            do k = 1, NumberOfTimePoints
                if (image(i,j,k,3) .lt. 0) then
                    image(i,j,k,3) = 0
                end if
            end do
        end do
    end do

    ! writes image arrays out into files if writeimages is set to .true.
    if (writeImages) then
        call writeImage(image, xPx, zPx, NumberOfTimePoints)
    end if

    ! writes angle distribution if testMods is set to .true.
    if (testMods) then
        print *, "writing angles"

        call writeAngleDistribution
    end if

    totalTraj = real(ncyc)*real(vectorsPerParticle)

    print "(ES8.1E2,a,a)", totalTraj, " ", "Total trajectories"

end program MCScattering