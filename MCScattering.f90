include "Constants/mathConstants.f90"
include "getInputs.f90"
include "Maths/getSpeeds.f90"
include "Maths/getDirections.f90"
include "Maths/sheetIntersection.f90"
include "Maths/imaging.f90"

program MCScattering
    use getInputs
    use getSpeeds
    use getDirections
    use mathConstants
    use sheetIntersection
    use imaging
 
    implicit none

    ! Variables concerning input parameters
    integer :: ncyc
    real(kind=r14) :: x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
    real(kind=r14) :: colPos, skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed
    logical :: scattering

    integer :: i, j, k, vectorsPerParticle, acceptedCounter, totalTraj, NumberOfTimePoints, xPx, zPx, startTimePoint, endTimePoint
    real(kind=r14) :: tWheel
    real(kind=r14) :: t0, mostLikelyProbability, speed, scatteredSpeed, startTime, endTime, runTime, acceptanceRatio, entryTime, exitTime
    real(kind=r14), dimension(3) :: sheetDimensions, sheetCentre, topInter, bottomInter, frontInter, backInter
    ! particle vectors, speeds and start times are given in these arrays with (1,:) for ingoing and (2,:) for scattered for use in do loop
    real(kind=r14), dimension(2,3) :: particleVector, particleStartPos
    real(kind=r14), dimension(2) :: particleSpeed, particleTime
    ! intersection of planes for top (1,:) bottom (2,:) front (3,:) and back (4,:)
    real(kind=r14), dimension(2,4,3) :: intersection
    logical, dimension(4) :: hitsSheet

    acceptedCounter = 0

    ! Without calling random seed, random number sequences can often be repeated
    call random_seed

    call cpu_time(startTime)

    ! Loads parameters from input file into main body of code for use in other functions
    call loadInputs(ncyc, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos, colPos, skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering)

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1

    ! TODO make pixel dimensions an input parameter
    xPx = 420
    zPx = 420

    ! allocates the image array, which is shared from the imaging class
    allocate(image(zPx,xPx,NumberOfTimePoints))

    image = 0
    
    ! Establishes array for sheet position and dimensions
    sheetCentre = 0D0
    sheetCentre(3) = sheetCentreZ

    sheetDimensions(1) = 0D0
    sheetDimensions(2) = halfSheetHeight*2D0
    sheetDimensions(3) = sheetWidth

    if (scattering == .TRUE.) then

        ! Calculates probability of most probable speed at given temperature for use in thermal desorption subroutines
        mostLikelyProbability = MBMostLikely(temp, mass)
        
        ! Sets the number of loops in later do loop depending on the number of vectors per particle
        vectorsPerParticle = 2

    else

        vectorsPerParticle = 1

    end if

    do i = 1, ncyc

        call ingoingSpeed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))
        call ingoingDirection(valveRad, valvePos, skimRad, skimPos, colRad, colPos, particleVector(1,:), particleStartPos(1,:))

        ! time taken to travel to the wheel (NOT time of origin for scattered particle)
        tWheel = abs(particleStartPos(1,3) / particleSpeed(1))
        
        ! Establishes scattered particle parameters based on ingoing beam particle
        particleTime(2) = particleTime(1) + tWheel
        particleStartPos(2,:) = particleStartPos(1,:) + (particleVector(1,:)*tWheel*particleSpeed(1))
        particleStartPos(2,3) = 0
       
        if (scattering == .TRUE.) then

            ! Obtains Maxwell Boltzmann speed as well as scattered direction
            call MBSpeed(maxSpeed, temp, mass, mostLikelyProbability, particleSpeed(2))
            call thermalDesorptionDirection(particleVector(2,:))

        end if

        ! Loops through ingoing trajectories (j=1) then scattered trajectories (j=2)
        do j = 1, vectorsPerParticle

            ! Finds coordinates of intersection with sheet planes and whether or not it lies within the sheet
            call getSheetIntersection(particleVector(j,:), particleStartPos(j,:), sheetCentre, sheetDimensions, intersection(j,:,:))
            call withinSheet(intersection(j,:,:), sheetCentre, sheetDimensions, hitsSheet)

            if (ANY(hitsSheet)) then

                ! If any sheet faces are hit, then intersection times are calculated
                call getIntersectionTime(hitsSheet, intersection(j,:,:), particleStartPos(j,:), particleVector(j,:), particleSpeed(j), particleTime(j), entryTime, exitTime)
                ! Finds corresponding image timepoints for entry and exit times
                call startEndTimePoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, tStep, startTimePoint, endTimePoint)
                ! Finds where in the sheet the particle is located and writes position to image array
                call getPosInProbe(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, particleTime(j), probeStart, tStep, particleSpeed(j), pxMmRatio, particleVector(j,:), particleStartPos(j,:), sheetDimensions)
                
                acceptedCounter = acceptedCounter + 1

            end if

        end do
    

    end do
    
    !print *, SUM(image)

    call writeImage(xPx, zPx, NumberOfTimePoints)

    call cpu_time(endTime)

    runTime = endTime - startTime

    totalTraj = ncyc*vectorsPerParticle
    acceptanceRatio = real(acceptedCounter)/((real(ncyc)*real(vectorsPerParticle)))

    print *, "Finished in", runTime, "seconds"
    print *, totalTraj, "Total trajectories"
    print *, acceptedCounter, "accepted trajectories"
    print "(a, F4.2, a)","  ", acceptanceRatio, " acceptance ratio"

end program MCScattering