include "../Constants/mathConstants.f90"
include "../Maths/imaging.f90"
include "../Maths/sheetIntersection.f90"

program testImaging
    use mathConstants
    use imaging
    use sheetIntersection
    implicit none

    integer :: i, j, NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, yPx
    real(kind=r14), dimension(3) :: particleVector, particleStartPos, sheetDimensions, sheetCentre
    real(kind=r14), dimension(4,3) :: intersection
    real(kind=r14) :: particleSpeed, t0, rand1, rand2, probeStart, probeEnd, tStep, pxMmRatio, entryTime, exitTime
    logical :: intersects
    logical, dimension(4) :: hitsSheet

    intersects = .true.

    open(unit=11,file='testParticlesFrontToBack.txt')

    read(11,*)

    sheetCentre(3) = 0.021D0 
    
    sheetDimensions(1) = 0
    sheetDimensions(2) = 0.004D0
    sheetDimensions(3) = 0.03D0 

    NumberOfTimePoints = 45 ; startTimePoint = 1 ; endTimePoint = 45 ; xPx = 420 ; zPx = 420 ; yPx = 420
    t0 = 0 ; probeStart = 4.0D-5 ; tStep = 1.0D-06 ; pxMmRatio = 0.00025D0
    probeEnd = probeStart + (NumberOfTimePoints*tStep)
    print *, 'probe end is', probeEnd


    allocate(image(zPx,xPx,NumberOfTimePoints))
    allocate(image2(yPx,xPx,NumberOfTimePoints))
    
    do i = 1, 100000

        read(11,*) j, particleStartPos(1), particleStartPos(2), particleStartPos(3), particleSpeed, particleVector(1), particleVector(2), particleVector(3)

        !print *, j, particleStartPos(1), particleStartPos(2), particleStartPos(3), particleSpeed, particleVector(1), particleVector(2), particleVector(3)

        if (intersects) then

            print *, '1'
            call getSheetIntersection(particleVector, particleStartPos, sheetCentre, sheetDimensions, intersection)
            print *, '2'
            call withinSheet(intersection, sheetCentre, sheetDimensions, hitsSheet)

            if (ANY(hitsSheet)) then
                
                print *, '3'
                call getIntersectionTime(hitsSheet, intersection, particleStartPos, particleVector, particleSpeed, t0, entryTime, exitTime)

                print *, '4'
                call startEndTimepoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, tStep, startTimePoint, endTimePoint)
                print *, '5'
                call getPosInProbe(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, t0, probeStart, tStep, particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)

            end if
        else

            !call getPosInProbeTest(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, yPx, zPx, t0, probeStart, tStep, particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)

        end if

    end do

    call writeImage(xPx, zPx, NumberOfTimePoints)
    call writeImageTest(xPx, yPx, NumberOfTimePoints)

end program testImaging