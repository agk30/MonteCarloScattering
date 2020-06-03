include "../Constants/mathConstants.f90"
include "../Maths/imaging.f90"

program testImaging
    use mathConstants
    use imaging
    implicit none

    integer :: i, j, NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, yPx
    real(kind=r14), dimension(3) :: particleVector, particleStartPos, sheetDimensions
    real(kind=r14) :: particleSpeed, t0, rand1, rand2, probeStart, tStep, pxMmRatio

    open(unit=11,file='testParticlesUp.txt')

    read(11,*)

    sheetDimensions(1) = 0
    sheetDimensions(2) = 0.004D0
    sheetDimensions(3) = 0.03D0 

    NumberOfTimePoints = 25 ; startTimePoint = 2 ; endTimePoint = 25 ; xPx = 420 ; zPx = 420 ; yPx = 420
    t0 = 0 ; probeStart = 5.0D-5 ; tStep = 1.0D-06 ; pxMmRatio = 0.00025D0

    allocate(image(zPx,xPx,NumberOfTimePoints))
    allocate(image2(yPx,xPx,NumberOfTimePoints))
    
    do i = 1, 100000

        read(11,*) j, particleStartPos(1), particleStartPos(2), particleStartPos(3), particleSpeed, particleVector(1), particleVector(2), particleVector(3)

        !print *, j, particleStartPos(1), particleStartPos(2), particleStartPos(3), particleSpeed, particleVector(1), particleVector(2), particleVector(3)

        call getPosInProbeTest(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, yPx, zPx, t0, probeStart, tStep, particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)

    end do

    call writeImage(xPx, zPx, NumberOfTimePoints)
    call writeImageTest(xPx, yPx, NumberOfTimePoints)

end program testImaging