program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image
    double precision, allocatable, dimension(:,:,:,:,:) :: roi
    double precision, dimension(2) :: centrePx, roiCentre
    double precision, allocatable, dimension(:) :: angle, angleRadians, roiRadius
    double precision :: pi
    integer :: i, j, radii, angles, m, startDelay, endDelay, numImg, roiSize, numRoi, numAngles, tstep
    character(300) :: fileNumber, c, roiNumber, partialPath

    pi = 3.141592653589793D0

    centrePx(1) = 210.0 ; centrePx(2) = 294.0
    
    roiCentre(1) = centrePx(1)

    startDelay = 74

    endDelay = 178

    !time step between images in microseconds
    tstep = 2

    numImg = (endDelay - startDelay) / tstep

    ! (length of roi box - 1) / 2
    roiSize = 4

    numRoi = 5

    numAngles = 5

    allocate(angle(numAngles))
    allocate(roiRadius(numRoi))
    allocate(angleRadians(numAngles))

    !ROI Angles in degrees
    angle(1) = -45
    angle(2) = -22.5
    angle(3) = 0
    angle(4) = 22.5
    angle(5) = 45

    roiRadius(1) = 50.0
    roiRadius(2) = 60.0
    roiRadius(3) = 70.0
    roiRadius(4) = 80.0
    roiRadius(5) = 90.0

    angleRadians = angle*((2D0*pi)/360D0)

    allocate(roi(-roiSize:roiSize,-roiSize:roiSize,numRoi,numAngles,numImg))

    partialPath = "<Path>"
 
    do i = startDelay, endDelay, tstep
        write(fileNumber,'(I0.3)') i
        open(10+i,file=trim(partialPath)//trim(fileNumber))
    end do

    do i = 1, numAngles
        write(roiNumber,'(I0.3)') angle(i)
        open(500+i,file="Data Files/Scattered at "//trim(roiNumber)//'.csv')
    end do

    ! loops through images
    do m = startDelay, endDelay, tstep

        do i = 1, 420    
            read(10+m,*) (image(j,i),j=1,420)
        end do

        ! loop through angles
        do angles = 1, numAngles
            ! loop through radii
            do radii = 1, 5
                roiCentre(1) = floor(sin(angleRadians(angles))*roiRadius(radii))
                roiCentre(2) = floor(cos(angleRadians(angles))*roiRadius(radii))
                do i = -roiSize, roiSize   
                    do j = -roiSize, roiSize
                        roi(i,j,radii,angles,m) = image((int(centrePx(1)) + int(roiCentre(1)) + i),&
                         (int(centrePx(2)) - int(roiCentre(2) + j)))
                    end do
                end do
            end do
        end do
    end do

    do angles = 1, numAngles
        write(500+angles,*) "Time Step,","ROI 1,","ROI 2,","ROI 3,","ROI 4,","ROI 5"
            do i = startDelay, endDelay, tstep
                write(500+angles,'(I3)',advance='no') i
                write(500+angles,'(a)',advance='no') ", "
                do radii = 1, numRoi
                    write(500+angles,'(ES12.5)',advance='no') SUM(roi(:,:,radii,angles,i))
                    if (radii .lt. numRoi) then
                        write(500+angles,'(a)',advance='no') ", "
                    else
                        write(500+angles,'(a)',advance='no') new_line(c)
                    end if
                end do
            end do
    end do

end program roiAnalysis