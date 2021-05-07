program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image
    double precision, allocatable, dimension(:,:,:,:,:) :: roi
    double precision, dimension(2) :: centrePx, roiCentre
    double precision, dimension(5) :: roiRadius
    double precision, allocatable, dimension(:) :: angle
    double precision :: timeStep
    integer :: i, j, radii, angles, m, startImg, numimg, roiSize, numRoi, numAngles
    character(500) :: filename, c, roiNumber

    centrePx(1) = 206.0 ; centrePx(2) = 283.0
    roiRadius(1) = 78.0 ; roiRadius(2) = 98.0 ; roiRadius(3) = 118.0 ; roiRadius(4) = 138.0 ; roiRadius(5) = 158.0

    roiCentre(1) = centrePx(1)

    startImg = 74

    numimg = 53

    ! (length of roi box - 1) / 2
    roiSize = 4

    numRoi = 5

    numAngles = 5

    timeStep = 2D-6

    allocate(angle(numAngles))

    !do i = 1, numAngles
        !angle(i) = (i-1)*(180/numAngles) -90
        !angle(i) = angle(i)*((2D0*3.141592653589793D0)/360D0)
    !end do

    angle(1) = -45
    angle(2) = -22.5
    angle(3) = 0
    angle(4) = 22.5
    angle(5) = 45
    angle = angle*((2D0*3.141592653589793D0)/360D0)

    allocate(roi(-roiSize:roiSize,-roiSize:roiSize,numRoi,numAngles,startImg:startImg + ((numImg*2)-2)))
 
    do i = startImg, ((startImg + (2*numimg))-2), 2
        write(filename,'("C:\Users\Maks\Desktop\PhD\Image sequences for ROI analysis\23042021_Q13_Ingoing Beam_BCKGRNDSUB_AVERAGE\23042021_Q13_Ingoing Beam_BCKGRNDSUB_AVERAGE_",I0.3)') i
        open(10+i,file=trim(filename))
    end do

    do i = 1, numAngles
        write(roiNumber,'(F3.1)') angle(i)
        open(500+i,file="Data Files/Scattered at "//trim(roiNumber)//'.csv')
    end do

    ! loops through images
    do m = startImg, ((startImg + (2*numimg))-2), 2

        do i = 1, 420    
            read(10+m,*) (image(j,i),j=1,420)
        end do

        ! loop through angles
        do angles = 1, numAngles
            ! loop through radii
            do radii = 1, 5
                roiCentre(1) = floor(sin(angle(angles))*roiRadius(radii))
                roiCentre(2) = floor(cos(angle(angles))*roiRadius(radii))
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
        print *, (500+angles)
        write(500+angles,*) "Time Step,","ROI 1,","ROI 2,","ROI 3,","ROI 4,","ROI 5"
            do i = startImg, ((startImg + (2*numimg))-2), 2
                write(500+angles,'(I3)',advance='no') i
                write(500+angles,'(a)',advance='no') ","
                do radii = 1, numRoi
                    write(500+angles,'(ES12.5)',advance='no') SUM(roi(:,:,radii,angles,i))
                    if (radii .lt. numRoi) then
                        write(500+angles,'(a)',advance='no') ","
                    else
                        write(500+angles,'(a)',advance='no') new_line(c)
                    end if
                end do
            end do
    end do

end program roiAnalysis