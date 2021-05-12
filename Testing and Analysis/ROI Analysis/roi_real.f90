program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image
    double precision, allocatable, dimension(:,:,:,:,:) :: roi
    double precision, dimension(2) :: centrePx, roiCentre
    double precision, dimension(5) :: roiRadius
    double precision, allocatable, dimension(:) :: angle
    double precision :: timeStep
    integer :: i, j, radii, angles, m, startImg, numimg, roiSize, numRoi, numAngles, endImg
    character(200) :: filename, c, roiNumber

    !centrePx(1) = 210.0 ; centrePx(2) = 294.0
    centrePx(1) = 206.0 ; centrePx(2) = 280.0
    roiRadius(1) = 50.0 ; roiRadius(2) = 60.0 ; roiRadius(3) = 70.0 ; roiRadius(4) = 80.0 ; roiRadius(5) = 90.0

    roiCentre(1) = centrePx(1)

    angle = angle*((2D0*3.141592653589793D0)/360D0)

    startImg = 68

    endImg = 178

    numimg = 45

    ! (length of roi box - 1) / 2
    roiSize = 4

    numRoi = 5

    numAngles = 36

    timeStep = 1D-6

    allocate(angle(numAngles))

    do i = 1, numAngles
        angle(i) = (i-1)*(180/numAngles) -90
        angle(i) = angle(i)*((2D0*3.141592653589793D0)/360D0)
    end do

    angle(1) = -0.785398
    angle(2) = 0.785398

    !allocate(roi(-roiSize:roiSize,-roiSize:roiSize,numRoi,numAngles,numimg))
    allocate(roi(-roiSize:roiSize,-roiSize:roiSize,numRoi,numAngles,endImg))
    !do i = 1, numimg
        !if ((i+startImg) .lt. 100) then
            !write(filename,'("../../Blurred Images 1/Image",I3,".txt")') (i+startImg)

        !else
            !write(filename,'("../../Blurred Images 1/Image",I3,".txt")') (i+startImg)    
        !end if
        !open(10+i,file=trim(filename))
    !end do

    do i = startImg, endImg, 2
        if (i .lt. 100) then
            write(filename,'("C:\Users\adam\Downloads\OneDrive_7_11-03-2021\16022021_Q12_Squalane_SUM\16022021_Q12_Squalane_SUM_0"&
            ,I0)') (i-68)+startImg
        else
            write(filename,'("C:\Users\adam\Downloads\OneDrive_7_11-03-2021\16022021_Q12_Squalane_SUM\16022021_Q12_Squalane_SUM_"&
            ,I3)') (i-68)+startImg
        end if
        open(10+(i-67),file=trim(filename))
    end do

    do i = 1, numAngles
        write(roiNumber,'(I2)') i
        open(500+i,file="Data Files/"//trim(roiNumber)//'.csv')
    end do

    ! loops through images
    !do m = 1, numimg
    do m = startImg, endImg, 2
        do i = 1, 420    
            !read(10+m,*) (image(j,i),j=1,420)
            read(10+(m-67),*) (image(j,i),j=1,420)
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

    do angles = 1, 2 !numAngles
        write(500+angles,*) "Time Step,","ROI 1,","ROI 2,","ROI 3,","ROI 4,","ROI 5"
            !do i = 1, numimg
            do i = startImg, endImg, 2
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