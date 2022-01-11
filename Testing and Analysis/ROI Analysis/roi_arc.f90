program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image, background
    double precision, allocatable, dimension(:,:,:,:,:) :: roi
    double precision, dimension(2) :: centrePx, roiCentre, ProbePx
    double precision, dimension(5) :: roiRadius
    double precision, allocatable, dimension(:) :: angle, angledeg, arcradii, wedges, wedgesdeg
    double precision, allocatable, dimension(:,:,:) :: ArcROI
    double precision :: timeStep, radiiGap, wedgeAngle, hypotenuse, probeangle
    integer :: acceptedArc, n
    integer :: i, j, radii, angles, m, startImg, numimg, roiSize, numRoi, numAngles, mode, numArcs, numWedges, probex, probey, y, x
    character(500) :: filename, c, roiNumber

    !CHOOSE MODE FIRST
    mode = 2

    !Universal variables
    centrePx(1) = 206.0 ; centrePx(2) = 283.0   !Centrepoint of ingoing beam and liquid surface area of interaction
    startImg = 10                               !Discharge-probe delay of the first image in the sequence
    numimg = 28                                 !Number of images in the sequence
    timeStep = 2D-6                             !Timestep between images in the sequence in seconds

    !Mode 1 variables
    roiRadius(1) = 78.0 ; roiRadius(2) = 98.0 ; roiRadius(3) = 118.0 ; roiRadius(4) = 138.0 ; roiRadius(5) = 158.0
    roiCentre(1) = centrePx(1)
    
    ! (length of roi box - 1) / 2
    roiSize = 4
    numRoi = 5
    numAngles = 5

    allocate(angle(numAngles))
    allocate(angledeg(numAngles))
    angle(1) = -45
    angle(2) = -22.5
    angle(3) = 0
    angle(4) = 22.5
    angle(5) = 45
    angledeg = angle
    angle = angle*((2D0*3.141592653589793D0)/360D0)
    
    allocate(roi(-roiSize:roiSize,-roiSize:roiSize,numRoi,numAngles,startImg:startImg + ((numImg*2)-2)))
    roi = 0

    !Mode 2 variables
    ProbePx(1) = 60.0 ; ProbePx(2) = 120.0  !top left corner of the laser sheet probe
    probex = 300.0 ; probey = 130.0         !x and y dimension of the laser sheet probe
    numArcs = 7E0                           !number of arcs, i.e. the distances from the centrepoint 
    radiiGap = (real(probey)/real(numArcs)) !height of the ROI
    wedgeAngle = 15E0                       !angular dimension of the ROI
    numWedges = int(180/wedgeAngle)         !number of wedges, i.e. the divisions across the arcs

    allocate(arcradii(numArcs))
    allocate(wedges(numWedges))
   
    do i = 1, numArcs
        arcradii(i) = radiiGap*i            !distance of each arc from the centrepoint
    end do

    print *, arcradii * 0.00025

    do j = 1, numWedges
        wedges(j) = wedgeAngle*j            !angle of each wedge
    end do

    wedgesdeg = wedges-90                               !saves the wedges in degrees before conversion to radians for file naming purposes
    wedges = wedges-(90+(wedgeAngle/2))                 !this ensures the wedges run from -ve to +ve angles and shifts them by half-measure to centre the ROI on specific angles
    wedges = wedges*((2D0*3.141592653589793D0)/360D0)   !converts from degrees to radians
    
    allocate(ArcROI(numArcs,numWedges,startImg:startImg + ((numImg*2)-2)))
    ArcROI = 0
    acceptedArc = 0

    !Create and open output csv files
    if (mode .eq. 1) then
        do i = 1, numAngles
            write(roiNumber,'(F6.2)') angledeg(i)
            open(500+i,file="Data Files/Final angle "//trim(roiNumber)//'.csv')
        end do
    else

        do i = 1, numWedges
            write(roiNumber,'(F6.2)') wedgesdeg(i)
            open(500+i,file="Data Files/Final angle "//trim(roiNumber)//'.csv')
        end do
    end if

    !Open image file 
    do i = startImg, ((startImg + (2*numimg))-2), 2
        write(filename,'("D:\Scattering Images\2021-11-23_150521\Blurred Images",I0.3)') i
        open(10+i,file=trim(filename))
    end do
    !open(10+startImg, file="C:\Users\adam\Documents\Data\05072021_1_Q11_IB_TOF Profile\05072021_1_Q11_IB_TOF Profile_ChC098")
    do i = 1, 420    
        read(10+startImg,*) (background(j,i),j=1,420)
    end do

    rewind 10+startImg
  
    ! loops through images
    do m = startImg, ((startImg + (2*numimg))-2), 2
        
        !Reads individual pixel's intensities
        do i = 1, 420    
            read(10+m,*) (image(j,i),j=1,420)
        end do

        image = image - background
            
            !Mode 1 is square ROIs at certain angle and distance
            if (mode .eq. 1) then
            
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
            
            !Mode 2 is arc ROIs
            else

                do y = int(probePx(2)), int((probey+probePx(2)))

                    do x = int(probePx(1)), int((probex+probePx(1)))
                        hypotenuse = SQRT(((real(x) - real(centrePx(1)))**2) + ((real(y) - real(centrePx(2)))**2))
                
                        if (hypotenuse .lt. arcradii(numArcs)) then
                
                            do i = numArcs, 1, -1
                
                                if (hypotenuse .lt. arcradii(i)) then
                
                                    if (hypotenuse .gt. arcradii(i-1)) then
                                        acceptedArc = i
                
                                        EXIT
                                    end if
                                end if
                            end do
                
                            probeangle =  acos((real(centrePx(2)) - real(y))/hypotenuse)

                            if (x .lt. centrePx(1)) then
                                probeangle = -probeangle
                            end if
            
                            if ((probeangle .gt. wedges(1)) .and. (probeangle .lt. wedges(numWedges))) then
                
                                do i = 1, numWedges - 1
                
                                    if  (probeangle .gt. wedges(i)) then
                
                                        if (probeangle .lt. wedges(i+1)) then
                                            ArcROI(acceptedArc, i, m) = ArcROI(acceptedArc, i, m) + image(x,y)
                
                                        end if
                                    end if
                                end do
                            end if
                        end if
                    end do
                end do
            end if
    end do
    
    !Output csv from mode 1
    if (mode .eq. 1) then
        do angles = 1, numAngles
            write(500+angles,*) "Time Step,","ROI 1,","ROI 2,","ROI 3,","ROI 4,","ROI 5"
                do i = startImg, ((startImg + (2*numimg))-2), 2
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
    
    !Output csv from mode 2
    else
        do n = 1, numWedges
            write(500+n,*) "Time Step,","ROI 1,","ROI 2,","ROI 3,","ROI 4,","ROI 5,","ROI 6,","ROI 7," 
                do i = startImg, ((startImg + (2*numimg))-2), 2
                    write(500+n,'(I3)',advance='no') i
                    write(500+n,'(a)',advance='no') ", "
                    do j = 1, numArcs
                        write(500+n,'(ES12.5)',advance='no') ArcROI(j,n,i)
                        !print *, SUM(ArcROI(j,:,i))
                        if (j .lt. numArcs) then
                            write(500+n,'(a)',advance='no') ", "
                        else
                            write(500+n,'(a)',advance='no') new_line(c)
                        end if
                    end do
                end do
        end do
    end if
end program roiAnalysis