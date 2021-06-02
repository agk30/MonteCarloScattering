program ralss
    implicit none

    double precision, dimension(420,420) :: image
    double precision, allocatable, dimension(:,:,:,:,:) :: roi
    double precision, dimension(2) :: centrePx, roiCentre, ProbePx
    double precision, dimension(5) :: roiRadius
    double precision, allocatable, dimension(:) :: angle, angledeg, arcradii, wedges, wedgesdeg
    double precision, allocatable, dimension(:,:,:,:) :: ArcROI
    double precision, allocatable, dimension(:,:) :: factor
    double precision :: timeStep, radiiGap, wedgeAngle, hypotenuse, probeangle, signal
    integer :: acceptedArc, arc_index, wedge_index, delay, chc_index
    integer :: i, j, k, m, startImg, numimg, roiSize, numRoi, numAngles, mode, numArcs, numWedges, probex, probey, y, x
    character(500) :: filename, sin_path, sout_path
    character(4) :: transition_str, surface_str
    integer :: num_entries, start_delay, stop_delay, header_lines, num_delays
    double precision, allocatable, dimension(:,:,:,:) :: corr_data, SurfaceIn, SurfaceOut
    double precision :: tolerance, lower, upper

    ! Retrieve command line arguements
    call get_command_argument(1, sin_path)
    call get_command_arguemnt(2, sout_path)

    call parse_path(sin_path, transition_str, surface_str)

    !CHOOSE MODE FIRST
    mode = 2

    !Universal variables
    centrePx(1) = 206.0 ; centrePx(2) = 283.0   !Centrepoint of ingoing beam and liquid surface area of interaction
    startImg = 74                               !Discharge-probe delay of the first image in the sequence
    numimg = 53                                 !Number of images in the sequence
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

    do j = 1, numWedges
        wedges(j) = wedgeAngle*j            !angle of each wedge
    end do

    wedgesdeg = wedges-90                               !saves the wedges in degrees before conversion to radians for file naming purposes
    wedges = wedges-(90+(wedgeAngle/2))                 !this ensures the wedges run from -ve to +ve angles and shifts them by half-measure to centre the ROI on specific angles
    wedges = wedges*((2D0*3.141592653589793D0)/360D0)   !converts from degrees to radians
    
    allocate(ArcROI(numArcs,numWedges,startImg:startImg + ((numImg*2)-2),2))
    ArcROI = 0
    acceptedArc = 0
  
    !loop for surface in, loop for surface out
    do k = 1, 2
        do i = startImg, ((startImg + (2*numimg))-2), 2
            if (k == 1) then

                chc_index =  index(trim(sin_path), "ChC", back=.TRUE.)
                
                write(filename,sin_path(1:chc_index + 2)//'(I0.3)') i
                open(10+i,file=trim(filename))
            else

                chc_index =  index(trim(sout_path), "ChC", back=.TRUE.)

                write(filename,sout_path(1:chc_index + 2)//'(I0.3)') i
                open(10+i,file=trim(filename))
            end if
        end do

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
                                    ArcROI(acceptedArc, i, m, k) = ArcROI(acceptedArc, i, m, k) + image(x,y)
                                end if
                            end if
                        end do
                    end if
                end if
            end do
        end do
        do i = startImg, ((startImg + (2*numimg))-2), 2
            close(10+i)
        end do
    end do

     !These parameters HAVE to be specified by the user
    header_lines = 4                            !Number of header lines in the input files
    start_delay = 82                            !Delay from which to start calculating the correction
    stop_delay = 96                             !Delay at which to finish calculating the correction
    num_delays = ((96-80)/2)+1                  !Number of delays (for 2 us timesteps)
    num_entries = 0 - header_lines              !Number of timpeoints, i.e. images, counted below
    tolerance = 1E-12                           !Tolerance of the least squares fit
    upper = 1.5                                 !Upper limit of the brackets for least squares fit 
    lower = 0.5                                 !Lower limit of the brackets for least squares fit
    
    allocate(SurfaceIn(numArcs,numWedges,2,numimg))
    allocate(SurfaceOut(numArcs,numWedges,2,numimg))
    allocate(corr_data(numArcs,numWedges,2,numimg))
    corr_data = 0
    SurfaceIn = 0
    SurfaceOut = 0

    allocate(factor(numArcs,numWedges))

    do arc_index = 1, numArcs
        do wedge_index = 1, numWedges

    !Construct Surface in and surface out arrays in form of "delay, intensity" on each row
    
            do i = 1, numimg
                SurfaceIn(arc_index,wedge_index,1,i) = (i*startImg) + ((i-1)*2)
                SurfaceIn(arc_index,wedge_index,2,i) = ArcROI(arc_index,wedge_index,i,1)
            end do

            do i = 1, numimg
                SurfaceOut(arc_index,wedge_index,1,i) = (i*startImg) + ((i-1)*2)
                SurfaceOut(arc_index,wedge_index,2,i) = ArcROI(arc_index,wedge_index,i,2)
            end do               

            allocate(corr_data(numArcs,numWedges,2,numimg))
            corr_data = 0

            !Calls the subroutine that does the least squares analysis
            call least_squares_fit(SurfaceIn(arc_index,wedge_index,2,:), SurfaceOut(arc_index,wedge_index,2,:), lower, upper, tolerance, factor(arc_index,wedge_index))

            print *, factor(arc_index,wedge_index)

            corr_data(arc_index, wedge_index,1,:) = SurfaceIn(arc_index, wedge_index,1,:)
            corr_data(arc_index, wedge_index,2,:) = SurfaceIn(arc_index, wedge_index,2,:) - (factor(arc_index, wedge_index)*SurfaceOut(arc_index, wedge_index,2,:))

        end do
    end do

    do arc_index = 1, numArcs
        do wedge_index = 1, numWedges
            do i = 1, numimg
                open(unit=10,file="<path>"//"_"//trim(transition_str)//"_"//trim(surface_str)//".csv")
                write(10,'(a)') "Delay, Signal,"
                delay = nint(corr_data(arc_index, wedge_index,1,i))
                signal = corr_data(arc_index, wedge_index,2,i)
                write(10,'(I3,a,ES12.5,a)') delay, ", ", signal, ","
            end do
        end do
    end do
    
    contains

        subroutine least_squares_fit (data_array, background_array, low_limit, high_limit, tolerance, factor)
            implicit none
        
            integer :: i, array_size
            double precision, intent(in) :: low_limit, high_limit
            double precision, intent(out) :: factor, tolerance
            double precision, intent(in), dimension(:) :: data_array, background_array
            double precision, dimension(3) :: bracket
            double precision :: sum_square_lower, sum_square_upper, upper, lower
        
            array_size = SIZE(data_array)
        
            bracket(3) = high_limit
            bracket(1) = low_limit
        
            bracket(2) = bracket(3) - ((bracket(3) - bracket(1)) / 2.0)
        
            sum_square_lower = 0
            sum_square_upper = 0
        
            lower = bracket(2) - ((bracket(2) - bracket(1)) / 2.0)
            upper = bracket(3) - ((bracket(3) - bracket(2)) / 2.0)
        
            do
                sum_square_lower = 0
                sum_square_upper = 0
        
                do i = 1, array_size
                    sum_square_lower = sum_square_lower + ((data_array(i) - (lower*background_array(i)))**2.0)
                    sum_square_upper = sum_square_upper + ((data_array(i) - (upper*background_array(i)))**2.0)
                end do
        
                if (sum_square_lower .lt. sum_square_upper) then
                    if (abs(bracket(3) - bracket(1)) .lt. tolerance) then
                        factor = lower
                        EXIT
                    else
                        bracket(3) = bracket(2)
                        bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)
        
                        upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                        lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
                    end if
                else
                    if (abs(bracket(3) - bracket(1)) .lt. tolerance) then
                        factor = upper
                        EXIT
                    else
                        bracket(1) = bracket(2)
                        bracket(2) = bracket(3) - ((bracket(3) - bracket(1))/2.0)
        
                        upper = bracket(3) - ((bracket(3) - bracket(2))/2.0)
                        lower = bracket(2) - ((bracket(2) - bracket(1))/2.0)
                    end if
                end if
            end do
        end subroutine least_squares_fit

        subroutine parse_path(path, transition_str, surface_str)
            implicit none

            integer :: transition_index, surface_index
            character(500), intent(in) :: path
            character(4), intent(out) :: transition_str, surface_str
        
            transition_index = index(trim(path), "Q1", back=.TRUE.)
        
            transition_str = path(transition_index:transition_index+2)
        
            surface_index = transition_index + 4
            select case (path(surface_index:surface_index))
                case("I")
                    surface_str = "IB"
                case("S")
                    if (path(surface_index+2:surface_index+2) == "A") then
                        surface_str = "SQA"
                    else
                        surface_str = "SQE"
                    end if
                case("P")
                    surface_str = "PFPE"
            end select
        end subroutine
end program ralss