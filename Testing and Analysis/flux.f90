include "../src/pchip.f90"

program flux
    use pchip_module

    implicit none

    integer :: i, j, k, image, numimg, startImg, numAngles, sampleGap, radius, roiSize, ic(2), n, &
        incfd, nwk, ierr, ia, ib
    double precision, allocatable, dimension(:) :: angle
    integer, dimension(2) :: centrePx, roi
    double precision, allocatable, dimension(:,:) :: fluxPoints
    double precision, allocatable, dimension(:,:,:) :: imageMatrix
    double precision, allocatable, dimension(:,:) :: roiSquare
    double precision, allocatable :: x(:), f(:,:), d(:,:), wk(:,:)
    double precision :: vc(2), value
    character(100) :: filename
    logical :: skip

    ! define the centre of the surface relative to image. These numbers should match current simulation
    centrePx(1) = 210; centrePx(2) = 294

    ! define radius of region of interest in pixels (currently 1 px = 0.00025 m)
    radius = 100

    ! define number of angle with which to sample
    ! works across 180 degree range, so 36 angles corresponds to 5 degrees between data points
    ! MUST BE A FACTOR OF 180
    numAngles = 36

    ! number of images to sample across from appearance profile set
    numimg = 40

    ! image number to start with in sequence (example: start at "Image 21.txt" would mean setting this to 21)
    startImg = 70

    sampleGap = 1

    numimg = numimg / sampleGap

    ! while technically ROI size, this is really ((length of ROI box) - 1) / 2
    ! ROI size of 4 means box is 9 pixels wide
    roiSize = 4

    ! allocate arrays
    allocate(angle(numAngles))
    ! fluxPoints is the array containing the data points we will integrate with.
    ! will contain the dependant variables for each given angle
    allocate(fluxPoints(numAngles,numimg))
    ! working image storage (image files will be loaded to this array)
    allocate(imageMatrix(420,420,numimg))
    ! array containing the data for the ROI. Pixels will be extracted from images and put into this array for summation
    allocate(roiSquare(-roiSize:roiSize,-roiSize:roiSize))

    ! sets array to 0 to ensure actual values are being displayed
    fluxPoints = 0

    ! populate list of angles to loop through
    ! important: angles end up as radians for calculation later
    do i = 1, numAngles
        angle(i) = (i-1)*(180/numAngles) -90
        angle(i) = angle(i)*((2D0*3.141592653589793D0)/360D0)
    end do

    ! loads images into memory
    do i = 0, (numimg - 1)
        if ((i+startImg) .lt. 100) then
            write(filename,'("../Images2/Image",I3,".txt")') ((i*sampleGap)+startImg)
        else
            write(filename,'("../Images2/Image",I3,".txt")') ((i*sampleGap)+startImg)    
        end if
        open(10+i,file=trim(filename))
        do j = 1, 420    
            read(10+i,*) (imageMatrix(k,j,i+1),k=1,420)
        end do
    end do

    ! loops across all angles chosen
    do i = 1, numAngles
        ! finds pixel cooridinates of ROI centre
        roi(1) = floor(sin(angle(i))*radius)
        roi(2) = floor(cos(angle(i))*radius)
        !print *, i, angle(i)
        do image = 1, numimg
            do j = -roiSize, roiSize
                do k = -roiSize, roiSize
                    if (((centrePx(1)+roi(1)+j) .lt. 420) .and. ((centrePx(2)-roi(2)+k) .lt. 420)) then
                        roiSquare(j,k) = imageMatrix( centrePx(1)+roi(1)+j , centrePx(2)-roi(2)+k ,image)
                    else
                        print *, "Referenced pixel outside of image"
                    end if
                end do
            end do
            ! for each angle, sequence of data points are obtained representing ROI intensity at each time point
            fluxPoints(i,image) = SUM(roiSquare)
            if (fluxPoints(i,image) == 0) then
                !print *,  (centrePx(1)+roi(1)+j), (centrePx(2)-roi(2)+k), angle(i)
            end if
        end do
    end do

    ! interpolate the function between points
    ! Please be aware that these variables are parameters for the PCHIP package
    ! Consult pchip.f90 for full descriptions of these parameters
    ic(1) = 0 ; ic(2) = 0
    n = numimg
    incfd = 1
    nwk = 2*numimg
    ! ia and ib are the limits for the definite integral
    ! these are probably the only things you might want to change
    ia = 1
    ib = numimg
    skip = .false.

    allocate(x(numimg))
    allocate(f(incfd,numimg))
    allocate(d(incfd,numimg))
    allocate(wk(2,numimg))

    ! x is independant variable
    ! simple integration being done here assuming each x point is simply 1 unit
    do i = 1, numimg
        x(i) = i*sampleGap
    end do

    ! performs the interpolation and then integration
    ! please consult pchips.f90 for description of parameters, it's a bit messy
    do j =1, numAngles
        do i = 1, numimg
            f(1,i) = fluxPoints(j,i)
        end do

        ! interpolation: finds function with which we can perform integration
        call dpchsp (ic, vc, n, x, f, d, incfd, wk, nwk, ierr)

        ! integration, value = definite integral
        value = dpchid (n, x, f, d, incfd, skip, ia, ib, ierr)

        ! returns integrated singal for each angle
        print *, value
    end do

end program flux