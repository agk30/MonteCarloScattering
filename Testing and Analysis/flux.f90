include "../src/pchip.f90"

program flux
    use pchip_module

    implicit none

    integer :: i, j, k, image, numimg, startImg, numAngles, radius, roiSize, ic(2), n, &
        incfd, nwk, ierr, ia, ib
    double precision, dimension(36) :: angle
    integer, dimension(2) :: centrePx, roi
    double precision, dimension(36,40) :: fluxTotal
    double precision, dimension(420,420,40) :: imageMatrix
    double precision, dimension(-4:4,-4:4) :: roiSquare
    double precision :: vc(2), x(40), f(1,40), d(1,40), wk(2,40), value
    character(100) :: filename
    logical :: skip

    fluxTotal = 0

    centrePx(1) = 210; centrePx(2) = 294

    radius = 90

    do i = 1, 36
        angle(i) = (i-1)*5 -90
        angle(i) = angle(i)*((2D0*3.141592653589793D0)/360D0)
        !roi (2) = floor(cos(angle(i))*radius)
        !if (roi(2) .gt. 24) then
        !    print *, i
        !end if
    end do

    !angle = angle*((2D0*3.141592653589793D0)/360D0)

    numAngles = 36

    numimg = 40

    startImg = 70

    roiSize = 4

    ! loads images into memory
    do i = 1, numimg
        if ((i+startImg) .lt. 100) then
            write(filename,'("../Images2/Image",I3,".txt")') (i+startImg)
        else
            write(filename,'("../Images2/Image",I3,".txt")') (i+startImg)    
        end if
        open(10+i,file=trim(filename))
        do j = 1, 420    
            read(10+i,*) (imageMatrix(k,j,i),k=1,420)
        end do

    end do

    ! loops across all angles chosen
    do i = 1, numAngles
        roi(1) = floor(sin(angle(i))*radius)
        roi(2) = floor(cos(angle(i))*radius)

        !print *, roi(1), roi(2), angle(i)
        
        do image = 1, numimg
            do j = -4, 4
                do k = -4, 4
                    roiSquare(j,k) = imageMatrix( centrePx(1)+roi(1)+j , centrePx(2)-roi(2)+k ,image)
                end do
            end do
            ! for each angle, sequence of data points are obtained representing ROI intensity at each time point
            fluxTotal(i,image) = SUM(roiSquare)
        end do
    end do

    ! interpolate the function between points
    ic(1) = 0 ; ic(2) = 0

    n = 40

    incfd = 1

    nwk = 80

    ia = 1

    ib = 40

    skip = .false.

    ! x is independant variable
    do i = 1, 40
        x(i) = i
    end do

    ! performs the interpolation and then integration
    ! please consult pchips.f90 for description of parameters, it's a bit messy
    do j =1, 36
        do i = 1, 40
            f(1,i) = fluxTotal(j,i)
            !print *, f(1,i)
        end do

        ! interpolation: finds function with which we can perform integration
        call dpchsp (ic, vc, n, x, f, d, incfd, wk, nwk, ierr)

        !print *, ierr

        ! integration, value = definite integral
        value = dpchid (n, x, f, d, incfd, skip, ia, ib, ierr)

        print *, value

    end do

end program flux