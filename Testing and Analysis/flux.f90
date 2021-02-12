program flux
    implicit none

    integer :: i, j, k, image, numimg, startImg, numAngles, radius, roiSize
    double precision, dimension(36) :: angle
    integer, dimension(2) :: centrePx, roi
    double precision, dimension(36) :: fluxTotal
    double precision, dimension(420,420,40) :: imageMatrix
    double precision, dimension(-4:4,-4:4) :: roiSquare
    character(100) :: filename

    fluxTotal = 0

    centrePx(1) = 210; centrePx(2) = 294

    radius = 50

    do i = 1, 36
        angle(i) = (i-1)*5 -90
        angle = angle*((2D0*3.141592653589793D0)/360D0)
        !roi (2) = floor(cos(angle(i))*radius)
        !if (roi(2) .gt. 24) then
        !    print *, i
        !end if
    end do

    angle = angle*((2D0*3.141592653589793D0)/360D0)

    numAngles = 36

    numimg = 40

    startImg = 70

    roiSize = 4

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

    do i = 1, numAngles
        roi(1) = floor(sin(angle(i))*radius)
        roi(2) = floor(cos(angle(i))*radius)

        !print *, i
        
        do image = 1, numimg
            do j = -4, 4
                do k = -4, 4
                    roiSquare(j,k) = imageMatrix( centrePx(1)+roi(1)+j , centrePx(2)-roi(2)+k ,image)
                end do
            end do
            fluxTotal(i) = fluxTotal(i) + SUM(roiSquare)
        end do
    end do

    do i = 1, 36
        print *, fluxTotal(i)
    end do

end program flux