program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image
    double precision, dimension(-4:4,-4:4,5,5,40) :: roi
    double precision, dimension(2) :: centrePx, roiCentre
    double precision, dimension(5) :: roiRadius
    double precision, dimension(5) :: angle
    integer :: i, j, k, l, m, startImg, numimg, roiSize
    character(100) :: filename

    centrePx(1) = 210.0 ; centrePx(2) = 323.0
    roiRadius(1) = 50.0 ; roiRadius(2) = 60.0 ; roiRadius(3) = 70.0 ; roiRadius(4) = 80.0 ; roiRadius(5) = 90.0
    angle(1) = -45D0 ; angle(2) = -22.5 ; angle(3) = 0 ; angle(4) = 22.5D0 ; angle(5) = 45

    roiCentre(1) = centrePx(1)

    angle = angle*((2D0*3.141592653589793D0)/360D0)

    startImg = 29

    numimg = 40

    roiSize = 0
 
    do i = 1, numimg
        if ((i+startImg) .lt. 100) then
            write(filename,'("Images/Image",I3,".txt")') (i+startImg)
        else
            write(filename,'("Images/Image",I3,".txt")') (i+startImg)    
        end if
        open(10+i,file=trim(filename))
    end do
    


    open(100,file='roidata.txt')


    ! loops through images
    do m = 1, numimg

        do i = 1, 420
            
            read(10+m,*) (image(j,i),j=1,420)

        end do

        ! loop through angles
        do l = 1, 5

            ! loop through radii
            do k = 1, 5

                roiCentre(1) = floor(sin(angle(l))*roiRadius(k))
                roiCentre(2) = floor(cos(angle(l))*roiRadius(k))

                do i = -roiSize, roiSize   
                    do j = -roiSize, roiSize
                        roi(i,j,k,l,m) = image((int(centrePx(1)) + int(roiCentre(1)) + i),&
                         (int(centrePx(2)) - int(roiCentre(2) + j)))
                        end if
                    end do
                end do

            end do

        end do

    end do
    do k = 1, 5

        do j = 1, 5

            write(100,*)  "Angle", angle(k), "ROI", j

            do i = 1, numimg

                ! j for radius, k for angle
                write(100,*) SUM(roi(:,:,j,k,i))

            end do

        end do 
    
    end do

end program roiAnalysis