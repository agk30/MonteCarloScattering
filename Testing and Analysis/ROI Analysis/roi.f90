program roiAnalysis
    implicit none

    double precision, dimension(420,420) :: image
    double precision, dimension(-4:4,-4:4,5,5,12) :: roi
    double precision, dimension(2) :: centrePx, roiCentre
    double precision, dimension(5) :: roiRadius
    double precision, dimension(5) :: angle
    integer :: i, j, k, l, m

    centrePx(1) = 210.0 ; centrePx(2) = 297.0
    roiRadius(1) = 30.0 ; roiRadius(2) = 40.0 ; roiRadius(3) = 54.0 ; roiRadius(4) = 68.0 ; roiRadius(5) = 85.0
    angle(1) = -45 ; angle(2) = -22.5 ; angle(3) = 0 ; angle(4) = 22.5 ; angle(5) = 45

    roiCentre(1) = centrePx(1)
 
    open(11,file='Images/Image 55.txt')
    open(12,file='Images/Image 60.txt')
    open(13,file='Images/Image 65.txt')
    open(14,file='Images/Image 70.txt')
    open(15,file='Images/Image 75.txt')
    open(16,file='Images/Image 80.txt')
    open(17,file='Images/Image 85.txt')
    open(18,file='Images/Image 90.txt')
    open(19,file='Images/Image 95.txt')
    open(20,file='Images/Image100.txt')
    open(21,file='Images/Image105.txt')
    open(22,file='Images/Image110.txt')

    ! loops through images
    do m = 1, 12

        do i = 1, 420
            
            read(10+m,*) (image(j,i),j=1,420)

        end do

        !print *, image(210,212)

        ! loop through angles
        do l = 1, 5

            ! loop through radii
            do k = 1, 5

                roiCentre(1) = floor(sin(angle(l))*roiRadius(k))
                roiCentre(2) = floor(cos(angle(l))*roiRadius(k))

                do i = -4, 4
                    do j = -4, 4
                        roi(i,j,k,l,m) = image((int(centrePx(1)) + int(roiCentre(1)) + i),&
                         (int(centrePx(2)) - int(roiCentre(2) + j)))
                    end do
                end do

                if ((m == 1) .and. (l == 2) .and. (k == 1)) then

                    !print *, roiCentre(1), roiCentre(2), centrePx(1), centrePx(2)
                    !print *, ((int(centrePx(1)) - int(roiCentre(1))))
                    !print *, (int(centrePx(2)) - int(roiCentre(2)))

                end if

            end do

        end do

    end do

    do k = 1, 5

        do j = 1, 5

            print *,  "Angle", angle(k), "ROI", j

            do i = 1, 12

                ! j for radius, k for angle
                print *, SUM(roi(:,:,j,k,i))

            end do

        end do 
    
    end do

end program roiAnalysis