module mod_tests
    use mathConstants

    integer, dimension(-90:90) :: angleDist

    contains

        ! produces a file of binned angle distributions. 
        subroutine angleDistribution(outgoing)
            implicit none

            double precision, dimension(3), intent(in) :: outgoing
            double precision, dimension(3) :: ingoing
            integer :: i
            double precision :: angle

            ingoing(1) = 0 ; ingoing(2) = 0 ; ingoing(3) = 1

            angle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360/(2*pi))
            
            if (outgoing(1) .lt. 0) then

                angle = -angle
                
            end if

            do i = -90, 90

                if ((angle .gt. i) .and. (angle .lt. i+1)) then

                    angleDist(i) = angleDist(i) + 1

                    exit

                end if

            end do

        end subroutine angleDistribution

        subroutine writeAngleDistribution
            implicit none

            integer :: i

            open(unit=500,file='angledist.txt')

            do i = -90, 90

                write(500,*) angleDist(i), (i)

            end do

        end subroutine writeAngleDistribution

        subroutine angle_speed_distribution(vector, speed, distributionBin)
            implicit none

            double precision, intent(in) :: vector(3), speed
            double precision :: normal(3), normalxz(2), normalyz(2), vectorxz(2), vectoryz(2)
            double precision :: anglexz, angleyz, maxSpeed
            integer :: distributionBin(:,:)
            integer :: i, j, angleBinSize, speedBinSize, intSpeed

            angleBinSize = 2

            speedBinSize = 10

            maxSpeed = 2500

            normal(1) = 0 ; normal(2) = 0 ; normal(3) = 1.0
            normalxz(1) = 0 ; normalxz(2) = 1
            normalyz(1) = 0 ; normalyz(2) = 1

            vectorxz(1) = vector(1) ; vectorxz(2) = vector(3)
            vectoryz(1) = vector(2) ; vectoryz(2) = vector(3)

            if ((speed .lt. maxSpeed) .and. (vector(1) .gt. 0)) then

                anglexz = acos((dot_product(vectorxz,normalxz)/(norm2(vectorxz)*norm2(normalxz))))
                angleyz = acos((dot_product(vectoryz,normalyz)/(norm2(vectoryz)*norm2(normalyz))))
                anglexz = ceiling((anglexz*(360.0/(pi*2)))/angleBinSize)
                angleyz = ceiling((angleyz*(360.0/(pi*2)))/angleBinSize)

                

                if ((anglexz .ge. 1) .and. (anglexz .le. 45)) then

                    intSpeed = ceiling(speed/speedBinSize)

                    !print *, speed

                    distributionBin(int(anglexz),intSpeed) = distributionBin(int(anglexz),intSpeed) + 1
                end if

            end if

        end subroutine angle_speed_distribution

        subroutine write_angle_speed(bin)

            integer, dimension(:,:) :: bin
            integer :: i, j, angleBinSize, speedBinSize, maxSpeed
            character :: c

            angleBinSize = 2

            speedBinSize = 10

            maxSpeed = 2500

            open(600,file='angle_speed_distribution.csv')

            write(600,'(a)',advance='no') "Speed m/s,"

           ! do i = 1, (90/angleBinSize)
              !  write(600,'(I2)',advance='no') (i*angleBinSize)

              !  if (i .lt. (90/angleBinSize)) then
               !     write(600,'(a)',advance='no') ","
               ! else
                  !  write(600,'(a)',advance='no') new_line(c)
                !end if
            !end do
            do j = 1, (90/angleBinSize)
                do i = 1, (maxSpeed/speedBinSize)
               
                    write(600,'(I4,a)',advance='no') i*speedBinSize,","
                    write(600,'(I4,a)',advance='no') j*angleBinSize,","
                    write(600,'(I5,a)',advance='no') bin(j,i),","
                    write(600,'(a)',advance='no') new_line(c)
                    !if (j .gt. (90/angleBinSize)) then
                    !    write(600,'(a)',advance='no') new_line(c)
                    !end if
                end do
            end do

        end subroutine write_angle_speed

end module mod_tests
