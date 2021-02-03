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

            double precision :: vector(3), speed
            double precision :: normal(3), normalxz(2), normalyz(2), vectorxz(2), vectoryz(2)
            double precision :: anglexz, angleyz, maxSpeed
            integer :: distributionBin(:,:)
            integer :: i, j, angleBinSize, speedBinSize

            angleBinSize = 2

            speedBinSize = 10

            maxSpeed = 2500

            normal(1) = 0 ; normal(2) = 0 ; normal(3) = 1.0
            normalxz(1) = 0 ; normalxz(2) = 1
            normalyz(1) = 0 ; normalyz(2) = 1

            vectorxz(1) = vector(1) ; vectorxz(2) = vector(3)
            vectoryz(1) = vector(2) ; vectoryz(2) = vector(3)

            if ((speed .lt. maxSpeed) .and. (vector(2) .gt. 0)) then

                anglexz = cos((dot_product(vectorxz,normalxz)/(norm2(vectorxz)*norm2(normalxz))))
                angleyz = cos((dot_product(vectoryz,normalyz)/(norm2(vectoryz)*norm2(normalyz))))
                anglexz = floor((anglexz*(360.0/(pi*2)))/angleBinSize) + 1
                angleyz = floor((angleyz*(360.0/(pi*2)))/angleBinSize) + 1

                if ((anglexz .ge. 1) .and. (anglexz .le. 45)) then

                    speed = floor(speed/speedBinSize) + 1

                    !print *, speed, anglexz

                    distributionBin(int(anglexz),int(speed)) = distributionBin(int(anglexz),int(speed)) + 1
                end if

            end if

        end subroutine angle_speed_distribution

        subroutine write_angle_speed(bin)

            integer, dimension(:,:) :: bin

            !print *, bin

        end subroutine write_angle_speed

end module mod_tests
