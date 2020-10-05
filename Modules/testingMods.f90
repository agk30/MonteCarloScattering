module mod_tests
    use mathConstants

    integer, dimension(-90:90) :: angleDist

    contains

        ! produces a file of binned angle distributions. 
        subroutine angleDistribution(outgoing)
            implicit none

            real(kind=r14), dimension(3), intent(in) :: outgoing
            real(kind=r14), dimension(3) :: ingoing
            integer :: i
            real(kind=r14) :: angle

            ingoing(1) = 0 ; ingoing(2) = 0 ; ingoing(3) = 1

            angle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360/2*pi)
            
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

end module mod_tests
