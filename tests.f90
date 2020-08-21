module tests
    use mathConstants

    integer, dimension(-90:90) :: angleDist

    contains

     subroutine angleDistribution(position)
        implicit none

        real(kind=r14), dimension(3), intent(in) :: position
        integer :: angle

        angle = floor(atand(position(1) / position(3)))

        if ((angle .gt. -90) .and. (angle .lt. 90)) then

            angleDist(angle) = angleDist(angle) + 1

        end if

     end subroutine angleDistribution

     subroutine writeAngleDistribution
        implicit none

        integer :: i

        open(unit=500,file='angledist.txt')

        do i = -90, 90

            write(500,*) angleDist(i), i

        end do

     end subroutine writeAngleDistribution

end module tests
