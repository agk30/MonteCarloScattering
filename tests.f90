module tests
    use mathConstants

    integer, dimension(-9:9) :: angleDist

    contains

    ! produces a file of binned angle distributions. 
    subroutine angleDistribution(position)
        implicit none

        real(kind=r14), dimension(3), intent(in) :: position
        integer :: i
        real(kind=r14) :: angle

        angle = atand(position(1) / position(3)) / 10

        do i = -9, 9

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

        do i = -9, 9

            write(500,*) angleDist(i), i

        end do

     end subroutine writeAngleDistribution

end module tests
