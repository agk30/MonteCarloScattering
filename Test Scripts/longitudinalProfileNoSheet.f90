program longitudinalProfile
    implicit none

    integer :: row, column, xPx, yPx
    integer, dimension (420,420) :: intensityCount

    xPx = 420
    yPx = 420

    open(unit=11,file='No Sheet Images/Image 9.txt')
    open(unit=12,file='longitudinalProfileNoSheet.txt')

    do row = 1, yPx

        intensityCount = 0

        read(11,*) (intensityCount(row,column),column=1,xPx)

        write(12,*) row, SUM(intensityCount)

    end do

end program longitudinalProfile