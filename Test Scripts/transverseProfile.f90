program transverseProfile
    implicit none

    integer :: row, column, xPx, yPx
    integer, dimension (420,420) :: intensityCount

    xPx = 420
    yPx = 420

    open(unit=11,file='Images/Image 9.txt')
    open(unit=12,file='transverseProfile.txt')

    do row = 1, xPx

        read(11,*) (intensityCount(row,column),column=1,yPx)

    end do

    do column = 1, yPx

    
        write(12,*) column, SUM(intensityCount(220:250,column))


    end do

end program transverseProfile