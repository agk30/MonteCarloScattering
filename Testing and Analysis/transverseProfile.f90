program transverseProfile
    implicit none

    integer :: row, column, xPx, yPx, pxSelection, currPx, i
    real, dimension (420,420) :: intensityCount
    character*100 :: filename, ifname, rowString

    xPx = 420
    yPx = 420

    open(unit=11,file= '../Images2/Image 44.txt')
    !open(unit=11,file= '06112019_1_Q11_IB TOF Profile_ChC106')

    open(unit=12,file= 'transverseProfile.txt')

    do row = 1, xPx
        read(11,*) (intensityCount(row,column),column=1,yPx)
    end do

    do column = 1, 420
        write((12),*) sum(intensityCount(205:215,column))
    end do

end program transverseProfile