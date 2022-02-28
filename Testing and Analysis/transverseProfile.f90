program transverseProfile
    implicit none

    integer :: row, column, xPx, yPx, pxSelection, currPx, i
    real, dimension (420,420) :: intensityCount, background
    character*100 :: filename, ifname, rowString

    xPx = 420
    yPx = 420

    open(unit=11,file= "C:\Users\adam\Documents\Data\05072021_1_Q11_IB_TOF Profile\05072021_1_Q11_IB_TOF Profile_ChC134")
    !open(unit=11,file= '06112019_1_Q11_IB TOF Profile_ChC106')

    open(unit=20,file= "C:\Users\adam\Documents\Data\05072021_1_Q11_IB_TOF Profile\05072021_1_Q11_IB_TOF Profile_ChC098")

    open(unit=12,file= 'transverseProfile.txt')

    do row = 1, xPx
        read(11,*) (intensityCount(row,column),column=1,yPx)
    end do

    do row = 1, xPx
        read(20,*) (background(row,column),column=1,yPx)
    end do

    intensityCount = intensityCount - background

    do column = 1, 420
        write((12),*) sum(intensityCount(175:185,column))
    end do

end program transverseProfile