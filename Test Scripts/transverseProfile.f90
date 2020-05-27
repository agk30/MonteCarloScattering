program transverseProfile
    implicit none

    integer :: row, column, xPx, yPx
    real, dimension (420,420) :: intensityCount

    xPx = 420
    yPx = 420

    open(unit=11,file='Real Images/01102019_1_Q11_IB TOF Profile_ChC102')
    open(unit=12,file='realTransverseProfile.txt')

    do row = 1, xPx

        read(11,*) (intensityCount(row,column),column=1,yPx)

    end do

    do column = 1, yPx

    
        write(12,*) column, SUM(intensityCount(220:250,column))


    end do

end program transverseProfile