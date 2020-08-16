program transverseProfile
    implicit none

    integer :: row, column, xPx, yPx, pxSelection, currPx, i
    real, dimension (420,420) :: intensityCount
    character*100 :: filename, ifname, rowString

    xPx = 420
    yPx = 420

    filename = "CC_007x007_001x001"

    print *, "Enter name of IF file"

    read(*,*) ifname

    open(unit=11,file= 'Smoothed Images/'//trim(ifname))

    do row = 1, xPx

        read(11,*) (intensityCount(row,column),column=1,yPx)

    end do
    
    do i = 1, 9

        currPx = (80 + ((i-1)*30))

        write(rowString,'(I3)') currPx
    
        open(unit=(12+i),file= 'Slices/'//trim(rowString)//'_'//'longi_slice_'//trim(ifname)//'.txt')


        do row = 1, yPx

        
            write((12+i),*) intensityCount(row,currPx)


        end do

    end do

end program transverseProfile