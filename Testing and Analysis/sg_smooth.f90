include "../src/hwlib/SGArray.f90"

program sg_smooth
    use sgconv
    implicit none

    integer :: row, column
    double precision, dimension(420,420) :: input, output
    character(100) :: matrixPath

    matrixPath = "../SG Matrices/CC_027x027_003x003.dat"

    open(11,file="../Real Images/")

    ! reads IF image into array
    do row = 1, 420      
        read(11,*) (input(row,column),column=1,420)       
    end do

    call sg_array (420, 420, 27, trim(matrixPath), input, output)

end program sg_smooth