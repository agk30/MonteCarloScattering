include "../src/hwlib/SGArray.f90"

program sg_smooth
    use sgconv
    implicit none

    integer :: row, column, i ,j, image
    double precision, dimension(420,420) :: input, ifImage, emptyImage, output
    character(200) :: matrixPath

    input = 0
    ifImage = 0

    matrixPath = "../SG Matrices/CC_027x027_003x003.dat"

    open(10,file="../Real Images/16022021_Q12_Squalane_SUM/16022021_Q12_Squalane_SUM_068")
    open(500,file="deconv_image.txt")
    open(11,file="../Real Images/09032021_Q12_Instrument Function_SUM_00")
    open(12,file="../Real Images/16022021_Q12_Squalane_SUM/16022021_Q12_Squalane_SUM_094")
    open(600,file="smoothed_IF.txt")

    do row = 1, 420      
        read(10,*) (emptyImage(row,column),column=1,420)       
    end do

    !do image = 68, 178, 2

        ! reads IF image into array
        do row = 1, 420      
            read(11,*) (input(row,column),column=1,420)       
        end do

        input = input - emptyImage

        call sg_array (420, 420, 27, matrixPath, input, ifImage)

        ifImage = ifImage / maxval(ifImage)

        do i = 1, 420
            do j = 1, 420
                write(600,'(ES12.5,a)',advance='no') ifImage(i,j)," "
            end do

            write(600,*)
        end do

        do row = 1, 420      
            read(12,*) (output(row,column),column=1,420)       
        end do

        output = output - emptyImage

        do i = 1, 420
            do j = 1, 420
                if (ifImage(i,j) .gt. 0.3) then
                    output(i,j) = output(i,j) / ifImage(i,j)
                else
                    output(i,j) = 0
                end if
            end do
        end do
    
        do i = 1, 420
            do j = 1, 420
                write(500,'(ES12.5,a)',advance='no') output(i,j)," "
            end do
            write(500,*)
        end do
    
    !end do

end program sg_smooth