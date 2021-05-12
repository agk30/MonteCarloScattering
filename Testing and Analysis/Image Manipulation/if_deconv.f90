program if_deconv
    implicit none

    double precision, dimension(420,420) :: image, ifImage
    double precision :: background
    integer :: i, j

    open(10,file="../../Real Images/16022021_Q12_Squalane_SUM_094")
    open(11,file="../smoothed_IF.txt")
    open(12,file="deconv_image.txt")

    do i = 1, 420    
        read(10,*) (image(j,i),j=1,420)
    end do

    background = sum(image(1:75,1:75)) / 5625D0

    image = image - background

    do i = 1, 420    
        read(11,*) (ifImage(j,i),j=1,420)
    end do

    ifImage = ifImage/maxval(ifImage)

    do i = 1, 420
        do j = 1, 420
            if (ifImage(i,j) .gt. 0.4) then
                image(i,j) = image(i,j) / ifImage(i,j)
            else
                image(i,j) = 0
            end if
        end do
    end do

    image(206,280) = 16000

    do i = 1, 420
        do j = 1, 420
            write(12,'(ES12.5,a)',advance='no') image(i,j)," "
        end do
        write(12,*)
    end do

end program if_deconv