program fftransverse
    implicit none

    integer :: i, j, start(2)
    double precision :: image(420,420), ifimage(420,420), profile(250)

    !open(11,file='18112020_1_Q11_IB TOF Profile_ChC092')
    open(11,file='Image  9.txt')
    open(12,file='smoothedIF.txt')

    start(1) = 0
    start(2) = 332

    profile = 0

    !print *, profile(213)

    do i = 1, 420
        read(11,*) (image(j,i),j=1,420)
    end do

    do i = 1, 420
        !read(12,*) (ifimage(j,i),j=1,420)
    end do

    !image = image / ifimage

    print *, profile(213)

    do i = 1, 250
        do j = 1, 5
            profile(i) = profile(i) + image(start(1)+i+j,start(2)-i+j)

            !print *, profile(213)
            
            !print *, profile(i), start(1)+i+j,start(2)-i+j
        end do
    end do

    do i = 1, 250
        print *, profile(i)
    end do

end program fftransverse