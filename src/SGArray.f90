!TODO fix this damn thing. Bad hardcoded variables - get it done pls

module sgconv
    use mathConstants

    contains

    ! builds the kernel used in 3D SG smoothing subroutine 
    ! returns 1D array in wraparound order from the 2D padded input array   
    subroutine constructKernel(i, j, ksize, nd, columnKernel, padInput)
            implicit none
            
            integer :: k, l
            integer, intent(in) :: ksize, nd, i, j
            double precision, intent(inout), dimension(:,:) :: padInput
            double precision, intent(inout), dimension(:) :: columnKernel

            
            do k = 1, ksize
                columnKernel(((k-1)*ksize)+1 : k*ksize) = padInput((i+nd+nd-k+1), j:(j+ksize-1))
            end do

        end subroutine constructKernel

        subroutine sgarray(xPx, zPx, ksize, polyOrder, output)
        
            implicit none
        
            integer :: xPx, zPx, row, column, i, j, k, l, m, ksize, nd, nk, polyOrder
            double precision, allocatable, dimension(:,:) :: input, diffinput
            double precision, allocatable, dimension(:,:), intent (out) :: output
            double precision, allocatable, dimension(:,:) :: padInput
            double precision, allocatable, dimension(:) :: sgmatrix
            double precision, allocatable, dimension(:,:) :: kernel
            double precision, allocatable, dimension(:) :: columnKernel
            double precision :: dotprod
            
            nd = (ksize-1)/2
            nk = ((ksize**2) - 1)/2
        
            allocate(input(xPx,zPx))
            allocate(output(xPx,zPx))
            allocate(diffinput(xPx,zPx))
            allocate(padInput(xPx+ksize-1,zPx+ksize-1))
            allocate(sgmatrix(ksize**2))
            allocate(kernel(-nd:nd,-nd:nd))
            allocate(columnKernel(ksize**2))

            ! TODO ask for input of real IF file on start up, fix concatenation
           ! open(11,file='../Real Images/02032020_1_Q11_IF')
            open(11,file='../Real Images/01102020_4_Q11_IF')
           ! open(12,file='SG Matrices/CC_027x027_00'//char(polyOrder)//'x00'//char(polyOrder)//'.dat')
            open(12,file='../SG Matrices/CC_027x027_003x003.dat')
        
            read(12,*) sgmatrix(1:(((ksize**2)-1)/2)+1)

            print *, "1"
        
            ! input SG matrices only contain half + 1 of the required array, symmetric about the last entry. this loop builds
            ! the rest of the array
            do i = 1, nk     
                sgmatrix(i + nk + 1) = sgmatrix(nk+1-i)    
            end do
        
            ! reads IF image into array
            do row = 1, xPx       
                read(11,*) (input(row,column),column=1,zPx)       
            end do

            ! TODO make this a variable. Sets the baseline of the image to 0, or close to 0
            input = input - 950
            padInput = 0D0

            print *, "2"
        
            ! builds a padded array. the image array requires padding outside of its bounds with the addition of nd fields on 
            ! either side of the image and also above and below it
            do i = 1, 420     
                do j = 1, 420
        
                    padInput(i+nd,j+nd) = input(i,j)
        
                end do       
            end do

            print *, "3"
      
            do i = 1, 420        
                do j = 1, 420      
                    call constructKernel(i, j, ksize, nd, columnKernel, padInput)
       
                    ! convolutes image with SG cooefficients
                    dotprod = dot_product(columnKernel, sgmatrix)
        
                    output(i,j) = dotprod

                    !ensures no negative numbers in output (not really needed, should be fixed)
                    if (output(i,j) .lt. 0) then
                        !output(i,j) = 0
                    end if  
                end do  
            end do

            print *, "4"

            ! extra info. shows differences in input image and the output image.
            !do i = 1, 420   
                !do j = 1, 420    
                 ! diffinput(i,j) = input(i,j) - output(i,j)    
                !end do  
            !end do
        end subroutine sgarray

end module sgconv



