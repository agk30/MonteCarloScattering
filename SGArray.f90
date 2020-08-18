!TODO fix this damn thing. Bad hardcoded variables - get it done pls

module sgconv

    contains

        subroutine constructKernel(i, j, ksize, nd, columnKernel, padInput)
            implicit none
            
            integer :: k, l
            integer, intent(in) :: ksize, nd, i, j
            double precision, intent(inout), dimension(:,:) :: padInput
            double precision, intent(inout), dimension(:) :: columnKernel

            
            do k = 1, ksize

                !do l = 1, ksize
                
                !print *, k


                columnKernel(((k-1)*ksize)+1 : k*ksize) = padInput(i+nd-k+nd, j:j+ksize)

                !end do

            end do

        end subroutine constructKernel

        subroutine sgarray(output)
        
            implicit none
        
            integer :: xPx, yPx, row, column, i, j, k, l, m, ksize, nd, nk
            double precision, allocatable, dimension(:,:) :: input, diffinput
            double precision, allocatable, dimension(:,:), intent (out) :: output
            double precision, allocatable, dimension(:,:) :: padInput
            double precision, allocatable, dimension(:) :: sgmatrix
            double precision, allocatable, dimension(:,:) :: kernel
            double precision, allocatable, dimension(:) :: columnKernel
            double precision :: dotprod
            character*100 :: filename, ifname
            
        
            xPx = 420
            yPx = 420
            ksize = 27
            nd = (ksize-1)/2
            nk = ((ksize**2) - 1)/2
        
        
            allocate(input(xPx,yPx))
            allocate(output(xPx,yPx))
            allocate(diffinput(xPx,yPx))
            allocate(padInput(xPx+ksize-1,yPx+ksize-1))
            allocate(sgmatrix(ksize**2))
            allocate(kernel(-nd:nd,-nd:nd))
            allocate(columnKernel(ksize**2))
        
            
            open(11,file='Real Images/02032020_1_Q11_IF')
            open(12,file='SG Matrices/CC_027x027_003x003.dat')
        
        
            read(12,*) sgmatrix(1:(((ksize**2)-1)/2)+1)
        
            
        
            do i = 1, nk
        
                sgmatrix(i + nk + 1) = sgmatrix(nk+1-i)
                !sgmatrix((ksize**2 - 1)/2 + 2: ksize**2) = sgmatrix((ksize**2 - 1)/2 : 1)
        
            end do
        
            do row = 1, xPx
        
                read(11,*) (input(row,column),column=1,yPx)
        
            end do

            input = input - 950
            padInput = 0D0
        
            do i = 1, 420
        
                do j = 1, 420
        
                    padInput(i+nd,j+nd) = input(i,j)
        
                end do
        
            end do
        
            do i = 1, 420
        
                do j = 1, 420
        
                    call constructKernel(i, j, ksize, nd, columnKernel, padInput)
        
                    dotprod = dot_product(columnKernel, sgmatrix)
        
                    output(i,j) = dotprod

                    if (output(i,j) .lt. 0) then

                        output(i,j) = 0

                    end if
        
                end do
        
            end do
        
            do i = 1, 420
        
                do j = 1, 420
        
                    diffinput(i,j) = input(i,j) - output(i,j)
        
                end do
        
            end do
        end subroutine sgarray

end module sgconv



