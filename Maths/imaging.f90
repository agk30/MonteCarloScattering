module imaging    
    use mathConstants
    use tests

    ! Array shared by entire class - REMEMBER ALWAYS TO ALLOCATE BEFORE USE
    ! image2 is used in the writing of images viewed from along z axis (see writeImageTest)
    integer*8, dimension(:,:,:), allocatable :: image2!, image
    logical :: anglestats

    contains

            ! Uses the entry time and exit time to find the corresponding timepoint for imaging
        subroutine startEndTimePoints(NumberOfTimepoints, entryTime, exitTime, probeStart, probeEnd, &
             tStep, startTimePoint, endTimePoint)
            implicit none

            integer, intent(in) :: NumberOfTimePoints
            real(kind=r14), intent(in) :: entryTime, exitTime, probeStart, probeEnd, tStep
            integer, intent(out) :: startTimePoint, endTimePoint
            
            ! If entry time is less than probe start time, then imaging for that particle starts from the beginning of the process
            if (entryTime .lt. probeStart) then

                startTimePoint = 1

            else

                ! Finds nearest timepoint above entry time, but + 1 is added due to the fact that logically,
                ! the image sequence starts from image 1 rather than image 0
                startTimePoint = ceiling((entryTime - probeStart) / tStep) + 1

            end if

            ! If exit time is greater than end of probe, then end timepoint then particle is imaged 
            ! all the way til the end of the probe time
            if (exitTime .gt. probeEnd) then

                endTimePoint = NumberOfTimepoints

            else

                ! Similarly to entry time, exit timepoint is found as the nearest timepoint above exit time + 1 for the same reasons
                endTimePoint = floor((exitTime - probeStart) / tStep) + 1

            end if

        end subroutine startEndTimePoints

        ! Finds the position a particle is in at any given timepoint, and finds its corresponding pixel position then writes it
        ! to the image array, adding intensity to that pixel region
        subroutine getPosInProbe(image, NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, t0, probeStart, tStep, &
             particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)
            implicit none

            real(kind=r14), intent(inout), dimension(:,:,:) :: image
            integer, intent(in) :: NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx
            real(kind=r14), intent(in) :: probeStart, tStep, particleSpeed, pxMmRatio, t0
            real(kind=r14), dimension(3), intent(in) :: particleVector, particleStartPos, sheetDimensions
            integer :: t, posInProbexPx, posinProbeyPx, posInProbezPx, sheetCentrePx, yPx, i
            real(kind=r14) :: currentTime, angle
            real(kind=r14), dimension(3) :: posInProbe
            logical :: zImage
            
            ! for testing angle distribution. To be moved to separate module.
            anglestats = .true.

            ! testing purpose: should be left as false for normal operation, other inputs would be required in normal 
            ! operation to enable this properly.
            zImage = .false.

            ! for sake of creating image viewed along z axis, y pixels are the same as z but can be changed if needed
            yPx = zPx

            ! Loops from entry timepont to exit timepoint to avoic wasting cycles when particle is not within sheet
            do t = startTimePoint, endTimePoint

                ! currentTime refers to the time it has taken the particle to travel from its starting point
                ! to the point in space at the given timepoint
                currentTime = probeStart + (t-1)*tStep - t0

                ! Real space position for particle
                posInProbe(1) = particleStartPos(1) + (particleVector(1)*particleSpeed*currentTime)
                posInProbe(2) = particleStartPos(2) + (particleVector(2)*particleSpeed*currentTime)
                posInProbe(3) = particleStartPos(3) + (particleVector(3)*particleSpeed*currentTime)

                ! Relative pixel position for particle
                ! Note: the subtraction at the end of each statement alters the position of the particle within the image array.
                ! Altering the x-postion by half the width of the image centres the beam
                ! Similarly, the z position can be altered however a factor of 1.3 was found to centre the sheet within the middle
                ! of the image quite well
                posInProbexPx = (ceiling(posInProbe(1)/pxMmRatio) + floor(real(xPx/2)))
                posInProbeyPx = (ceiling(posInProbe(2)/pxMmRatio) + floor(real(yPx/2)))
                posInProbezPx = abs(ceiling(posInProbe(3)/pxMmRatio) - floor(real(zPx/1.3)))

                if ((anglestats .eqv. .true.) .and. (t == 83)) then

                    if (((SQRT(posInProbe(1)**2 + posInProbe(3)**2)) .lt. 0.015) .and. ((SQRT(posInProbe(1)**2 + posInProbe(3)**2)) &
                    .gt. 0.013)) then

                     call angleDistribution(posInProbe)

                    end if

                end if


                ! Only writes to array if particle is within bounds of the image
                if ((posInProbexPx .lt. xPx) .and. (posInProbexPx .gt. 0)) then

                    image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + 1D0
                    
                    if (zImage) then
                    
                        image2(posInProbeyPx,posInProbexPx,t) = image2(posInProbeyPx,posInProbexPx,t) + 1D0

                    end if

                end if

            end do

        end subroutine getPosInProbe

        ! Writes out image array into a sequence of images
        subroutine writeImage(image, xPx, zPx, NumberOfTimePoints)
            implicit none

            real(kind=r14), intent(inout), dimension(:,:,:,:) :: image
            integer :: t, i, j, k, NumberOfTimePoints, xPx, zPx
            character(30) :: fileName

            print *, 'entering write'

            do k = 1, 3
            
                do t = 1, NumberOfTimePoints

                    if (k == 1) then
                    
                        write(fileName,'("Images/Image",I2,".txt")')t

                    else if (k == 2) then

                        write(fileName,'("Images2/Image",I2,".txt")')t

                    else

                        write(fileName,'("Images3/Image",I2,".txt")')t

                    end if

                    open(unit=20+t,file=filename)


                    do i = 1, zPx

                        do j = 1, xPx

                            write(20+t,'(ES12.5)',advance='no') image(i,j,t,k)

                        end do

                        write(20+t,*)

                    end do

                end do

            end do

        end subroutine writeImage

        ! for testing purposes only. Be careful using this in main code as you may need to adjust what
        ! variables are being input and output
        subroutine getPosInProbeTest(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, yPx, zPx, t0, probeStart, &
             tStep, particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)
            implicit none

            integer, intent(in) :: NumberOfTimePoints, startTimePoint, endTimePoint, xPx, yPx, zPx
            real(kind=r14), intent(in) :: probeStart, tStep, particleSpeed, pxMmRatio, t0
            real(kind=r14), dimension(3), intent(in) :: particleVector, particleStartPos, sheetDimensions
            integer :: t, posInProbexPx, posInProbeyPx, posInProbezPx, sheetCentrePx
            real(kind=r14) :: currentTime
            real(kind=r14), dimension(3) :: posInProbe

            ! Loops from entry timepont to exit timepoint to avoic wasting cycles when particle is not within sheet
            do t = startTimePoint, endTimePoint

                ! currentTime refers to the time it has taken the particle to travel from its starting point
                ! to the point in space at the given timepoint
                currentTime = probeStart + (t-1)*tStep - t0

                ! Real space position for particle
                posInProbe(1) = particleStartPos(1) + (particleVector(1)*particleSpeed*currentTime)
                posInProbe(2) = particleStartPos(2) + (particleVector(2)*particleSpeed*currentTime)
                posInProbe(3) = particleStartPos(3) + (particleVector(3)*particleSpeed*currentTime)
                
                if ((posInProbe(3) .ge. (0.021D0 - sheetDimensions(3)/2)) .and. ((posInProbe(3) .lt. &
                 (0.021D0 + sheetDimensions(3)/2)) )) then
                    if((posinProbe(2) .gt. (sheetDimensions(2)/-2)) .and. (posinProbe(2) .lt. (sheetDimensions(2)/2))) then
                        
                        ! Relative pixel position for particle
                        ! Note: the subtraction at the end of each statement alters the position of the
                        ! particle within the image array.
                        ! Altering the x-postion by half the width of the image centres the beam
                        ! Similarly, the z position can be altered however a factor of 1.3 was found to centre the sheet
                        ! within the middle of the image quite well
                        posInProbexPx = abs(ceiling(posInProbe(1)/pxMmRatio) - floor(real(xPx/2)))
                        posInProbeyPx = abs(ceiling(posInProbe(2)/pxMmRatio) - floor(real(yPx/2)))
                        posInProbezPx = abs(ceiling(posInProbe(3)/pxMmRatio) - floor(real(zPx/1.3)))

                        ! Only writes to array if particle is within bounds of the image
                        if ((posInProbexPx .lt. xPx) .and. (posInProbexPx .gt. 0)) then

                            !image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + 1
                            image2(posInProbeyPx,posInProbexPx,t) = image2(posInProbeyPx,posInProbexPx,t) + 1


                        end if
                    end if
                end if

            end do

        end subroutine getPosInProbeTest

        ! for testing purposes only: produces a set of images viewed along z axis
        subroutine writeImageTest(xPx, yPx, NumberOfTimePoints)
            implicit none

            integer :: t, i, j, NumberOfTimePoints, xPx, yPx
            character(30) :: fileName

            print *, 'entering write'

            do t = 1, NumberOfTimePoints


                write(fileName,'("Images2/Image",I2,".txt")')t
                open(unit=20+t,file=filename)

                do i = 1, yPx

                    do j = 1, xPx

                        write(20+t,'(i7)',advance='no') image2(i,j,t)

                    end do

                    write(20+t,*)

                end do

            end do

        end subroutine writeImageTest

        subroutine convim(imin,nx,ny,gaussdev,imout)
            !Convolutes input image imin with a gaussian of st. dev. gaussdev (in pixels), to produce imout.
                implicit none
                real(kind=r14), dimension(:,:), intent(in) :: imin(nx,ny)
                real(kind=r14), dimension(:,:), intent(out) :: imout(nx,ny)
                integer, intent(in) :: nx,ny
                double precision, intent(in) :: gaussdev
                
                double precision, dimension(:), allocatable :: gauss
                double precision, dimension(:,:) :: immid(nx,ny)
                
                double precision, parameter :: sqrt2 = 1.414213562d0
                
                integer :: gsize
                double precision :: cent
                integer :: icent
                integer :: j,k,l
                double precision :: lim1,lim2
                double precision :: dblej
                        
            !Calculate size of gaussian array required. This can just be a 1D gaussian since 2D Gaussian is separable
            ! and hence convolution can be applied stepwise for each dimension.
                if ((gaussdev .eq. 0.0d0) .and. (sum(imin) .gt. 0)) then
                        imout = imin
                else
                        gsize = floor(6.0d0*gaussdev)+1
                        if (mod(gsize,2) .eq. 0) gsize = gsize + 1 !Gaussian array will always have an odd-numbered size
                        allocate(gauss(gsize))
            !Calculate gaussian:
                        cent  = (dble(gsize)+1.0d0)*0.50d0
                        icent = int(cent)
                        
                        do j = 1,gsize
                            if (j .eq. icent) then
                                    lim1 = 0.5d0/(sqrt2*gaussdev)
                                    gauss(j) = erf(lim1)
                            else
                                    dblej = dble(j)
                                    lim1 = (abs(dblej - cent) +0.50d0)/(sqrt2*gaussdev)
                                    lim2 = lim1 - 0.5d0/(sqrt2*gaussdev)
                                    gauss(j) = erf(lim1)-erf(lim2)
                            end if
                        end do
            
            !Convolute along x.
                        do k = 1, ny
                            do j = 1, nx
                                    immid(j,k) = 0.0d0
                                    do l = 1, gsize
                                        if (((j+l-icent) .gt. 0) .and. ((j+l-icent) .le. nx)) then
                                                immid(j,k) = immid(j,k) + imin(j+l-icent,k)*gauss(l)
                                        end if
                                    end do
                            end do
                        end do
            !Convolute along y.
                        do j = 1, nx
                            do k = 1, ny
                                    imout(j,k) = 0.0d0
                                    do l = 1, gsize
                                        if (((k+l-icent) .gt. 0) .and. ((k+l-icent) .le. ny)) then
                                                imout(j,k) = imout(j,k) + immid(j,k+l-icent)*gauss(l)
                                        end if
                                    end do
                            end do
                        end do
            !Normalize output image to have same intensity as input image:
                        

                    imout = imout*sum(imin)/sum(imout)

                        
                end if
            end subroutine convim
        

end module imaging