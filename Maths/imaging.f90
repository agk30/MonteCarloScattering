module imaging    
    use mathConstants

    ! Array shared by entire class
    integer*8, dimension(:,:,:), allocatable :: image

    contains

            ! Uses the entry time and exit time to find the corresponding timepoint for imaging
        subroutine startEndTimePoints(NumberOfTimepoints, entryTime, exitTime, probeStart, probeEnd, tStep, startTimePoint, endTimePoint)
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

           ! If exit time is greater than end of probe, then end timepoint then particle is imaged all the way til the end of the probe time
            if (exitTime .gt. probeEnd) then

                endTimePoint = NumberOfTimepoints

            else

                ! Similarly to entry time, exit timepoint is found as the nearest timepoint above exit time + 1 for the same reasons
                endTimePoint = floor((exitTime - probeStart) / tStep) + 1

            end if

        end subroutine startEndTimePoints

        ! Finds the position a particle is in at any given timepoint, and finds its corresponding pixel position then writes it
        ! to the image array, adding intensity to that pixel region
        subroutine getPosInProbe(NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx, t0, probeStart, tStep, particleSpeed, pxMmRatio, particleVector, particleStartPos, sheetDimensions)
            implicit none

            integer, intent(in) :: NumberOfTimePoints, startTimePoint, endTimePoint, xPx, zPx
            real(kind=r14), intent(in) :: probeStart, tStep, particleSpeed, pxMmRatio, t0
            real(kind=r14), dimension(3), intent(in) :: particleVector, particleStartPos, sheetDimensions
            integer :: t, posInProbexPx, posInProbezPx, sheetCentrePx
            real(kind=r14) :: currentTime
            real(kind=r14), dimension(3) :: posInProbe

            ! Loops from entry timepont to exit timepoint to avoic wasting cycles when particle is not within sheet
            do t = startTimePoint, endTimePoint

                ! currentTime refers to the time it has taken the particle to travel from its starting point
                ! to the point in space at the given timepoint
                currentTime = probeStart + (t-1)*tStep - t0

                ! Real space position for particle
                posInProbe(1) = particleStartPos(1) + (particleVector(1)*particleSpeed*currentTime)
                posInProbe(3) = particleStartPos(3) + (particleVector(3)*particleSpeed*currentTime)

                ! Relative pixel position for particle
                ! Note: the subtraction at the end of each statement alters the position of the particle within the image array.
                ! Altering the x-postion by half the width of the image centres the beam
                ! Similarly, the z position can be altered however a factor of 1.3 was found to centre the sheet within the middle of the image quite well
                posInProbexPx = abs(ceiling(posInProbe(1)/pxMmRatio) - floor(real(xPx/2)))
                posInProbezPx = abs(ceiling(posInProbe(3)/pxMmRatio) - floor(real(zPx/1.3)))

                ! Only writes to array if particle is within bounds of the image
                if ((posInProbexPx .lt. xPx) .and. (posInProbexPx .gt. 0)) then

                    image(posInProbezPx,posInProbexPx,t) = image(posInProbezPx,posInProbexPx,t) + 1

                end if

            end do

        end subroutine getPosInProbe

        ! Writes out image array into a sequence of images
        subroutine writeImage(xPx, zPx, NumberOfTimePoints)
            implicit none

            integer :: t, i, j, NumberOfTimePoints, xPx, zPx
            character(30) :: fileName

            print *, 'entering write'

            do t = 1, NumberOfTimePoints


                write(fileName,'("Images/Image",I2,".txt")')t
                open(unit=20+t,file=filename)

                do i = 1, zPx

                    do j = 1, xPx

                        write(20+t,'(i7)',advance='no') image(i,j,t)

                    end do

                    write(20+t,*)

                end do

            end do

        end subroutine writeImage

end module imaging