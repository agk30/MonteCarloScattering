module sheetIntersection
    use mathConstants

    contains

        ! Returns intersection coordinates for top, bottom, front and back planes of the sheet
        subroutine getSheetIntersection (particleVector, startPos, sheetCentre, sheetDimensions, intersection)

            double precision, intent(in), dimension(3) :: particleVector, startPos, sheetCentre, sheetDimensions
            ! intersection contains x y z coordinates of intersection for each plane, in order of top, bottom,
            ! front and back from 1 - 4
            double precision, intent(out), dimension(4,3) :: intersection
            double precision, dimension(4,3) ::  pointOnPlane, planeVector
            double precision :: topOfFraction, bottomOfFraction, D
            integer :: row

            !since the parametric equation of a line can be given from a vector and a point, the eqution may be expressed as
            ! x = x0 + at or x = startPos(1) + (vector(1)*t)
            !equation for a plane is also required, using input values. the plane equation takes the form of Ax + By + Cz + D = 0
            !A B and C correspond to x y and z components of the vector orthogonal to the plane. D = -(A*x0 + By0 + Cz0)
            !intersection of a line with a plane can be found with equation x = x0 - (a*(Ax0 + By0 + Cz0 + D)/(a*A + b*B + c*C))
            !                                                               y = y0 - (a*(Ax0 + By0 + Cz0 + D)/(a*A + b*B + c*C))
            !                                                               z = z0 - (a*(Ax0 + By0 + Cz0 + D)/(a*A + b*B + c*C))
            !top and bottom share same orthogonal vector, as do front and back

            !TODO add in some way of inputting desried plane properties

            planeVector = 0D0 ; pointOnPlane = 0D0
            ! top plane
            planeVector(1,2) = 1D0 ; pointOnPlane(1,2) = sheetCentre(2) + (sheetDimensions(2)/2D0)
            ! bottom plane
            planeVector(2,2) = 1D0 ; pointOnPlane(2,2) = sheetCentre(2) - (sheetDimensions(2)/2D0)
            ! front plane
            planeVector(3,3) = 1D0 ; pointOnPlane(3,3) = sheetCentre(3) + (sheetDimensions(3)/2D0)
            ! back plane
            planeVector(4,3) = 1D0 ; pointOnPlane(4,3) = sheetCentre(3) - (sheetDimensions(3)/2D0)

            do row = 1, 4

                D = -((planeVector(row,1)*pointOnPlane(row,1)) + (planeVector(row,2)*pointOnPlane(row,2)) &
                 + (planeVector(row,3)*pointOnPlane(row,3)))
                
                topOfFraction = (planeVector(row,1)*startPos(1)) + (planeVector(row,2)*startPos(2)) &
                 + (planeVector(row,3)*startPos(3)) + D
                bottomOfFraction = (planeVector(row,1)*particleVector(1)) + (planeVector(row,2)*particleVector(2)) &
                 + (planeVector(row,3)*particleVector(3))
                
                intersection(row,1) = startPos(1) - (particleVector(1)*(topOfFraction/bottomOfFraction))
                intersection(row,2) = startPos(2) - (particleVector(2)*(topOfFraction/bottomOfFraction))
                intersection(row,3) = startPos(3) - (particleVector(3)*(topOfFraction/bottomOfFraction))
            
            end do

        end subroutine getSheetIntersection

        ! Returns a logical array showing which faces of the sheet were intersected if any
        subroutine withinSheet (intersection, sheetCentre, sheetDimensions, within)

            double precision, intent(in), dimension(4,3) :: intersection
            double precision, intent(in), dimension(3) :: sheetCentre, sheetDimensions
            ! 1 2 3 and 4 correspond to top, bottom, front and back
            logical, intent(out), dimension(4) :: within
            double precision :: sheetFront, sheetBack, sheetTop, sheetBottom
            integer :: i

            within = .FALSE.

            ! Establishes z coordinates for front and back faces, then y coordinates for top and bottom faces
            sheetFront = sheetCentre(3) + (sheetDimensions(3)/2D0)
            sheetBack = sheetCentre(3) - (sheetDimensions(3)/2D0)
            sheetTop = sheetCentre(2) + (sheetDimensions(2)/2D0)
            sheetBottom = sheetCentre(2) - (sheetDimensions(2)/2D0)

            ! Top face
            if ((intersection(1,3) .gt. sheetBack) .and. (intersection(1,3) .lt. sheetFront)) then

                within(1) = .TRUE.

            end if

            ! Bottom face
            if ((intersection(2,3) .gt. sheetBack) .and. (intersection(2,3) .lt. sheetFront)) then
                
                within(2) = .TRUE.

            end if

            ! Front face
            if ((intersection(3,2) .gt. sheetBottom) .and. (intersection(3,2) .lt. sheetTop)) then

                within(3) = .TRUE.

            end if

            ! Back face
            if ((intersection(4,2) .gt. sheetBottom) .and. (intersection(4,2) .lt. sheetTop)) then

                within(4) = .TRUE.

            end if

        end subroutine withinSheet

        ! Returns time of entry and exit from the sheet
        subroutine getIntersectionTIme(hitsSheet, intersection, particleStartPos, particleVector, &
             particleSpeed, particleTime, entryTime, exitTime)
             
            logical, dimension(4), intent(in) :: hitsSheet
            double precision, dimension(4,3), intent(in) :: intersection
            double precision, dimension(3), intent(in) :: particleStartPos, particleVector
            double precision, intent(in) :: particleSpeed, particleTime
            double precision, dimension(4) :: intersectionTime
            double precision, intent(out) :: entryTime, exitTime
            integer :: i

            intersectionTime = 0

            ! Loops across the number of faces of the sheet according to established index
            ! of 1, 2, 3 and 4 corresponding to top, bottom, front and back faces
            do i = 1, 4

                if (hitsSheet(i)) then

                    intersectionTime(i) = particleTime + (abs(particleStartPos(3) - intersection(i,3))) &
                     / abs(particleVector(3)*particleSpeed)

                end if

            end do

            ! Finds the maximum time of intersection and minimum
            entryTime = minval(intersectionTime, mask = intersectionTime .gt. 0)
            exitTime = maxval(intersectionTime)

        end subroutine getIntersectionTime

end module sheetIntersection