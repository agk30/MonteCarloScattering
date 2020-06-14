module getDirections
    use mathConstants

    contains

        ! Finds a trajectory that passes though skimmer and collimator
        subroutine ingoingDirection(valveRad, valvePos, skimRad, skimPos, colRad, colPos, ingoingUnitVector, valve)

            ! valve(1), valve(2) and valve(3) correspond to x y and z coordinates of valve position with which to draw line.
            ! Likewise with skimmer() and collimator()
            real(kind=r14), dimension(3) :: skimmer, collimator
            real(kind=r14), intent (in) :: valveRad, valvePos, skimRad, skimPos, colRad, colPos
            real(kind=r14), intent(out), dimension(3) :: ingoingUnitVector, valve
            real(kind=r14) :: mx, my, cx, cy, z
            logical :: hit

            hit = .FALSE.

            ! Loops until a suitable trajectory is found
            do while (hit .eqv. .FALSE.)

                ! Finds random point on valve for particle origin
                call discPick(valve(1),valve(2))
                valve = valve*valveRad
                valve(3) = valvePos

                ! Finds random point on skimmer for particle to pass through
                call discPick(skimmer(1),skimmer(2))
                skimmer = skimmer*skimRad
                skimmer(3) = skimPos

                ! Finds linear properties for line between valve and collimator positions
                call fitline(valve(1), valve(3), skimmer(1), skimmer(3), mx, cx)
                call fitline(valve(2), valve(3), skimmer(2), skimmer(3), my, cy)

                ! Caclulates positoin of particle at collimator and decides if it passes through or not
                collimator(1) = mx*(valvePos - colPos) + cx
                collimator(2) = my*(valvePos - colPos) + cy
                collimator(3) = colPos
                ! z is the hypotenuse of the triangle formed using x and y coordinates of the particle's collimator position
                z = SQRT(collimator(1)**2 + collimator(2)**2)
                
                if (z .lt. colRad) then

                    call unitVector(mx, my, ingoingUnitVector)

                    hit = .TRUE.

                end if

            end do

        end subroutine ingoingDirection

        subroutine rotation(oldVector, theta, newVector)

            real(kind=r14), intent(in), dimension(3) :: oldVector
            real(kind=r14), intent(out), dimension(3) :: newVector
            real(kind=r14), dimension(3,3) :: rotationMatrix
            real(kind=r14), intent(in) :: theta
            real(kind=r14) :: costheta, sintheta
            
            costheta = cosd(theta)
            sintheta = sind(theta)

            ! this matrix is for roation about the y axis only. Rotation about any other axis will require a different matrix.
            rotationMatrix = 0
            rotationMatrix(1,1) = costheta
            rotationMatrix(1,3) = sintheta
            rotationMatrix(2,2) = 1
            rotationMatrix(3,1) = -sintheta
            rotationMatrix(3,3) = costheta

            ! multiplies the rotation matrix by the vector in question
            newVector = MATMUL(rotationMatrix, oldVector)

        end subroutine rotation

        ! Finds thermal desorption trajetory based on a cos(theta) distribution of scattering angles
        subroutine thermalDesorptionDirection(scatteredDirection)

            integer :: i
            logical :: hit
            real(kind=r14) :: rand1, rand2, rand3, rand4, theta, phi, sin2theta
            real(kind=r14), dimension(3) :: scatteredDirection
    
            call random_number(rand1)
    
            phi = rand1 * pi / 2
    
            call random_number(rand2)
            ! see paper by J. Greenwood, Vacuum 2002 for explanation of how to generate angle distribution,
            ! not as simple as cos(theta)!
            theta = asin(SQRT(rand2))
    
            scatteredDirection(3) = (cos(theta))
    
            call random_number(rand3)
            
            if (rand3 .gt. 0.5D0) then

                scatteredDirection(1) = (cos(phi))*(sin(theta))
    
            else 
                scatteredDirection(1) = -(cos(phi))*(sin(theta))
    
            end if
    
            call random_number(rand4)
    
            if (rand4 .gt. 0.5D0) then
                
                scatteredDirection(2) = (sin(phi))*(sin(theta))
    
            else 
    
                scatteredDirection(2) = -(sin(phi))*(sin(theta)) 
                
            end if
    
        end subroutine thermalDesorptionDirection

        !Randomly pick a points from a unit radius circle.
        subroutine discPick(x,y)

            real(kind=r14), intent(out) :: x, y
            real(kind=r14) :: rand1, rand2

            !Random number for distance point is from centre of unit circle the value is square rooted so that the points alone
            ! this line will result in an even distribution of point on the unit circle rather than at a higher density at the centre
            call random_number(rand1)
            rand1 = sqrt(rand1)

            !Random number called for the angle this point is on the circle. The number is then converted into an angle in radians.
            call random_number(rand2)
            rand2 = 2.0D0*pi*rand2

            !trigonometry to turn the distance from centre of circle and angle into x and y coordinates on the unit circle
            x = rand1*cos(rand2)
            y = rand1*sin(rand2)

        end subroutine discPick

        !calculates the gradient and intercept of the line which connect the point on the skimmer and valve
        !note x and y are used here to mean y=mx+c not in reference to chamber coordinates (z chamber is x here)
        subroutine fitLine (y2, x2, y1, x1, m, c)

            real(kind=r14), intent(in) :: y2, x2, y1, x1 
            real(kind=r14), intent(out) :: m, c

            m = (y2-y1)/(x2-x1)
            c = y2 - (m*x2)

        end subroutine fitLine

        !Calcualtes unit vectors from gradients of the lines in the x and y directions. 
        subroutine unitVector (mx, my, v)

            real(kind=r14), intent(in) :: mx, my
            ! v(1), v(2), v(3) correspond to vector component in x, y, z direction
            real(kind=r14), dimension(3), intent(out) :: v
            real(kind=r14) :: magnitude

            magnitude = sqrt(mx**2.0D0 + my**2.0D0 + 1.0D0)
            v(1) = mx/magnitude
            v(2) = my/magnitude
            v(3) = -1.0D0/magnitude

        end subroutine unitVector

end module getDirections
