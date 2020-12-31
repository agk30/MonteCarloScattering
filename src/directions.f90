module directions
    use mathConstants
    use speeds

    contains

        ! Finds a trajectory that passes though skimmer and collimator
        subroutine ingoing_direction(valveRad, valvePos, skimRad, skimPos, colRad, colPos, ingoingUnitVector, valve)

            ! valve(1), valve(2) and valve(3) correspond to x y and z coordinates of valve position with which to draw line.
            ! Likewise with skimmer() and collimator()
            double precision, dimension(3) :: skimmer, collimator
            double precision, intent (in) :: valveRad, valvePos, skimRad, skimPos, colRad, colPos
            double precision, intent(out), dimension(3) :: ingoingUnitVector, valve
            double precision :: mx, my, cx, cy, z
            logical :: hit

            hit = .FALSE.

            ! Loops until a suitable trajectory is found
            do while (hit .eqv. .FALSE.)
                ! Finds random point on valve for particle origin
                call disc_pick(valve(1),valve(2))
                valve = valve*valveRad
                valve(3) = valvePos

                ! Finds random point on skimmer for particle to pass through
                call disc_pick(skimmer(1),skimmer(2))
                skimmer = skimmer*skimRad
                skimmer(3) = skimPos

                ! Finds linear properties for line between valve and collimator positions
                call fit_line(valve(1), valve(3), skimmer(1), skimmer(3), mx, cx)
                call fit_line(valve(2), valve(3), skimmer(2), skimmer(3), my, cy)

                ! Caclulates positoin of particle at collimator and decides if it passes through or not
                collimator(1) = mx*(valvePos - colPos) + cx
                collimator(2) = my*(valvePos - colPos) + cy
                collimator(3) = colPos
                ! z is the hypotenuse of the triangle formed using x and y coordinates of the particle's collimator position
                z = SQRT(collimator(1)**2 + collimator(2)**2)
                
                if (z .lt. colRad) then
                    call unit_vector(mx, my, ingoingUnitVector)
                    hit = .TRUE.
                end if
            end do

        end subroutine ingoing_direction

        !subroutine transverse_temp(temp, maxSpeed1D, mass, zPos, travelDistance, startTime, speed, startPoint, vector)
        subroutine transverse_temp(mean, sigma, zPos, travelDistance, startTime, speed, startPoint, vector)
            implicit none

            double precision, dimension(3) :: startPoint
            double precision, dimension(3), intent(inout) :: vector
            double precision :: zPos, travelDistance, startTime, speed, maxSpeed1D, temp, mass, MBSpeed, rand, z2, mean, sigma

            startTime = (startTime + abs(travelDistance/(vector(3)*speed)))
            startPoint(1) = startPoint(1) + (startTime*speed*vector(1))
            startPoint(2) = startPoint(2) + (startTime*speed*vector(2))
            startPoint(3) = zPos
    
            vector = vector*speed

            call random_number(rand)

            !call one_dim_MB_speed(maxSpeed1D, temp, mass, MBSpeed)
            if (rand .gt. 0.5) then
                call lorentzian_distribution(MBSpeed)
            else
                call gaussian_distribution(mean, sigma, MBSpeed, z2)
            end if

            vector(1) = MBSpeed

            call random_number(rand)

            if (rand .gt. 0.5) then
                call lorentzian_distribution(MBSpeed)
            else
                call gaussian_distribution(mean, sigma, MBSpeed, z2)
            end if

            vector(2) = MBSpeed

            vector = vector/norm2(vector)

        end subroutine transverse_temp

        subroutine rotation(oldVector, theta, newVector)

            double precision, intent(in), dimension(3) :: oldVector
            double precision, intent(out), dimension(3) :: newVector
            double precision, dimension(3,3) :: rotationMatrix
            double precision, intent(in) :: theta
            double precision :: costheta, sintheta
            
            costheta = cos(theta*((2*pi)/360.0D0))
            sintheta = sin(theta*((2*pi)/360.0D0))

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

        ! It must be noted that this cosine distribution of scattering angles is simply a function that fits data,
        ! not necessarily the absolute correct distribution one would expect to see
        subroutine cosine_distribution(cosinePower, scatteredDirection)

            implicit none
        
            double precision :: rand1, rand2, phi, theta, x
            integer, intent(in) :: cosinePower
            double precision, dimension(3), intent(out) :: scatteredDirection
        
            call random_number(rand1)
            call random_number(rand2)

            phi = rand1*2*pi
        
            x = rand2**(1.0/dble(cosinePower+1.0))
        
            theta = dacos(x)

            scatteredDirection(1) = sin(theta)*cos(phi)
            scatteredDirection(2) = sin(theta)*sin(phi)
            scatteredDirection(3) = cos(theta)
            
        end subroutine cosine_distribution

        !Randomly pick a points from a unit radius circle.
        subroutine disc_pick(x,y)

            double precision, intent(out) :: x, y
            double precision :: rand1, rand2

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

        end subroutine disc_pick

        !calculates the gradient and intercept of the line which connect the point on the skimmer and valve
        !note x and y are used here to mean y=mx+c not in reference to chamber coordinates (z chamber is x here)
        subroutine fit_line (y2, x2, y1, x1, m, c)

            double precision, intent(in) :: y2, x2, y1, x1 
            double precision, intent(out) :: m, c

            m = (y2-y1)/(x2-x1)
            c = y2 - (m*x2)

        end subroutine fit_line

        !Calcualtes unit vectors from gradients of the lines in the x and y directions. 
        subroutine unit_vector (mx, my, v)

            double precision, intent(in) :: mx, my
            ! v(1), v(2), v(3) correspond to vector component in x, y, z direction
            double precision, dimension(3), intent(out) :: v
            double precision :: magnitude

            magnitude = sqrt(mx**2.0D0 + my**2.0D0 + 1.0D0)
            v(1) = mx/magnitude
            v(2) = my/magnitude
            v(3) = -1.0D0/magnitude

        end subroutine unit_vector

end module directions
