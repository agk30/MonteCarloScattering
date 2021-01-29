module speeds
    use mathConstants

    contains

        ! Calculates speed of ingoing particle based on cumulative integral function
        subroutine ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, speed, t0)

            ! variables relating to cumulative integral function of arrival times
            double precision, intent(in) :: x0, aMax, aMin, h, s
            double precision, intent(in) :: dist, pulseLength
            double precision, intent(out) :: speed, t0
            double precision :: t, x, arrivalTime

            ! Calculate random time of creation
            call random_number(t)
            t0 = (t*pulseLength) - (pulseLength/2.0)

            ! CaLculate TOF based on cumulative integral function from real data anf fit by Origin.
            ! Function in Origin is called Logistics5.
            call random_number(x)
            arrivalTime = x0/(((aMax-aMin)/(x-aMin))**(1.0D0/s)-1.0D0)**(1.0D0/h)
            speed = dist/(arrivalTime*1D-6)
        end subroutine ingoing_speed

        ! Calculates speed based on Maxwell-Boltzmann Distribution of speeds
        subroutine MB_speed(maxSpeed, temp, mass, mostLikelyProbability, scatteredSpeed)

            logical :: hit
            double precision, intent(in) :: maxSpeed, temp, mass, mostLikelyProbability
            double precision, intent(inout) :: scatteredSpeed
            double precision :: rand1, rand2, probability, normalisedProbability

            hit = .FALSE.

            do while (hit .eqv. .FALSE.)
                call random_number(rand1)
                scatteredSpeed = rand1*maxSpeed

                probability = MB_probability(temp, scatteredSpeed, mass)

                ! Calculates the probability of the speed with respect to the most probable speed equalling 1.
                ! The Maxwell-Boltzmann distribution is already normalised to 1, meaning that the sum of all
                ! probabilities from zero to infinity will equal 1.
                ! It is possible to avoid this step, however, it would take a very long time to
                ! achieve a hit due to the small value of probability.
                normalisedProbability = probability/mostLikelyProbability

                call random_number(rand2)

                if (normalisedProbability .gt. rand2) then
                    hit = .TRUE.
                end if
            end do
        end subroutine MB_speed

        subroutine lorentzian_distribution(speed)
            implicit none

            double precision :: speed, gamma, rand1

            call random_number(rand1)

            gamma = 40.0D0

            speed = gamma*tan(pi*(rand1-0.5D0))

        end subroutine lorentzian_distribution

        ! This method apparently is very efficient, however it generates two Gaussian distributed numbers at a time
        ! Make sure this is incorporated somehow to avoid wasting cycles
        subroutine gaussian_distribution(mean, sigma, z1, z2)
            implicit none

            double precision :: rand1, rand2, v1, v2, rSquared, z1, z2, mean, sigma

            do
                call random_number(rand1)
                call random_number(rand2)

                v1 = (2.0*rand1) - 1.0
                v2 = (2.0*rand2) - 1.0

                rSquared = (v1**2.0) + (v2**2.0)

                if (rSquared .lt. 1) then
                    z1 = v1*SQRT((-2.0*log(rSquared))/rSquared)
                    z2 = v2*SQRT((-2.0*log(rSquared))/rSquared)

                    z1 = mean + sigma*z1
                    z2 = mean + sigma*z2

                    EXIT
                end if
            end do
        end subroutine

        ! Finds probability of particle travelling at given speed
        function MB_probability(temp, speed, mass) result(probability)

            double precision :: part1, part2, part3, speed, temp, mass, probability

            !part 1, 2, 3 correspond to individual parts of the maxwell-boltzmann distribution
            ! formula for calculating probability of a given speed
            part1 = 4.0D0*pi*speed*speed
            part2 = (mass/(2*pi*boltzmannConstant*temp))**(3.0D0/2.0D0)
            part3 = DEXP((-mass*speed*speed)/(2.0D0*boltzmannConstant*temp))

            probability = part1*part2*part3

        end function MB_probability

        ! Finds the most probable speed and its probability to use in normalisation 
        function MB_most_likely (temp, mass) result(mostLikelyProbability)

            double precision :: temp, mass, mostProbableSpeed, mostLikelyProbability

            mostProbableSpeed = sqrt((2.0D0*boltzmannConstant*temp)/mass)
            mostLikelyProbability = MB_probability(temp, mostProbableSpeed, mass)

        end function MB_most_likely

        subroutine deflection_angle(ingoing, outgoing, deflectionAngle)
            implicit none

            double precision, intent(in), dimension(3) :: ingoing, outgoing
            double precision, intent(out) :: deflectionAngle

            ! since this dot product finds the angle between the two vectors, it necessarily finds the deflection angle
            ! this is because the vectors are assumed to begin at the same point, and this is not the case with
            ! the ingoing and outgoing vectors, so the step where the angle is subtracted from 180 is not necessary
            deflectionAngle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360.0D0/(2*pi))

        end subroutine deflection_angle

        subroutine soft_sphere_speed(mass, internalRatio, surfaceMass, initialSpeed, deflectionAngle, finalSpeed)
            implicit none

            double precision :: initialEnergy, finalEnergy, massRatio, surfaceMass &
            ,part1, part2, part3, part4, part5, internalRatio, energyDiff, mass
            double precision, intent(in) :: initialSpeed, deflectionAngle
            double precision, intent(out) :: finalSpeed
            
            massRatio = mass/surfaceMass*1000.0D0
            initialEnergy = 0.5D0 * mass * initialSpeed * initialSpeed

            part1 = (2.0D0*massRatio)/((1+massRatio)**2.0D0)

            part2 = 1 + (massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2.0))
        
            part3 = cos(deflectionAngle*(2*pi/360.0D0))

            part4 = SQRT(1 - (massRatio*massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2)) - internalRatio*(massRatio + 1))
        
            part5 = internalRatio*((massRatio + 1.0)/(2.0*massRatio))

            energyDiff = part1*(part2 - (part3 * part4) + part5) * initialEnergy

            finalEnergy = initialEnergy - energyDiff

            finalSpeed = SQRT(2*finalEnergy/mass)

        end subroutine soft_sphere_speed

end module speeds