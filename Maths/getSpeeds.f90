module getSpeeds
    use mathConstants

    contains

        ! Calculates speed of ingoing particle based on cumulative integral function
        subroutine ingoingSpeed(x0, aMax, aMin, h, s, dist, pulseLength, speed, t0)

            ! variables relating to cumulative integral function of arrival times
            real(kind=r14), intent(in) :: x0, aMax, aMin, h, s
            real(kind=r14), intent(in) :: dist, pulseLength
            real(kind=r14), intent(out) :: speed, t0
            real(kind=r14) :: t, x, arrivalTime

            ! Calculate random time of creation
            call random_number(t)
            t0 = t*pulseLength

            ! CaLculate TOF based on cumulative integral function from real data anf fit by Origin.
            ! Function in Origin is called Logistics5.
            call random_number(x)
            arrivalTime = x0/(((aMax-aMin)/(x-aMin))**(1.0D0/s)-1.0D0)**(1.0D0/h)
            speed = dist/(arrivalTime*1D-6)

        end subroutine ingoingSpeed

        ! Calculates speed based on Maxwell-Boltzmann Distribution of speeds
        subroutine MBSpeed(maxSpeed, temp, mass, mostLikelyProbability, scatteredSpeed)

            logical :: hit
            real(kind=r14), intent(in) :: maxSpeed, temp, mass, mostLikelyProbability
            real(kind=r14), intent(inout) :: scatteredSpeed
            real(kind=r14) :: speed, rand1, rand2, probability, normalisedProbability

            hit = .FALSE.

            do while (hit .eqv. .FALSE.)

                call random_number(rand1)
                scatteredSpeed = rand1*maxSpeed

                probability = MBProbability(temp, scatteredSpeed, mass)

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

        end subroutine MBSpeed

        ! Finds probability of particle travelling at given speed
        function MBProbability (temp, speed, mass) result(probability)

            real(kind=r14) :: part1, part2, part3, speed, temp, mass, probability

            !part 1, 2, 3 correspond to individual parts of the maxwell-boltzmann distribution
            ! formula for calculating probability of a given speed
            part1 = 4.0D0*pi*speed*speed
            part2 = (mass/(2*pi*boltzmannConstant*temp))**(3.0D0/2.0D0)
            part3 = DEXP((-mass*speed*speed)/(2.0D0*boltzmannConstant*temp))

            probability = part1*part2*part3

        end function MBProbability

        ! Finds the most probable speed and its probability to use in normalisation 
        function MBMostLikely (temp, mass) result(mostLikelyProbability)

            real(kind=r14) :: temp, mass, mostProbableSpeed, mostLikelyProbability

            mostProbableSpeed = sqrt((2.0D0*boltzmannConstant*temp)/mass)
            mostLikelyProbability = MBProbability(temp, mostProbableSpeed, mass)

        end function MBMostLikely

        subroutine softSphereSpeed(mass, internalRatio, surfaceMass, initialSpeed, ingoing, outgoing, finalSpeed)
            implicit none

            real(kind=r14) :: initialEnergy, finalEnergy, massRatio, particleMass, surfaceMass &
            ,part1, part2, part3, part4, part5, deflectionAngle, internalEnergyLoss, internalRatio, energyDiff, mass
            real(kind=r14), intent(in) :: initialSpeed
            real(kind=r14), intent(in), dimension(3) :: ingoing, outgoing
            real(kind=r14), intent(out) :: finalSpeed

            ! since this dot product finds the angle between the two vectors, it necessarily finds the deflection angle
            ! this is because the vectors are assumed to begin at the same point, and this is not the case with
            ! the ingoing and outgoing vectors, so the step where the angle is subtracted from 180 is not necessary
            deflectionAngle = acosd(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing)))
            
            massRatio = mass/surfaceMass*1000.0D0
            initialEnergy = 0.5D0 * mass * initialSpeed * initialSpeed

            part1 = (2.0D0*massRatio)/((1+massRatio)**2.0D0)

            part2 = 1 + (massRatio*(sind(deflectionAngle)**2.0))
        
            part3 = cosd(deflectionAngle)
        
            part4 = SQRT(1 - (massRatio*massRatio*(sind(deflectionAngle)**2)) - internalRatio*(massRatio + 1))
        
            part5 = internalRatio*((massRatio + 1.0)/(2.0*massRatio))

            energyDiff = part1*(part2 - (part3 * part4) + part5) * initialEnergy

            finalEnergy = initialEnergy - energyDiff

            finalSpeed = SQRT(2*finalEnergy/mass)

            !print *, finalSpeed

        end subroutine softSphereSpeed

end module getSpeeds