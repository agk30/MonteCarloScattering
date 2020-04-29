module class_BoltzmannDistribution
    implicit none

    private
    double precision :: pi = 3.14159265359
    double precision :: e = 2.7182818284
    double precision :: boltzmannConstant = 1.38064852E-23

    public :: MBProbability

    contains
        ! return type    function  name         (inputs, for function) result(output var name)
        double precision function MBProbability (T, speed, mass) result(probability)
            implicit none
            double precision :: part1, part2, part3, speed, T, mass

            !part 1, 2, 3 correspond to individual parts of the maxwell-boltzmann distribution formula for calculating probability of a given speed
            part1 = 4*pi*speed*speed
            part2 = (mass/(2*pi*boltzmannConstant*T))**(3.0/2.0)
            part3 = e**((-mass*speed*speed)/(2*boltzmannConstant*T))

            probability = part1*part2*part3

        end function MBProbability

end module class_BoltzmannDistribution