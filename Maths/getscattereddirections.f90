module class_getscattereddirections
    implicit none

    private
    double precision :: pi = 3.14159265359

    public :: getscattereddirections

    contains

    function getscattereddirections result(scattereddirection)

        integer :: i
        logical :: hit
        double precision :: rand1, rand2, rand3, rand4, theta, phi, sin2theta
        double precision, dimension(:) :: scattereddirection(3)

        call random_number(rand1)

        phi = rand1 * pi / 2

        call random_number(rand2)
        ! see paper by J. Greenwood, Vacuum 2002 for explanation of how to generate angle distribution, not as simple as cos(theta)!
        theta = asin(SQRT(rand2))

        scattereddirection(3) = (cos(theta))

        call random_number(rand3)
        
        if (rand3 .gt. 0.5) then
            scattereddirection(1) = (cos(phi))*(sin(theta))

        else 
            scattereddirection(1) = -(cos(phi))*(sin(theta))

        end if

        call random_number(rand4)

        if (rand4 .gt. 0.5) then
            
            scattereddirection(2) = (sin(phi))*(sin(theta))

        else 

            scattereddirection(2) = -(sin(phi))*(sin(theta)) 
        end if

    end function getscattereddirections

end module class_getscattereddirections