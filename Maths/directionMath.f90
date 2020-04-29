module class_directions
    implicit none

    private
    double precision :: pi = 3.14159265359

    contains

    !Randomly pick a points from a unit radius circle.
    function discPick result(x,y)
        implicit none
        double precision :: x, y, rand1, rand2

        !Random number for distance point is from centre of unit circle the value is square rooted so that the points alone this line will result in an even distribution of point on the unit circle rather than at a higher density at the centre
        call random_number(rand1)
        rand1 = sqrt(rand1)

        !Random number called for the angle this point is on the circle. The number is then converted into an angle in radians.
        call random_number(rand2)
        rand2 = 2*pi*rand2

        !trigonometry to turn the distance from centre of circle and angle into x and y coordinates on the unit circle
		x = rand1*cos(rand2)
		y = rand1*sin(rand2)

    end function discPick

    !calculates the gradient and intercept of the line which connect the point on the skimmer and valve
    !note x and y are used here to mean y=mx+c not in reference to chamber coordinates (z chamber is x here)
    function fitLine (y2, x2, y1, x1) result(m, c)
        implicit none

        double precision :: m, c, y2, x2, y1, x1

        m = (y2-y1)/(x2-x1)
        c = y2 - (m*x2)

    end function fitLine

    !Calcualtes unit vectors from gradients of the lines in the x and y directions. 
    function unitVector (mx, my) result(vx, vy, vz, a)
        implicit none

        double precision :: mx, my, vx, vy, vz, a

        a = sqrt(mx**2 + my**2 + 1**2)
        vx = mx/a
        vy = my/a

        if (z2.gt.z1) then

            vz = 1/a 

            else

            vz = -1/a 

        end if

    end function unitVector

end module class_directions
