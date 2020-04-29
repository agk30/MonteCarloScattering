program sheetIntersection
    implicit none
    double precision :: xValve, yValve, zValve, xPlusVector, yPlusVector, zPlusVector
    !a, b, and c are (x1 - x0), (y1 -y0) and (z1 - z0) in the parametric equation for the line
    double precision :: a, b, c, sheetLength, sheetHeight, sheetWidth, sheetCentre
    double precision :: xFrontIntersection, yFrontIntersection, zFrontIntersection
    double precision :: xBackIntersection, yBackIntersection, zBackIntersection
    double precision :: xTopIntersection, yTopIntersection, zTopIntersection
    double precision :: xBottomIntersection, yBottomIntersection, zBottomIntersection

    open(unit=11,file='sheetIntersection.txt')

    !sheet length is length along the x axis
    sheetLength = 0.05

    !sheet height is length along the y axis
    sheetHeight = 0.004

    !sheet witdh along z axis
    sheetWidth = 0.02125

    !sheet centre on z axis
    sheetCentre = 0.013685

    !establishes first point with which to draw line
    xValve = -0.139038E-03
    yValve = -0.222478E-03
    zValve = 0.153500E+00

    !establishes second point with which to draw line
    xPlusVector = xValve + 0.359194E-03
    yPlusVector = yValve + 0.622886E-02
    zPlusVector = zValve + -0.999981E+00

    !parametric equation for the line going through both points: x = x0 + (x1 - x0)t
    !if x1 = x0 + constant, then (x1 - x0) = constant, in this case, the unit vector

    a = 0.359194E-03
    b = 0.622886E-02
    c = -0.999981E+00

    !TODO replace raw numbers with variables

    !intersection with the front face of laser sheet; parametric equation for this plane is z - 0.2431 = 0
    xFrontIntersection = xValve - ((a*(zValve - 0.02431))/c)
    yFrontIntersection = yValve - ((b*(zValve - 0.02431))/c)
    zFrontIntersection = zValve - ((c*(zValve - 0.02431))/c)

    !intersection with the back face of laser sheet; parametric equation for the plane is z - 0.00306 = 0
    xBackIntersection = xValve - ((a*(zValve - 0.00306))/c)
    yBackIntersection = yValve - ((b*(zValve - 0.00306))/c)
    zBackIntersection = zValve - ((c*(zValve - 0.00306))/c)

    !intersection with the top face of laser sheet; parametric equation for this plane is y - 0.002 = 0
    xTopIntersection = xValve - ((a*(yValve - 0.002))/b)
    yTopIntersection = yValve - ((b*(yValve - 0.002))/b)
    zTopIntersection = zValve - ((c*(yValve - 0.002))/b)

    !intersection with the bottom face of laser sheet; parametric equation for this plane is y + 0.002 = 0
    xBottomIntersection = xValve - ((a*(yValve + 0.002))/b)
    yBottomIntersection = yValve - ((b*(yValve + 0.002))/b)
    zBottomIntersection = zValve - ((c*(yValve + 0.002))/b)

    if (abs(yFrontIntersection) .lt. sheetHeight/2 .and. abs(xFrontIntersection) .lt. sheetLength/2) then
        print *, 'enters front face'
    else if (abs(xTopIntersection) .lt. sheetlength/2 .and. zTopIntersection .gt. sheetCentre - sheetWidth/2 .and. zTopIntersection .lt. sheetCentre + sheetWidth/2) then
        print *, 'enters top face'
    else if (abs(xBottomIntersection) .lt. sheetlength/2 .and. zBottomIntersection .gt. sheetCentre - sheetWidth/2 .and. zBottomIntersection .lt. sheetCentre - sheetWidth) then
        print *, 'enters bottom face'
    else
        print *, 'does not enter sheet region'
    end if


    print *, xFrontIntersection, yFrontIntersection, zFrontIntersection
    print *, xBackIntersection, yBackIntersection, zBackIntersection
    print *, xTopIntersection, yTopIntersection, zTopIntersection
    print *, xBottomIntersection, yBottomIntersection, zBottomIntersection

end program sheetIntersection