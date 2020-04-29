program sheetIntersection
    implicit none
    integer :: i, j
    double precision :: xValve, yValve, zValve, xVector, yVector, zVector, xPlusVector, yPlusVector, zPlusVector
    !a, b, and c are (x1 - x0), (y1 -y0) and (z1 - z0) in the parametric equation for the line
    double precision :: a, b, c, s, ncyc, sheetLength, sheetHeight, sheetWidth, sheetCentre
    double precision :: xFrontIntersection, yFrontIntersection, zFrontIntersection
    double precision :: xBackIntersection, yBackIntersection, zBackIntersection
    double precision :: xTopIntersection, yTopIntersection, zTopIntersection
    double precision :: xBottomIntersection, yBottomIntersection, zBottomIntersection

    open(unit=11,file='sheetintersection.inp')
    !note: in getdirectionsv2, the file is called getdirectionsv2.txt, this is just a placeholder
    open(unit=12,file='getdirectionsv2.txt')
    open(unit=13,file='getdirections.inp')
    open(unit=14,file='sheetintersection.txt')
    write(14,*) "J, x pos, y pos, z pos, x vector, y vector, z vector"

    read(11,*) sheetLength
    read(11,*) sheetHeight
    read(11,*) sheetWidth
    read(11,*) sheetCentre
    read(12,*) !reads first line where labels are so it doesn't get in the way
    read(13,*) ncyc

    !establishes first point with which to draw line
    i = 1
    do j = 1, ncyc

        read(12,*) s, xValve, yValve, zValve, xVector, yVector, zVector
        

        !parametric equation for the line going through both points: x = x0 + (x1 - x0)t
        !if x1 = x0 + constant, then (x1 - x0) = constant, in this case, the unit vector

        a = xVector
        b = yVector
        c = zVector

        !TODO change parametric equations for actual sheet dimensions

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

        if ((abs(yFrontIntersection) .lt. sheetHeight/2 .and. abs(xFrontIntersection) .lt. sheetLength/2) .or. (abs(xTopIntersection) .lt. sheetlength/2 .and. zTopIntersection .gt. sheetCentre - sheetWidth/2 .and. zTopIntersection .lt. sheetCentre + sheetWidth/2) .or. (abs(xBottomIntersection) .lt. sheetlength/2 .and. zBottomIntersection .gt. sheetCentre - sheetWidth/2 .and. zBottomIntersection .lt. sheetCentre - sheetWidth)) then
         write(14,'(i7,6e20.8)') i, xValve, yValve, zValve, xVector, yVector, zVector
         i = i+1
        end if
    end do
print *, 'The number of selected trajectories is ' ,i
end program sheetIntersection