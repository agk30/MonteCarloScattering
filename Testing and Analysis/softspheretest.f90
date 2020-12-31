subroutine getDeflectionAngle(ingoing, outgoing, deflectionAngle)
    implicit none

    double precision, intent(in), dimension(3) :: ingoing, outgoing
    double precision, intent(out) :: deflectionAngle
    double precision :: pi = 3.141592653589793D0

    ! since this dot product finds the angle between the two vectors, it necessarily finds the deflection angle
    ! this is because the vectors are assumed to begin at the same point, and this is not the case with
    ! the ingoing and outgoing vectors, so the step where the angle is subtracted from 180 is not necessary
    deflectionAngle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360.0D0/(2*pi))

end subroutine getDeflectionAngle

subroutine softSphereSpeed(mass, internalRatio, surfaceMass, initialSpeed, deflectionAngle, finalSpeed)
    implicit none

    double precision :: initialEnergy, finalEnergy, massRatio, particleMass, surfaceMass &
    ,part1, part2, part3, part4, part5, internalEnergyLoss, internalRatio, energyDiff, mass, pi
    double precision, intent(in) :: initialSpeed, deflectionAngle
    double precision, intent(out) :: finalSpeed

    pi = 3.141592653589793D0
    
    massRatio = mass/surfaceMass*1000.0D0
    initialEnergy = 0.5D0 * mass * initialSpeed * initialSpeed

    part1 = (2.0D0*massRatio)/((1+massRatio)**2.0D0)

    part2 = 1 + (massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2.0))

    part3 = cos(deflectionAngle*(2*pi/360.0D0))

    part4 = SQRT(1 - (massRatio*massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2)) - internalRatio*(massRatio + 1))

    !print *, part4, massRatio, deflectionAngle, internalRatio

    part5 = internalRatio*((massRatio + 1.0)/(2.0*massRatio))

    energyDiff = part1*(part2 - (part3 * part4) + part5) * initialEnergy

    finalEnergy = initialEnergy - energyDiff

    finalSpeed = SQRT(2*finalEnergy/mass)

    !print *, finalSpeed

end subroutine softSphereSpeed

program sstest
    implicit none

    double precision :: ingoing(3), outgoing(3)
    double precision :: deflectionAngle, mass, internalRatio, surfaceMass, initialSpeed, finalSpeed
    double precision :: initialEnergy, finalEnergy

    ingoing(1) = 1D0 ; ingoing(2) = 0D0 ; ingoing(3) = -1D0
    outgoing(1) = 1D0 ; outgoing(2) = 0D0 ; outgoing(3) = 1D0

    mass = 17D-3
    surfaceMass = 100D0
    initialSpeed = 1000D0
    internalRatio = 0D0

    initialEnergy = 0.5D0*mass*(initialSpeed**2)

    call getDeflectionAngle(ingoing, outgoing, deflectionAngle)

    print *, deflectionAngle

    call softSphereSpeed(mass, internalRatio, surfaceMass, initialSpeed, deflectionAngle, finalSpeed)

    finalEnergy = 0.5D0*mass*(finalSpeed**2)

    print *, initialEnergy, finalEnergy

end program sstest