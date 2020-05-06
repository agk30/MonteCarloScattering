include "Constants/mathConstants.f90"
include "getInputs.f90"
include "Maths/getSpeeds.f90"
include "Maths/getDirections.f90"

program MCScattering
    use getInputs
    use getSpeeds
    use getDirections
    use mathConstants
 
    implicit none

    integer :: ncyc
    real(kind=r14) :: x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
    real(kind=r14) :: colPos, skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed
    logical :: scattering

    integer :: i
    real(kind=r14) :: mostLikelyProbability, speed, scatteredSpeed
    !ingoingUnit and scatteredUnit are the unit vectors of the ingoing and scattering particles respectively
    real(kind=r14), dimension(3) :: ingoingUnit, scatteredUnit
    ! Each array contains the actual vectors for the particle, ingoing and scattered respectively
    real(kind=r14), dimension(3) :: ingoingVector, scatteredVector

    call random_seed

    ! Loads parameters from input file into main body of code for use in other functions
    call loadInputs(ncyc, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos, colPos, skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering)

    if (scattering == .TRUE.) then

        mostLikelyProbability = MBMostLikely(temp, mass)

    end if

    do i = 1, ncyc

        call ingoingSpeed(x0, aMax, aMin, h, s, dist, pulseLength, speed)
        caLL ingoingDirection(valveRad, valvePos, skimRad, skimPos, colRad, colPos, ingoingUnit)
        ingoingVector = ingoingUnit*speed

        if (scattering == .TRUE.) then

            call MBSpeed(maxSpeed, temp, mass, mostLikelyProbability, scatteredSpeed)
            call thermalDesorptionDirection(scatteredUnit)
            scatteredVector = scatteredUnit*scatteredSpeed

        end if

    end do

end program MCScattering