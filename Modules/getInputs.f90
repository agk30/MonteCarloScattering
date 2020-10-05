module getInputs
    use mathConstants

    contains

        ! Loads input parameters into the main section of code, MCScattering.f90
        ! See inputs.inp for details on parameter definitions
        ! TODO pass over hash table instead of individual variables
         subroutine loadInputs (xPx, zPx, incidenceAngle, ncyc, x0, aMax, aMin, &
             h, s, dist, pulseLength, mass, massMol, energyTrans, surfaceMass, exitAngle, temp, skimPos, valvePos, colPos, &
             skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth,&
              probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering, gaussDev, ksize, polyOrder, testMods,&
               writeImages, fullSim, scatterFraction, scatterIntensity)
            implicit none

            integer, parameter :: r14 = selected_real_kind(14,30)
            integer, intent(out) :: ncyc, xPx, zPx, ksize, polyOrder
            real(kind=r14), intent(out) :: incidenceAngle, x0, aMax, aMin, &
            h, s, dist, pulseLength, mass, temp, skimPos, valvePos, gaussDev, massMol, energyTrans, surfaceMass, exitAngle
            real(kind=r14), intent(out) :: colPos, skimRad, valveRad, colRad, sheetCentre, &
             halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scatterFraction, scatterIntensity
            logical, intent(out) :: scattering, testMods, writeImages, fullSim

            open(unit=11,file="Inputs/experimentalInputs.inp")
            open(unit=12,file="Inputs/mathParameters.inp")
            open(unit=13,file="Inputs/imagingInputs.inp")

            read(11,*) skimPos
            read(11,*) valvePos
            read(11,*) colPos
            read(11,*) skimRad
            read(11,*) valveRad
            read(11,*) colRad
            read(11,*) sheetCentre
            read(11,*) halfSheetHeight
            read(11,*) sheetWidth
            read(11,*) pulseLength

            read(12,*) xPx
            read(12,*) zPx
            read(12,*) incidenceAngle
            read(12,*) x0
            read(12,*) aMax
            read(12,*) aMin
            read(12,*) h
            read(12,*) s
            read(12,*) dist
            read(12,*) mass
            read(12,*) massMol
            read(12,*) energyTrans
            read(12,*) surfaceMass
            read(12,*) exitAngle
            read(12,*) temp
            read(12,*) ncyc
            read(12,*) maxSpeed
            read(12,*) scatterFraction

            read(13,*) pxMmRatio
            read(13,*) probeStart
            read(13,*) probeEnd
            read(13,*) tStep
            read(13,*) gaussDev
            read(13,*) ksize
            read(13,*) polyOrder
            read(13,*) scattering
            read(13,*) fullSim
            read(13,*) testMods
            read(13,*) writeImages
            read(13,*) scatterIntensity        

        end subroutine loadInputs
        
end module getInputs