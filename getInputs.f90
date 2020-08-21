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
               writeImages, fullSim)
            implicit none

            integer, parameter :: r14 = selected_real_kind(14,30)
            integer, intent(out) :: ncyc, xPx, zPx, ksize, polyOrder
            real(kind=r14), intent(out) :: incidenceAngle, x0, aMax, aMin, &
            h, s, dist, pulseLength, mass, temp, skimPos, valvePos, gaussDev, massMol, energyTrans, surfaceMass, exitAngle
            real(kind=r14), intent(out) :: colPos, skimRad, valveRad, colRad, sheetCentre, &
             halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed
            logical, intent(out) :: scattering, testMods, writeImages, fullSim

            open(unit=11,file='Inputs/inputs.inp')
            
            read(11,*) xPx
            read(11,*) zPx
            read(11,*) incidenceAngle
            read(11,*) x0
            read(11,*) aMax
            read(11,*) aMin
            read(11,*) h
            read(11,*) s
            read(11,*) dist
            read(11,*) scattering
            read(11,*) mass
            read(11,*) massMol
            read(11,*) energyTrans 
            read(11,*) surfaceMass
            read(11,*) exitAngle
            read(11,*) temp
            read(11,*) ncyc
            read(11,*) maxSpeed
            read(11,*) skimPos
            read(11,*) valvePos
            read(11,*) colPos
            read(11,*) skimRad
            read(11,*) valveRad
            read(11,*) colRad
            read(11,*) sheetCentre
            read(11,*) halfSheetHeight
            read(11,*) sheetWidth
            read(11,*) pxMmRatio
            read(11,*) pulseLength
            read(11,*) probeStart
            read(11,*) probeEnd
            read(11,*) tStep
            read(11,*) gaussDev
            read(11,*) ksize
            read(11,*) polyOrder
            read(11,*) fullSim
            read(11,*) testMods
            read(11,*) writeImages

        end subroutine loadInputs
        
end module getInputs