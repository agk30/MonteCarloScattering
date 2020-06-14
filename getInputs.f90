module getInputs
    use mathConstants

    contains

        ! Loads input parameters into the main section of code, MCScattering.f90
        ! See inputs.inp for details on parameter definitions
        ! TODO pass over hash table instead of individual variables
         subroutine loadInputs (incidenceAngle, ncyc, x0, aMax, aMin, &
             h, s, dist, pulseLength, mass, temp, skimPos, valvePos, colPos, &
             skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth,&
              probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering)
            implicit none

            integer, parameter :: r14 = selected_real_kind(14,30)
            integer, intent(out) :: ncyc
            real(kind=r14), intent(out) :: incidenceAngle, x0, aMax, aMin, &
            h, s, dist, pulseLength, mass, temp, skimPos, valvePos
            real(kind=r14), intent(out) :: colPos, skimRad, valveRad, colRad, sheetCentre, &
             halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed
            logical, intent(out) :: scattering

            open(unit=11,file='Inputs/directionInputs.inp')
            open(unit=12,file='Inputs/chamberDimensions.inp')
            open(unit=13,file='Inputs/timingInputs.inp')

            
            read(11,*) incidenceAngle
            read(11,*) x0
            read(11,*) aMax
            read(11,*) aMin
            read(11,*) h
            read(11,*) s
            read(11,*) dist
            read(11,*) scattering
            read(11,*) mass
            read(11,*) temp
            read(11,*) ncyc
            read(11,*) maxSpeed

            read(12,*) skimPos
            read(12,*) valvePos
            read(12,*) colPos
            read(12,*) skimRad
            read(12,*) valveRad
            read(12,*) colRad
            read(12,*) sheetCentre
            read(12,*) halfSheetHeight
            read(12,*) sheetWidth
            read(12,*) pxMmRatio

            read(13,*) pulseLength
            read(13,*) probeStart
            read(13,*) probeEnd
            read(13,*) tStep

        end subroutine loadInputs
        
end module getInputs