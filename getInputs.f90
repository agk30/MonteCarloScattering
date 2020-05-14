module getInputs
    use mathConstants

    contains

        ! Loads input parameters into the main section of code, MCScattering.f90
        ! See inputs.inp for details on parameter definitions
        ! TODO pass over hash table instead of individual variables
         subroutine loadInputs (ncyc, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos, colPos, skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering)
            implicit none

            integer, parameter :: r14 = selected_real_kind(14,30)
            integer, intent(out) :: ncyc
            real(kind=r14), intent(out) :: x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
            real(kind=r14), intent(out) :: colPos, skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed
            logical, intent(out) :: scattering

            open(unit=11,file='inputs.inp')

            read(11,*) x0
            read(11,*) aMax
            read(11,*) aMin
            read(11,*) h
            read(11,*) s
            read(11,*) dist
            read(11,*) pulseLength
            read(11,*) scattering
            read(11,*) mass
            read(11,*) temp
            read(11,*) ncyc
            read(11,*) skimPos
            read(11,*) valvePos
            read(11,*) colPos
            read(11,*) skimRad
            read(11,*) valveRad
            read(11,*) colRad
            read(11,*) sheetCentre
            read(11,*) halfSheetHeight
            read(11,*) sheetWidth
            read(11,*) probeStart
            read(11,*) probeEnd
            read(11,*) tStep
            read(11,*) pxMmRatio
            read(11,*) maxSpeed

        end subroutine loadInputs
        
end module getInputs