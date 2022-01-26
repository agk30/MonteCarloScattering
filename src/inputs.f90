module inputs
    use mathConstants
    use m_config

    integer :: ncyc, ksize, polyOrder, cosinePowerTD, cosinePowerIS, xPx, zPx, surface_z
    double precision :: incidenceAngle, x0, aMax, aMin, h, s, dist, pulseLength, mass, temp, skimPos, valvePos
    double precision :: colPos, skimRad, valveRad, colRad, sheetCentreZ, halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, &
     pxMmRatio, maxSpeed, gaussDev, massMol, energyTrans, surfaceMass, exitAngle, scatterFraction, scatterIntensity, fLifeTime, &
      captureGateOpen, captureGateClose
    logical :: scattering, testMods, writeImages, fullSim
    logical :: correctDirection, normalRun, linux, hush
    character(200) :: imagePath, matrixPath, ifPath
        ! Fixed parameter variables
    logical :: fixedIngoingSpeed, fixedOutgoingSpeed, fixedStartPos, fixedScatterPos, fixedCreationTime, fixedScatterTime
    double precision :: speedIn, speedOut, creationTime, scatterTime
    double precision :: startx, starty, startz, scatterx, scattery, scatterz

    contains

        ! Loads input parameters into the main section of code, MCScattering.f90
        subroutine load_inputs
            implicit none
            logical :: file_exists
            type(CFG_t) :: inputs

            ! Build the variables used in the input file

            ! Experimental inputs
            call CFG_add(inputs, "skimPos", 0.1262D0 , "Position of Skimmer in z direction")
            call CFG_add(inputs, "valvePos", 0.1702D0 , "Position of Valve in z direction")
            call CFG_add(inputs, "colPos", 0.0787D0 , "Position of Collimator in z direction")
            call CFG_add(inputs, "skimRad", 0.001D0 , "Skimmer Radius")
            call CFG_add(inputs, "valveRad", 0.0015D0 , "Valve Radius")
            call CFG_add(inputs, "colRad", 0.0015D0 , "Collimator Radius")
            call CFG_add(inputs, "sheetCentre", 0.021D0 , "Distance from centre of sheet from to surface")
            call CFG_add(inputs, "halfSheetHeight", 0.002D0 , "Half height of laser sheet")
            call CFG_add(inputs, "sheetWidth", 0.0300D0 , "Full width of laser sheet")
            call CFG_add(inputs, "pulseLength", 10.0D-06 , "Discharge pulse length")
            call CFG_add(inputs, "surface_z", 294 , "z-offset of surface position in pixels from top of image")
            
            ! Imaging inputs
            call CFG_add(inputs, "pxMmRatio", 0.25D-03 , "Pixel to mm ratio for imaging")
            call CFG_add(inputs, "probeStart", 56.0D-06 , "Start of probe time (imaging time)")
            call CFG_add(inputs, "probeEnd", 160.0D-06 , "End of probe time (imaging time)")
            call CFG_add(inputs, "tStep", 1.0D-06 , "Time step between images")
            call CFG_add(inputs, "gaussDev", 2.0D0 , "Deviation parameter for gaussian blur routine")
            call CFG_add(inputs, "ksize", 27 , "Kernel size for SG routine")
            call CFG_add(inputs, "polyOrder", 3 , "Polynomial order for SG routine")
            call CFG_add(inputs, "scattering", .TRUE. , "Image scattering as well as ingoing beam?")
            call CFG_add(inputs, "fullSim", .TRUE. , "Image ingoing beam as well as scattering?")
            call CFG_add(inputs, "testMods", .FALSE. , "Including testing modules?")
            call CFG_add(inputs, "writeImages", .TRUE. , "Wirte images to files?")
            call CFG_add(inputs, "scatterIntensity", 3.0D0 , "Relative intensity of scattered signal to ingogin signal")
            call CFG_add(inputs, "fLifeTime", 700D-9 , "Fluorescent lifetime of the molecule")
            call CFG_add(inputs, "captureGateOpen", 350D-9 , "Time after probe fire at which imaging starts")
            call CFG_add(inputs, "captureGateClose", 10000D-9 , "Time after probe fire at which imaging ends")
            
            ! Mathematical Inputs
            call CFG_add(inputs, "xPx", 420 , "Number of image pixels in x direction")
            call CFG_add(inputs, "zPx", 420 , "Number of image pixels in z direction")
            call CFG_add(inputs, "incidenceAngle", 0D0 , "Incidence angle for beam")
            call CFG_add(inputs, "cosinePowerTD" , 1 , "Power of cosine weighted TD angle distribution (cos^n)")
            call CFG_add(inputs, "cosinePowerIS" , 4 , "Power of cosine weighted IS angle distribution (cos^n)")
            call CFG_add(inputs, "x0", 128.30344D0 , "Origin function parameter for speed generation")
            call CFG_add(inputs, "aMax", 1.00082D0 , "Origin function parameter for speed generation")
            call CFG_add(inputs, "aMin", -0.00851D0 , "Origin function parameter for speed generation")
            call CFG_add(inputs, "h", 24.98601D0 , "Origin function parameter for speed generation")
            call CFG_add(inputs, "s", 1.45285D0 , "Origin function parameter for speed generation")
            call CFG_add(inputs, "dist", 230.0D-03 , "Single-point LIF value-probe laser distance")
            call CFG_add(inputs, "mass", 2.8240519D-26 , "Molecule mass / kg")
            call CFG_add(inputs, "massMol", 0.017D0 , "Molecue mass / kg/mol")
            call CFG_add(inputs, "energyTrans", 0.0D0 , &
            "SS collision model - energy loss to internal surface motions")
            call CFG_add(inputs, "surfaceMass", 100D0 , "Effective liquid surface mass / g/mol")
            call CFG_add(inputs, "exitAngle", 0D0 , "Exit angle for SS model")
            call CFG_add(inputs, "temp", 298.0D0 , "Surface temp / K")
            call CFG_add(inputs, "ncyc", 10000000 , "Number of molcules to be sampled")
            call CFG_add(inputs, "maxSpeed", 3000.0D0 , "Max speed for MB speed calculation")
            call CFG_add(inputs, "scatterFraction", 0.5D0  , &
            "Fraction of molcules scattering in TD or IS. 0 for full TD, 1 for full IS")

            ! File Paths go here
            call CFG_add(inputs, "linux", .FALSE. , "Is this a linux system?")
            call CFG_add(inputs, "imagePath", "../Images/" , "Path to image files")
            call CFG_add(inputs, "matrixPath", "../SG Matrices/CC_027x027_003x003.dat",&
             "Path to SG matrix")
            call CFG_add(inputs, "ifPath", "../Real Images/09032021_Q12_Instrument Function_SUM_00",&
             "Path to instrument function image")

            ! Fixed parameters for molecules
            call CFG_add(inputs, "normalRun", .TRUE. , ".FALSE. to run with fixed parameters")

            call CFG_add(inputs, "fixedIngoingSpeed", .FALSE. , "Want to fix ingoing speed?")
            call CFG_add(inputs, "speedIn", 2000D0 , "Fixed ingoing speed")

            call CFG_add(inputs, "fixedOutgoingSpeed", .FALSE. , "Want to fix outoging speed?")
            call CFG_add(inputs, "speedOut", 1000D0 , "Fixed outgoing speed")

            call CFG_add(inputs, "fixedStartPos", .FALSE. , "Want to fix starting position?")
            call CFG_add(inputs, "startx", 0D0 , "Fixed start x position")
            call CFG_add(inputs, "starty", 0D0 , "Fixed start y position")
            call CFG_add(inputs, "startz", 0D0 , "Fixed start z position")

            call CFG_add(inputs, "fixedScatterPos", .FALSE. , "Want to fix scattering position?")
            call CFG_add(inputs, "scatterx", 0D0 , "Fixed scattering x position")
            call CFG_add(inputs, "scattery", 0D0 , "Fixed scattering y position")
            call CFG_add(inputs, "scatterz", 0D0 , "Fixed scattering z position")
            
            call CFG_add(inputs, "fixedCreationTime", .FALSE. , "Want to fix creation time?")
            call CFG_add(inputs, "creationTime", 0D0 , "Time at which molecules are created")

            call CFG_add(inputs, "fixedScatterTime", .FALSE. , "Want to fix scattering time?")
            call CFG_add(inputs, "scatterTime", 0D0 , "Time at which molecules scatter")

            call CFG_add(inputs, "hush", .FALSE., "shuts off console outputs except for filepath. This is to allow piping outputs to python code")

            ! Read .cfg file and parse for changes
            inquire(FILE="../Inputs/inputs.cfg", EXIST=file_exists)

            if (file_exists) then
                call CFG_read_file(inputs, "../Inputs/inputs.cfg")
            else
                inquire(FILE="inputs.cfg", EXIST=file_exists)
                if (file_exists) then
                    call CFG_read_file(inputs, "inputs.cfg")
                else
                    print '(a)', "Error: Input file not found. Running with default values."
                end if
            end if

            ! Parse any command line arguments
            call CFG_update_from_arguments(inputs)

            call CFG_get(inputs, "normalRun", normalRun)

            !Experimental inputs
            call CFG_get(inputs, "skimPos", skimPos)
            call CFG_get(inputs, "valvePos", valvePos)
            call CFG_get(inputs, "colPos", colPos)
            call CFG_get(inputs, "skimRad", skimRad)
            call CFG_get(inputs, "valveRad", valveRad)
            call CFG_get(inputs, "colRad", colRad)
            call CFG_get(inputs, "sheetCentre", sheetCentreZ)
            call CFG_get(inputs, "halfSheetHeight", halfSheetHeight)
            call CFG_get(inputs, "sheetWidth", sheetWidth)
            call CFG_get(inputs, "pulseLength", pulseLength)
            call CFG_get(inputs, "surface_z", surface_z)
        
            valvePos = valvePos! - 73E-3
        
            ! Imaging inputs
            call CFG_get(inputs, "pxMmRatio", pxMmRatio)
            call CFG_get(inputs, "probeStart", probeStart)
            call CFG_get(inputs, "probeEnd", probeEnd)
            call CFG_get(inputs, "tStep", tStep)
            call CFG_get(inputs, "gaussDev", gaussDev)
            call CFG_get(inputs, "ksize", ksize)
            call CFG_get(inputs, "polyOrder", polyOrder)
            call CFG_get(inputs, "scattering", scattering)
            call CFG_get(inputs, "fullSim", fullSim)
            call CFG_get(inputs, "testMods", testMods)
            call CFG_get(inputs, "writeImages", writeImages)
            call CFG_get(inputs, "scatterIntensity", scatterIntensity)
            call CFG_get(inputs, "fLifeTime", fLifeTime)
            call CFG_get(inputs, "captureGateOpen", captureGateOpen)
            call CFG_get(inputs, "captureGateClose", captureGateClose)
        
            ! Mathematical Inputs
            call CFG_get(inputs, "xPx", xPx)
            call CFG_get(inputs, "zPx", zPx)
            call CFG_get(inputs, "incidenceAngle", incidenceAngle)
            call CFG_get(inputs, "cosinePowerTD", cosinePowerTD)
            call CFG_get(inputs, "cosinePowerIS", cosinePowerIS)
            call CFG_get(inputs, "x0", x0)
            call CFG_get(inputs, "aMax", aMax)
            call CFG_get(inputs, "aMin", aMin)
            call CFG_get(inputs, "h", h)
            call CFG_get(inputs, "s", s)
            call CFG_get(inputs, "dist", dist)
            call CFG_get(inputs, "mass", mass)
            call CFG_get(inputs, "massMol", massMol)
            call CFG_get(inputs, "energyTrans", energyTrans)
            call CFG_get(inputs, "surfaceMass", surfaceMass)
            call CFG_get(inputs, "exitAngle", exitAngle)
            call CFG_get(inputs, "temp", temp)
            call CFG_get(inputs, "ncyc", ncyc)
            call CFG_get(inputs, "maxSpeed", maxSpeed)
            call CFG_get(inputs, "scatterFraction", scatterFraction)
        
            ! File paths go here
            call CFG_get(inputs, "linux", linux)
            call CFG_get(inputs, "imagePath", imagePath)
            call CFG_get(inputs, "matrixPath", matrixPath)
            call CFG_get(inputs, "ifPath", ifPath)
        
            ! Fixed parameters
        
            call CFG_get(inputs, "fixedIngoingSpeed", fixedIngoingSpeed)
            call CFG_get(inputs, "speedIn", speedIn)
        
            call CFG_get(inputs, "fixedOutgoingSpeed", fixedOutgoingSpeed)
            call CFG_get(inputs, "speedIn", speedIn)
        
            call CFG_get(inputs, "fixedStartPos", fixedStartPos)
            call CFG_get(inputs, "startx", startx)
            call CFG_get(inputs, "starty", starty)
            call CFG_get(inputs, "startz", startz)
        
            call CFG_get(inputs, "fixedScatterPos", fixedScatterPos)
            call CFG_get(inputs, "scatterx", scatterx)
            call CFG_get(inputs, "scattery", scattery)
            call CFG_get(inputs, "scatterz", scatterz)
            
            call CFG_get(inputs, "fixedCreationTime", fixedCreationTime)
            call CFG_get(inputs, "creationTime", creationTime)
        
            call CFG_get(inputs, "fixedScatterTime", fixedScatterTime)
            call CFG_get(inputs, "scatterTime", scatterTime)
        
            call CFG_get(inputs, "hush", hush)

            skimPos = skimPos + 0.08
            valvePos = valvePos  + 0.08
            colPos = colPos  + 0.08

        end subroutine load_inputs
        
end module inputs