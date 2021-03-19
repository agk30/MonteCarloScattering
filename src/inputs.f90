module inputs
    use mathConstants
    use m_config

    contains

        ! Loads input parameters into the main section of code, MCScattering.f90
        ! See inputs.inp for details on parameter definitions
        ! TODO pass over hash table instead of individual variables
         subroutine load_inputs (xPx, zPx, incidenceAngle, ncyc, x0, aMax, aMin, &
             h, s, dist, pulseLength, mass, massMol, energyTrans, surfaceMass, exitAngle, temp, skimPos, valvePos, colPos, &
             skimRad, valveRad, colRad, sheetCentre, halfSheetHeight, sheetWidth,&
              probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, scattering, gaussDev, ksize, polyOrder, testMods,&
               writeImages, fullSim, scatterFraction, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose, &
                cosinePowerTD, cosinePowerIS, runNumber, imagePath, blurredImagePath, ifImagePath, matrixPath, ifPath)
            implicit none

            integer, intent(out) :: ncyc, xPx, zPx, ksize, polyOrder, cosinePowerTD, cosinePowerIS, runNumber
            double precision, intent(out) :: incidenceAngle, x0, aMax, aMin, &
            h, s, dist, pulseLength, mass, temp, valvePos, gaussDev, massMol, energyTrans, surfaceMass, exitAngle
            double precision, intent(out) :: skimPos, colPos, skimRad, valveRad, colRad, sheetCentre, &
             halfSheetHeight, sheetWidth, probeStart, probeEnd, tStep, pxMmRatio, maxSpeed, &
              scatterFraction, scatterIntensity, fLifeTime, captureGateOpen, captureGateClose
            logical, intent(out) :: scattering, testMods, writeImages, fullSim
            character(200) :: imagePath, blurredImagePath, ifImagePath, matrixPath, ifPath

            type(CFG_t) :: my_cfg

            ! Build the variables used in the input file

            !Inputs for bash script running
            call CFG_add(my_cfg, "runNumber", 1 , "Position in run sequence")

            !Experimental inputs
            call CFG_add(my_cfg, "skimPos", 0.1730D0 , "Position of Skimmer in z direction")
            call CFG_add(my_cfg, "valvePos", 0.2150D0 , "Position of Valve in z direction")
            call CFG_add(my_cfg, "colPos", 0.1330D0 , "Position of Collimator in z direction")
            call CFG_add(my_cfg, "skimRad", 0.001D0 , "Skimmer Radius")
            call CFG_add(my_cfg, "valveRad", 0.0015D0 , "Valve Radius")
            call CFG_add(my_cfg, "colRad", 0.0015D0 , "Collimator Radius")
            call CFG_add(my_cfg, "sheetCentre", 0.021D0 , "Distance from centre of sheet from to surface")
            call CFG_add(my_cfg, "halfSheetHeight", 0.002D0 , "Half height of laser sheet")
            call CFG_add(my_cfg, "sheetWidth", 0.0300D0 , "Full width of laser sheet")
            call CFG_add(my_cfg, "pulseLength", 10.0D-06 , "Discharge pulse length")
            
            ! Imaging inputs
            call CFG_add(my_cfg, "pxMmRatio", 0.25D-03 , "Pixel to mm ratio for imaging")
            call CFG_add(my_cfg, "probeStart", 67.0D-06 , "Start of probe time (imaging time)")
            call CFG_add(my_cfg, "probeEnd", 180.0D-06 , "End of probe time (imaging time)")
            call CFG_add(my_cfg, "tStep", 1.0D-06 , "Time step between images")
            call CFG_add(my_cfg, "gaussDev", 2.0D0 , "Deviation parameter for gaussian blur routine")
            call CFG_add(my_cfg, "ksize", 27 , "Kernel size for SG routine")
            call CFG_add(my_cfg, "polyOrder", 3 , "Polynomial order for SG routine")
            call CFG_add(my_cfg, "scattering", .TRUE. , "Image scattering as well as ingoing beam?")
            call CFG_add(my_cfg, "fullSim", .TRUE. , "Image ingoing beam as well as scattering?")
            call CFG_add(my_cfg, "testMods", .FALSE. , "Including testing modules?")
            call CFG_add(my_cfg, "writeImages", .TRUE. , "Wirte images to files?")
            call CFG_add(my_cfg, "scatterIntensity", 3.0D0 , "Relative intensity of scattered signal to ingogin signal")
            call CFG_add(my_cfg, "fLifeTime", 700D-9 , "Fluorescent lifetime of the molecule")
            call CFG_add(my_cfg, "captureGateOpen", 350D-9 , "Time after probe fire at which imaging starts")
            call CFG_add(my_cfg, "captureGateClose", 10000D-9 , "Time after probe fire at which imaging ends")
            
            ! Mathematical Inputs
            call CFG_add(my_cfg, "xPx", 420 , "Number of image pixels in x direction")
            call CFG_add(my_cfg, "zPx", 420 , "Number of image pixels in z direction")
            call CFG_add(my_cfg, "incidenceAngle", 0D0 , "Incidence angle for beam")
            call CFG_add(my_cfg, "cosinePowerTD" , 1 , "Power of cosine weighted TD angle distribution (cos^n)")
            call CFG_add(my_cfg, "cosinePowerIS" , 4 , "Power of cosine weighted IS angle distribution (cos^n)")
            call CFG_add(my_cfg, "x0", 128.30344D0 , "Origin function parameter for speed generation")
            call CFG_add(my_cfg, "aMax", 1.00082D0 , "Origin function parameter for speed generation")
            call CFG_add(my_cfg, "aMin", -0.00851D0 , "Origin function parameter for speed generation")
            call CFG_add(my_cfg, "h", 24.98601D0 , "Origin function parameter for speed generation")
            call CFG_add(my_cfg, "s", 1.45285D0 , "Origin function parameter for speed generation")
            call CFG_add(my_cfg, "dist", 230.0D-03 , "Single-point LIF value-probe laser distance")
            call CFG_add(my_cfg, "mass", 2.8240519D-26 , "Molecule mass / kg")
            call CFG_add(my_cfg, "massMol", 0.017D0 , "Molecue mass / kg/mol")
            call CFG_add(my_cfg, "energyTrans", 0.0D0 , &
            "SS collision model - energy loss to internal surface motions")
            call CFG_add(my_cfg, "surfaceMass", 100D0 , "Effective liquid surface mass / g/mol")
            call CFG_add(my_cfg, "exitAngle", 0D0 , "Exit angle for SS model")
            call CFG_add(my_cfg, "temp", 298.0D0 , "Surface temp / K")
            call CFG_add(my_cfg, "ncyc", 10000000 , "Number of molcules to be sampled")
            call CFG_add(my_cfg, "maxSpeed", 3000.0D0 , "Max speed for MB speed calculation")
            call CFG_add(my_cfg, "scatterFraction", 0.5D0  , &
            "Fraction of molcules scattering in TD or IS. 0 for full TD, 1 for full IS")

            ! File Paths go here
            call CFG_add(my_cfg, "imagePath", "../Images/Images/" , "Path to image files")
            call CFG_add(my_cfg, "blurredImagePath", "../Images/Blurred Images/" ,&
             "Path to Gaussian blurred image files")
            call CFG_add(my_cfg, "ifImagePath", "../Images/IF Adjusted Images/",&
             "Path to IF adjusted image files")
            call CFG_add(my_cfg, "matrixPath", "../SG Matrices/CC_027x027_003x003.dat",&
             "Path to SG matrix")
            call CFG_add(my_cfg, "ifPath", "../Real Images/09032021_Q12_Instrument Function_SUM_00",&
             "Path to instrument function image")
            
            ! Read .cfg file and parse for changes
            call CFG_read_file(my_cfg, "../Inputs/inputs.cfg")

            ! Parse any command line arguments
            call CFG_update_from_arguments(my_cfg)

            ! Assign variables from .cfg
            call CFG_get(my_cfg, "runNumber", runNumber)

            call CFG_get(my_cfg, "skimPos", skimPos)
            call CFG_get(my_cfg, "valvePos", valvePos)
            call CFG_get(my_cfg, "colPos", colPos)
            call CFG_get(my_cfg, "skimRad", skimRad)
            call CFG_get(my_cfg, "valveRad", valveRad)
            call CFG_get(my_cfg, "colRad", colRad)
            call CFG_get(my_cfg, "sheetCentre", sheetCentre)
            call CFG_get(my_cfg, "halfSheetHeight", halfSheetHeight)
            call CFG_get(my_cfg, "sheetWidth", sheetWidth)
            call CFG_get(my_cfg, "pulseLength", pulseLength)

            call CFG_get(my_cfg, "pxMmRatio", pxMmRatio)
            call CFG_get(my_cfg, "probeStart", probeStart)
            call CFG_get(my_cfg, "probeEnd", probeEnd)
            call CFG_get(my_cfg, "tStep", tStep)
            call CFG_get(my_cfg, "gaussDev", gaussDev)
            call CFG_get(my_cfg, "ksize", ksize)
            call CFG_get(my_cfg, "polyOrder", polyOrder)
            call CFG_get(my_cfg, "scattering", scattering)
            call CFG_get(my_cfg, "fullSim", fullSim)
            call CFG_get(my_cfg, "testMods", testMods)
            call CFG_get(my_cfg, "writeImages", writeImages)
            call CFG_get(my_cfg, "scatterIntensity", scatterIntensity)
            call CFG_get(my_cfg, "fLifeTime", fLifeTime)
            call CFG_get(my_cfg, "captureGateOpen", captureGateOpen)
            call CFG_get(my_cfg, "captureGateClose", captureGateClose)

            call CFG_get(my_cfg, "xPx", xPx)
            call CFG_get(my_cfg, "zPx", zPx)
            call CFG_get(my_cfg, "incidenceAngle", incidenceAngle)
            call CFG_get(my_cfg, "cosinePowerTD", cosinePowerTD)
            call CFG_get(my_cfg, "cosinePowerIS", cosinePowerIS)
            call CFG_get(my_cfg, "x0", x0)
            call CFG_get(my_cfg, "aMax", aMax)
            call CFG_get(my_cfg, "aMin", aMin)
            call CFG_get(my_cfg, "h", h)
            call CFG_get(my_cfg, "s", s)
            call CFG_get(my_cfg, "dist", dist)
            call CFG_get(my_cfg, "mass", mass)
            call CFG_get(my_cfg, "massMol", massMol)
            call CFG_get(my_cfg, "energyTrans", energyTrans)
            call CFG_get(my_cfg, "surfaceMass", surfaceMass)
            call CFG_get(my_cfg, "exitAngle", exitAngle)
            call CFG_get(my_cfg, "temp", temp)
            call CFG_get(my_cfg, "ncyc", ncyc)
            call CFG_get(my_cfg, "maxSpeed", maxSpeed)
            call CFG_get(my_cfg, "scatterFraction", scatterFraction)

            ! File paths go here
            call CFG_get(my_cfg, "imagesPath", imagePath)
            call CFG_get(my_cfg, "blurredImagesPath", blurredImagePath)
            call CFG_get(my_cfg, "ifImagesPath", ifImagePath)
            call CFG_get(my_cfg, "matrixPath", matrixPath)
            call CFG_get(my_cfg, "ifPath", ifPath)

        end subroutine load_inputs
        
end module inputs