# MonteCarlo Scattering Simulation
Simulation of molecular beam scattering against a surface using Monte Carlo methods

## Concept
This program simulates a real experiment whereby a molecular beam is produced and aimed at a surface. The particles colliding with the surface will then scatter. The trajectories and speeds of the particles, both ingoing and scattered, are generated based on scientific models and real experimental data.

In the real experiment, the beam is imaged by passing through a laser sheet, allowing for detection of a particle that passes through. This simulation seeks to mimic this behaviour and will be used to compare scientific models to the real data obtained.

## Build instructions
The program may be built using the gfortran compiler and it is the recommended compiler to ensure functionality across different environments. Other compilers have not been validated, although ifort is likely to work. Once you have access to these compilers, clone or download this repository. In the repository's root directory, build the program:

`mkdir build`

`cmake ./ -B ./build/ -G "Unix Makefiles"`

`cd build`

`make`

Then to run, execute `MCScattering`

`./MCScattering`

## Inputs and Imaging
The program contains a default set of variables which may be overwritten using either the `inputs.cfg` file or through command line arguments (.cfg file overrides defaults, commmand line arguemnts overrides both defaults and .cfg file). Default input values and their descriptions are as follows:

**Batch Processing**
- `runNumber = 1` Index used for numbering runs. Typically this will be passed as a command line argument and will create different folders depending on the run number given.

**Experimental Apparatus Inputs**
- `skimPos = 0.1730D0` Position of Skimmer in z direction / m
- `valvePos = 0.2150D0` Position of Valve in z direction / m
- `colPos = 0.1330D0` Position of Collimator in z direction / m
- `skimRad = 0.001D0` Skimmer Radius / m
- `valveRad = 0.0015D0` Valve Radius / m
- `colRad = 0.0015D0` Collimator Radius / m
- `sheetCentre = 0.021D0` Distance from centre of sheet from to surface / m
- `halfSheetHeight = 0.002D0` Half height of laser sheet / m
- `sheetWidth = 0.0300D0` Full width of laser sheet / m
- `pulseLength = 10.0D-06` Discharge pulse length / s

**Imaging Inputs**
- `pxMmRatio = 0.25D-03` Pixel to mm ratio for imaging
- `probeStart = 67.0D-06` Start of probe time (imaging time) / s
- `probeEnd = 180.0D-06` End of probe time (imaging time) / s
- `tStep = 1.0D-06` Time step between images / s
- `gaussDev = 2.0D0` Deviation parameter for gaussian blur routine
- `ksize = 27` Kernel size for SG routine
- `polyOrder = 3` Polynomial order for SG routine
- `scattering = .TRUE.` Image scattering as well as ingoing beam?
- `fullSim = .TRUE.` Image ingoing beam as well as scattering?
- `testMods = .FALSE.` Including testing modules?
- `writeImages = .TRUE.` Write images to files?
- `scatterIntensity = 3.0D0` Relative intensity of scattered signal to ingoing signal
- `fLifeTime = 700D-9` Fluorescence lifetime of the molecule / s
- `captureGateOpen = 350D-9` Time after probe fire at which emission can be collected
- `captureGateClose = 360D-9` Time after probe fire at which emission will stop being collected

**Mathematical Inputs**
- `xPx = 420` Number of image pixels in x direction
- `zPx = 420` Number of image pixels in z direction
- `incidenceAngle = 0D0` Incidence angle for beam
- `x0 = 128.30344D0` Origin function parameter for speed generation
- `aMax = 1.00082D0` Origin function parameter for speed generation
- `aMin = -0.00851D0` Origin function parameter for speed generation
- `h = 24.98601D0` Origin function parameter for speed generation
- `s = 1.45285D0` Origin function parameter for speed generation
- `dist = 230.0D-03` Single-point LIF value-probe laser distance
- `mass = 2.8240519D-26` Molecule mass / kg
- `massMol = 0.017D0` Molecule mass / kg/mol
- `energyTrans = 0.0D0 ` SS collision model - energy loss to internal surface motions
- `surfaceMass = 100D0` Effective liquid surface mass / g/mol
- `exitAngle = 0D0` Exit angle for SS model
- `temp = 298.0D0` Surface temp / K
- `ncyc = 10000000` Number of molecules to be sampled
- `maxSpeed = 3000.0D0` Max speed for MB speed calculation
- `scatterFraction = 0.5D0` Fraction of molecules scattering in TD or IS. 0 for full TD, 1 for full IS

**File Paths**

File paths should be given as the directories to which images will be written. An additional file containing input paramters used in the run will be written to the subdirectory under this path as `/Run 1/input_values.cfg` for example.

- `imagePath` Path to output images
- `matrixPath` Path to SG matrix
- `ifPath` Path to input instrument function image

Command line arguments may be passed in the form of:

`-variableName='value'`

Once input parameters have been selected and the program is run, there will be a sequence of images generated in the folders specified by input paramters. These contain matrices of numbers corresponding to intensity at each given pixel region. Three different directories will be created which will contain contains raw images, Gaussian blurred images and the blurred images with the addition of the Savitzky-Golay filter. To view these images and to obtain a video of the sequence, the Image J program should be used. The program can be found [here](https://imagej.nih.gov/ij/download.html). Once the program is downloaded and executed, it you must direct it to the imaging macro, `ImportSeriesTextImages.txt`, contained within the `Imaging Macros` folder of this repository. In Image J, navigate to Plugins > Macros > Install, then direct the program to the imaging macro. Once this is done, you can now load the images into the program to view as a video. Now, navigate to Plugins > Macros > ImportSeriesTextImages. This should bring up a window showing the video of the images in succession. From here, I recommend going to Image > Lookup Tables then slecting Royal, but any of these would be fine to suit your needs. To save as a video, navigate to File > Save As > AVI... where you can save the video at a given framerate. Now you should have the desired output of this program.

## Contact

If you require any assistance building or running this program, contact me at agk30@hw.ac.uk

