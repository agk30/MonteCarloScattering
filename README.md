# MonteCarlo Scattering Simulation
Simulation of molecular beam scattering against a surface using Monte Carlo methods

## Concept
This program simulates a real experiment whereby a molecular beam is produced and aimed at a surface. The particles colliding with the surface will then scatter. The trajectories and speeds of the particles, both ingoing and scattered, are generated based on scientific models and real experimental data.

In the real experiment, the beam is imaged by passing through a laser sheet, allowing for detection of a particle that passes through. This simulation seeks to mimic this behaviour and will be used to compare scientific models to the real data obtained.

## Build instructions
Firstly, this program has been tested using the Intel Fortran Compiler 19.0 (XE 2019) on a Linux platform. This compiler is included in the Intel Parallel Studio (XE) package. Later versions of this compiler have not been tested. The program is also able to be compiled by gfortran in the normal way, but the `-fdec-math` is required. In Windows environments, ifort compiles the program, but the program fails to run currently, workingn towards a fix. Gfortran compiles the program properly and produces expected results. Once you have access to these compilers, clone or download this repository. Next, `cd` to the appropriate directory, ensuring the file MCScattering.f90 is contained within that directory. Simply compile the program using the command `ifort MCScattering.f90` or `gfortran -fdec-math MCScattering.f90` with any appropriate flags if you wish to rename the executable for example, and the .exe or a.out (or whatever you have named it) should be generated. No additional flags are necessary to compile this program.

## Inputs and Imaging
The inputs of the program are contained within the `inputs.inp` file and may be altered at will to suit your requirements. The default inputs have been tested to work correctly. See the file itself for definition of each parameter.

**Important:** Folders named `Images`, `Images 2` and `Images 3` must be created in the same directory as MCScattering.exe before execution for images to be obtained until a new version that can create the folder is released.

Once input parameters have been selected and the program is run, there will be a sequence of images generated in the `Images` folder. These contain matrices of integers corresponding to intensity at each given pixel region. To view these images and to obtain a video of the sequence, the Image J program should be used. The program can be found [here](https://imagej.nih.gov/ij/download.html). Once the program is downloaded and executed, it you must direct it to the imaging macro contained within the `Imaging Macros` folder of this repository. In Image J, navigate to Plugins > Macros > Install, then direct the program to the imaging macro. Once this is done, you can now load the images into the program to view as a video. Now, navigate to Plugins > Macros > ImportSeriesTextImages. This shoud bring up a window showing the video of the images in succession. From here, I recommend going to Image > Lookup Tables then slecting Royal, but any of these would be fine to suit your needs. To save as a video, navigate to File > Save As > AVI... where you can save the video at a given framerate. Now you should have the desired output of this program.

## Contact

If you require any assistance building or running this program, contact me at agk30@hw.ac.uk
