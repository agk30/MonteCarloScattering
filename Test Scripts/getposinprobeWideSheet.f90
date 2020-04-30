program getposinprobeWideSheet
implicit none

!Dimenisions of the cuboid of the laser sheet
double precision :: LScentre, LSheight, LSwidth, pxmmratio
!Timing related variables
double precision :: t0, TprobeIn, TprobeF, Tstep, tWheel
double precision, dimension(:),  allocatable :: Tdelay, TdelaySet
!Variables related to the movement of the molecule
double precision, dimension(:) :: colposvalve(3), direction(3), velocity (3), colposprobe(3), scatteredDirection(3), colposwheel(3)
integer, dimension (:) :: colposprobePX(3)
!Variables related to the movement of the molecule taken from getspeeds.txt
double precision :: speed, scatteredSpeed
!Other variables used by the program
integer :: ncyc, NumberOfTimepoints, i, j, k, s, z, t, Xpixels, Zpixels, a, b, c, d, yAtFrontSheet, yAtBackSheet
!The 3D array representing the text image
integer, dimension(:,:,:),  allocatable :: Image
!Simulated image filename
character(20) :: filename
double precision :: yDirection, zDirection, bottomEdge, topEdge

!Open/Create and open the following files
open(unit=11,file='../getdirectionsv2.txt')    !info of ncyc, colposvalve(1), colposvalve(2), colposvalve(3), direction(1), direction(2), direction(3)
open(unit=12,file='../getspeeds.txt')		  !info of times of creation and speeds
open(unit=13,file='../getposinprobe.inp')    !Size and position of the laser sheet, time delays and pixel-to-mm ratio
open(unit=14,file='getposinprobe.txt')    !Output data from getposinprobe.f90
open(unit=15,file='../getdirections.inp')    !used for ncyc 
open(unit=16,file='../scatteredgetdirections.txt')	!scattered vectors and time of collision with wheel
open(unit=17,file='../getdirections.inp')
open(unit=18,file='../MBSpeeds.txt')

!Read the lines from the aforementioned files, reads from top to bottom, line by line and assigns names to the data from those files
read(13,*) LScentre
read(13,*) LSheight
read(13,*) LSwidth
read(13,*) !Line with the variation of the angle of incidence of the molecular beam
read(13,*) TprobeIn
read(13,*) TprobeF
read(13,*) Tstep
read(13,*) pxmmratio
read(11,*) !Reads first line where labels are so it doesn't get in the way later on
read(12,*) !Reads first line where labels are so it doesn't get in the way later on
read(16,*) !Reads first line where labels are so it doesn't get in the way later on
read(18,*) !Reads first line where labels are so it doesn't get in the way later on
read(17,*) ncyc

LSwidth = 0.1

!Calculate the number of timepoints of the probe
NumberOfTimepoints = ((TprobeF-TprobeIn)/Tstep) + 1
!Number of pixels in the simulated image of the probe volume in the x-axis
Xpixels = 420
!Number of pixels in the simulated image of the probe volume in the z-axis
Zpixels = 420

!Allocate the number of elemens to the following arrays
allocate(Tdelay(NumberOfTimepoints))
allocate(TdelaySet(NumberOfTimepoints))
allocate(Image(Zpixels,Xpixels,NumberOfTimepoints))

Image = 0.0d0

write(14,'(a)') "J, Tdelay / s, Position in probe (x, y, z) /m, Position in probe (x,z) /px"

!Choose a trajectory
do j = 1, ncyc
	
	!Read the trajectory's origin position and direction unit vectors from getdirections.txt
	read(11,'(i,6e20.6)') s, colposvalve(1), colposvalve(2), colposvalve(3), direction(1), direction(2), direction(3)
	!Read the trajectory's time of creation and speed from getspeeds.txt
	read(12,'(i,2e20.6)') z, t0, speed
	!Read the direction unit vectors and time between creation and collision with wheel from getscattereddirections.txt
	read(16,'(i,6e20.6)') i, scatteredDirection(1), scatteredDirection(2), scatteredDirection(3)
	!Read the Maxwell-Boltzmann distribution of speeds to assign to scattering molecules
	read(18,'(i,6e20.6)') k, scatteredSpeed

	do t = 1, NumberOfTimepoints

		!calculate time delay
		Tdelay(t) = t0 + TprobeIn + (t-1)*Tstep
		TdelaySet(t) = TprobeIn + (t-1)*Tstep

		!print *, 'Tdelay(t) is', Tdelay(t)

		!Calculate the velocity per axis
		velocity(1) = direction(1)*speed
		velocity(2) = direction(2)*speed
		velocity(3) = direction(3)*speed

		!print *, 'initial velocity is', velocity(3)
				
		!Calculate the final position of a trajectory
		colposprobe(1) = colposvalve(1)+(velocity(1)*Tdelay(t))
		colposprobe(2) = colposvalve(2)+(velocity(2)*Tdelay(t))
		colposprobe(3) = colposvalve(3)+(velocity(3)*Tdelay(t))

		yDirection = colposprobe(2)
		zDirection = colposprobe(3)
		bottomEdge = (LScentre - (LSwidth/2))
		topEdge = (LScentre + (LSwidth/2))
		
	
		! check if particle within laser sheet along z-axis
		if ((zDirection .ge. bottomEdge) .and. (zDirection .le. topEdge)) then
			! check for negative z direction if no surface present
			if (zDirection .gt. 0) then
				! check if particle within laser sheet along y-axis
				if ((yDirection .gt. -1*LSHeight) .and. (yDirection .le. LSheight)) then
						

					!calculates the pixel corresponding to the position in x and z axes
					colposprobePX(1) = ceiling(abs(((colposprobe(1)-0.036)/pxmmratio)))
					!0.036 refers to the width of the sheet + the distance between the sheet and the surface
					colposprobePX(3) = ceiling(abs(((colposprobe(3)-0.07)/pxmmratio)))

					!print *, 'colposprobePX(1) is', colposprobePX(1)

					if ((colposprobePX(1) .lt. Xpixels) .and. (colposprobePX(3) .lt. Zpixels)) then
						write(14,'(i7,4e20.6, 2i)') j, TdelaySet(t), colposprobe(1), colposprobe(2), colposprobe(3), colposprobePX(1), colposprobePX(3)

						Image(colposprobePX(3),colposprobePX(1),t) = Image(colposprobePX(3),colposprobePX(1),t) + 1
					end if
				end if
			end if
		end if
	
	end do
end do

do t =1, NumberOfTimepoints

	write(filename,"(A5,I2,A4)") "Image", t, ".txt"
	open(unit=20+t,file='Images/'//filename)

	do a = 1, Zpixels
	
		do b = 1, Xpixels
	
			write(20+t,'(i7)',advance='no') Image(a,b,t)

		end do
	
		write(20+t,*)
	
	end do
	
end do
	
end program getposinprobeWideSheet