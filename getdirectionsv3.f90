! A fortran95 program for G95
! By WQY

include "Maths/getscattereddirections.f90"
include "Maths/directionMath.f90"

program getdirections

	use class_getscattereddirections
	use class_directions

	!Calculates a set of directions for the collider intersecting the wheel surface from a specific experimental geometry. Saves the results as formatted file for reading in by the program getvelocities.
	implicit none
	double precision :: valvepos,skimpos,collimatorpos,rvalve,rskim,rcollimator !Details of position and radius of valve and skimmer orifices in the chamber.
	Double precision, dimension(:) :: colposskim(3), colposvalve(3), direction(3), scattereddirection(3), eulercolposvalve(3),eulerdirection(3) !XYZ positions of the collider at the skimmer, valve and direction of the trajectory. x=along heriot cell, y=perpendicular to heriot cell, z=direction of molecular beam
	double precision, dimension(:) :: colposcollimator(3)
	Double precision :: mx, cx, my,cy,cz, a, z !variables to describe the lines which connect the point on the skimmer and valve. a is used to calculate the unit vectors see subroutine called unitvector.
	Logical :: hit, Hdata,scattered !allows the molecular beam to be weighted by cos4theta, should it save heriot cell data?
	Double precision :: theta !angle of velocity the z axis which is used to weight the trajectories by cos4theta
	integer :: j, i, ncyc,MoleBeamAng,totaltrajectories ! indexing and number of cycles and molecular beam angle to the surface normal
	Double precision :: xHC, yHC, HCpos, rHC !position of trajectories at the position of the heriot cell so one can determine the size of the molecular beam at the heriot cell.
	
	open(unit=11,file='heriotcellpositions.txt') !file writing trajectory positions at a given distance from the skimmer
	open(unit=12,file='getdirectionsv2.txt') !writes initial position on valve and unit vectors to a txt file for read into program getvelocities
	open(unit=13,file='getdirections.inp') !input file which contained the geometries of the chamber ncyc and a position to output trajectory positions
	open(unit=14,file='rHC.txt') !output radius of the beam at the heriot cell
	open(unit=15,file='scatteredgetdirections.txt') !output radius of the beam at the heriot cell

	!read in chamber geometries and number of cycles 
	read(13,*) ncyc
	read(13,*) scattered
	read(13,*) HCpos
	read(13,*) skimpos
	read(13,*) valvepos
	read(13,*) collimatorpos
	read(13,*) rskim
	read(13,*) rvalve
	read(13,*) rcollimator
	read(13,*) Hdata
	read(13,*) MoleBeamAng

	Write(12,'(7a)') "J,xvalve,yvalve,zvalve,vx,vy,vz"
	Write(15,'(7a)') "J,vx, vy ,vz"

	totaltrajectories = 0
		!Ramdomly choose xy coordinates of two points in the skimmer and valve orifices for which the collider will pass through. The zcoodinate comes from positions of valve and skimmer stored in getdirection.inp 
		do j = 1, ncyc
	
			hit = .false.! makes sure we start from not hit so that the angle of the trajectoy can be weighted by cos4theta later in the code
	
			!Makes sure that the cos4theta density of the molecular beam is accounted for by the If statement in the subsoutiene find angle
			do while (.not. hit) 
				call discpick(colposcollimator(1),colposcollimator(2))
				colposcollimator=colposcollimator*rcollimator
				colposcollimator(3)=collimatorpos

				call discpick(colposvalve(1),colposvalve(2))
				colposvalve=colposvalve*rvalve
				colposvalve(3)=valvepos

				!calculates the gradient and intercept of the line which connects the two points in the skimmer and valve
				call fitline(mx,cx,colposvalve(1),colposvalve(3),colposcollimator(1),colposcollimator(3))
				call fitline(my,cy,colposvalve(2),colposvalve(3),colposcollimator(2),colposcollimator(3))
				cz=0

				!Calculates the position of the collider at the distance of the collimator
				!colposcollimator(1) = mx*collimatorpos + cx
				!colposcollimator(2) = my*collimatorpos + cy
				!colposcollimator(3) = collimatorpos
				!z = sqrt(colposcollimator(1)**2 + colposcollimator(2)**2)

				totaltrajectories = totaltrajectories+1

				!checks whether the collider passes through the collimator or not
				!if (z .le. rcollimator) then

				!Convert lines into unit vectors for the directions
				Call unitvector(mx,my,direction(1),direction(2),direction(3),a)
		
				call findangle(mx,my,theta,hit,colposvalve(3),colposcollimator(3))

				!end if

			end do
			
			!If molecular beam is not normal to the surface then we have to perform an euler transformation 
			if (MoleBeamAng.NE.0) then	
				!second angle is 0 because you are only roating about the y axis as I do for the laser beams.
				call euler(ncyc,MoleBeamAng,0, eulercolposvalve,eulerdirection, colposvalve, direction)	
				colposvalve=eulercolposvalve
				direction=eulerdirection
				!Need to recalcualte where molecular beam now intercts the wheel for later
				call fitline(mx,cx,colposvalve(1),colposvalve(3),colposskim(1),colposskim(3))
				call fitline(my,cy,colposvalve(2),colposvalve(3),colposskim(2),colposskim(3))
			end if
			
			write(12,'(i7,6e20.8)') j ,colposvalve(1) ,colposvalve(2), colposvalve(3), direction(1), direction(2),direction(3)
		
			if (Hdata) then
				! calcualte the size of beam at heriot cell
				xHC=mx*HCpos+cx
				yHC=my*HCpos+cy
				write(11,'(i7,2e20.8)') j, xHC, yHC
			
				if (xHC .gt. rHC) then
					rHC = xHC
				end if
			
				if (j==ncyc) then
					write(14,*) rHC, "!radius of beam at Heriot Cell"
				end if 
			end if

			if (scattered) then

				scattereddirection = getscattereddirections() !generate scattering direction
	
				write(15,'(i7,6e20.8)') j,scattereddirection(1), scattereddirection(2), scattereddirection(3) !writes intersect with wheel to file
				
			end if

		end do

		!For each trajectory if the wheel is in then assign directions for each trajectory at random with a costheta distribution
	
	
	!Deletes beam info at heriot cell if not required see getdirections.inp to change
	if (.not.Hdata) then
	CLOSE(unit=11,status='delete')
	CLOSE(unit=14,status='delete')
	end if
	
	if (.not.scattered) then
			CLOSE(unit=16,status='delete')
	end if	
print *, 'total trajectories = ',totaltrajectories
end program getdirections
