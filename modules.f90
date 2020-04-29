subroutine scatteringangle(xscatterdirection,yscatterdirection,zscatterdirection)
implicit none

logical :: hit
double precision :: x4,x5,x6,x7,x8
double precision :: costheta, theta, phi,pi,directionxy
double precision, intent(out) :: xscatterdirection,yscatterdirection,zscatterdirection

pi=3.14159265359

	!randomly assigns a dirrection for each of these trajectories to scatter at with a costheta distribution
			hit=.false.
		
		!Pick theta for angle for leavning the surface and weight this by costheta
		Do while (.not. hit)
			! randomly selec thetax
			call random_number(x4)
			theta=x4*pi/2
			costheta= (cos(theta))

			call random_number(x5)
			
		
			!Hit always starts as false due to program getdirections. Here it is changed to true if trajectory has cos4theta > random number. The trajectory is accepted. If not the loop runs again.
			If (costheta .gt. x5) then
			hit = .true.
			end if 
		
		end do
		
		!Calculate scattering direction 
		zscatterdirection=1*sin(theta)
		directionxy=1*cos(theta)
		call random_number(x6)
		phi=x6*pi/2
		
		call random_number(x7)
		if (x7.gt.0.5) then
		xscatterdirection=directionxy*sin(phi)
		else
		xscatterdirection=-directionxy*sin(phi)
		end if
		
		call random_number(x8)
		if (x8.gt.0.5) then
		yscatterdirection=directionxy*cos(phi)
		else
		yscatterdirection=-directionxy*cos(phi)
		end if

end subroutine scatteringangle

subroutine discpick(x,y)
!Randomly pick a points from a unit radius circle.
        implicit none

		double precision, intent(out) :: x,y !x and y positions on unit radius circle
		double precision :: rand1,rand2 !two random numbers called for the distance of point from centre of circle and angle of this point from the centre

		!Random number for distance point is from centre of unit circle the value is square routed so that the points alone this line will result in an even distribution of point on the unit circle rather than at a higher density at the centre
		call random_number(rand1)
		rand1 = sqrt(rand1)
		
		!Random number called for the angle this point is on the circle. The number is then converted into an angle in radians.
		call random_number(rand2)
		rand2 = 2*3.14159265359*rand2
		
		!trigonometry to turn the distance from centre of circle and angle into x and y coordinates on the unit circle
		x = rand1*cos(rand2)
		y = rand1*sin(rand2)

end subroutine discpick

subroutine fitline(m,c,y2,x2,y1,x1) 
!calculates the gradient and intercept of the line which connect the point on the skimmer and valve
!note x and y are used here to mean y=mx+c not in reference to chamber coordinates (z chamber is x here)
        implicit none

		double precision, intent(out) :: m,c ! gradient and intercept of line which describes trajectory
		double precision, intent(in) :: y2,x2,y1,x1 !x and y coordinates of points on valve and skimmer
		
		!calculates the gradient and intercept of straight line which describes the trajectory
		m=(y2-y1)/(x2-x1)
		c=y2-(m*x2)

end subroutine fitline
	
subroutine unitvector(mx,my,vx,vy,vz,a)
!Calcualtes unit vectors from gradients of the lines in the x and y directions. 

	double precision, intent(in) :: mx, my
	double precision, intent(out) :: vx, vy,vz, a
	
	a=sqrt(mx**2+my**2+1**2)
	vx = mx / a
	vy = my / a 
	
	!If statement makes sure that the vz is the correct sign
	if (z2.gt.z1) then
	vz = 1/a 
	else
	vz = -1/a 
	end if
	
end subroutine unitvector

subroutine findangle(mx,my,theta,hit,z1,z2) 
!calculates the angle of the trajectory from the molecular beam axis then accepts the trajectory if cos4theta is larger than a random number so the molecular beam has the form cos4theta. 
        implicit none
		
		double precision, intent(in) :: mx,my,z1,z2 !gradients/vector of trajectory in the x and y directions
		Double precision :: mxy, rand3,cos4theta !gradient/vector of trajectory about the z axis, random number to compare theta to, cos4theta
		Logical , intent(out) :: hit !set to true if trajectory is accepted and false is rejected
		Double precision, intent(out) :: theta !angle of trajectory from z axis
		
		!calculates the vector of trajectory about the z axis and then angle of vector from z axis then cos4theta 
		mxy = sqrt((mx**2)+(my**2))
		theta=atan(mxy/1)
		cos4theta= (cos(theta))**4

		call random_number(rand3)
		
		!Hit always starts as false due to program getdirections. Here it is changed to true if trajectory has cos4theta > random number. The trajectory is accepted. If not the loop runs again.
		If (cos4theta .gt. rand3) then
		hit = .true.
		end if 
		
	end subroutine findangle

		
subroutine rnm(rmean,sd, rand)
!Produces normal deviates.	
		double precision :: random,rmean,sd,rand
		real*4 :: x,u,s,t,c0,c1,c2,d1,d2,d3,t2,t3!,rand,rmean,sd	

 	
		call random_number(x)
	  
		if(x.gt.0.999999)then
			x=0.999999
		endif
		
		if(x.lt.1.0e-6)then
			u=1.0e-06
		else
			u=x
		endif
		s = 1.
	       
		if (u .gt. 0.5) then
			u = 1. - u
			s = -1.
		end if
      
		t = sqrt(- (2 * log(u)))
		c0 = 2.515517
		c1 = 0.802853
		c2 = 0.010328
		d1 = 1.432788
		d2 = 0.189269
		d3 = 0.001308
		t2 = t * t
		t3 = t * t2
		x = s * (t - (((c0 + (c1 * t)) + (c2 * t2)) / (((1. + (d1 * t)) + &
	(d2 * t2)) + (d3 * t3))))
	
		rand = rmean + (dble(x) * sd)
		return 
	end subroutine rnm
	

 subroutine euler(ncyc,theta1,theta2,valvepos1, velocity1, valvepos, velocity)
		implicit none
		double precision, intent(in) :: theta1,theta2 !theta 1 is y in chamber axis coordinates
		Integer :: i, j
		Double precision:: t0, ctheta1, ctheta2, stheta1, stheta2
		Double precision, dimension(:), intent(out) :: valvepos1(3), velocity1(3)
		Double precision, dimension(:), intent(in) :: valvepos(3), velocity(3)
		Integer, intent(in):: ncyc
	
		double precision :: cphi,ctheta,cchi,sphi,stheta,schi
	
	!Calculate cos and sin of theta1 and theta2

		ctheta1 = cos(theta1)
		ctheta2 = cos(theta2)
		stheta1 = sin(theta1)
		stheta2 = sin(theta2)
		
		
	!rotates about the z axis
	velocity1(1)=ctheta1*ctheta2*velocity(1)+stheta1*velocity(2)+ctheta1*stheta2*velocity(3)
	velocity1(2)=-stheta1*ctheta2*velocity(1)+ctheta1*velocity(2)-stheta1*stheta2*velocity(3)
	velocity1(3)=-stheta2*velocity(1)+ctheta2*velocity(3)
	
	!Rotates initial position
	valvepos1(1)=ctheta1*ctheta2*valvepos(1)+stheta1*valvepos(2)+ctheta1*stheta2*valvepos(3)
	valvepos1(2)=-stheta1*ctheta2*valvepos(1)+ctheta1*valvepos(2)-stheta1*stheta2*valvepos(3)
	valvepos1(3)=-stheta2*valvepos(1)+ctheta2*valvepos(3)
	
end subroutine euler

subroutine beampositioncalculation(npasses,xpass, ypass, xtheta, ytheta)
implicit none

Double precision, dimension(:),  allocatable :: mirrorspotx,mirrorspoty
integer :: nspots, j
integer, intent(in) :: npasses
double precision, dimension(:),  intent(out):: xpass(npasses), ypass(npasses), xtheta(npasses), ytheta(npasses)

open(unit=20, file='mirrorspotposx.inp')
open(unit=21, file='mirrorspotposy.inp')
open(unit=22, file='mirrorspotpos.txt')
open(unit=23, file='beampositioncalculation.txt')

nspots=npasses+1
allocate(mirrorspotx(nspots))
allocate(mirrorspoty(nspots))


!First read in spot positions on the two mirrors
do j=1, nspots 
read(20,*) mirrorspotx(j)
read(21,*) mirrorspoty(j)
write(22,'(2e20.6)') mirrorspotx(j), mirrorspoty(j)
end do

write(23,*) 'xpass(j)/mm, ypass(j)/mm, xtheta(j)/rad, ytheta(j)/rad'
! calcualtes the crossing point of laser on z axis (i.e. centre of molecular beam) and the two angles for Euler transformation and writes them to file
do j=1, npasses
xpass(j)=(((mirrorspotx(j+1)-mirrorspotx(j))/400)*200+mirrorspotx(j))
ypass(j)=((mirrorspoty(j+1)-mirrorspoty(j))/400)*200+mirrorspoty(j)
xtheta(j)= TANH((mirrorspotx(j+1)-mirrorspotx(j))/400)
ytheta(j)= TANH((mirrorspoty(j+1)-mirrorspoty(j))/400)


write(23,'(4e20.6)') xpass(j), ypass(j), xtheta(j), ytheta(j)
end do

close(23, Status='KEEP')

end subroutine beampositioncalculation

subroutine circlelineintersect(intersect,initialposition, velocity,posentry, posexit, poscentrecylinder, r,option)
!inputs are initials positions and velocity of trajectory along with laser beam position and outputs are if there was an intersection and what the entry and exit points were.
	implicit none
	double precision ::x1,y1,x2,y2, a,b,c !coefficients of quadratic equation
	double precision :: x3,y3 !x and y coordinates for line and circle intercept
	Double precision :: u1,u2 !routes for entry and exit points of trajectory in laser beam
	Logical :: intersect
	Double precision, dimension(:), intent(in) :: initialposition(3),velocity(3) ! initial position and velocity of trajectory
	Double precision, dimension(:), intent(in) :: poscentrecylinder(3) ! position centre of cylinder
	Double precision, intent(in) :: r ! radius of cyclinder
	Double precision, dimension(:), intent(out) :: posentry(3), posexit(3) ! positions entry and exit 
	Integer :: option !0 if circle is in the YZ plane (FM set up) and 1 if in the XZ plane(LIF set up) inout from detectbeam.inp
	
	!If statement places circle in the YZ plane (FM set up) and 1 if in the XZ plane(LIF set up)
	If (option.eq.0) then
	!Calcuate x1,y1,x2,y2 from the initial skimmer positions and velocity of the trajectory 
		x1=initialposition(2)
		y1=initialposition(3)
		x2=x1+velocity(2)
		y2=y1+velocity(3)
		x3= poscentrecylinder(2)
		y3= poscentrecylinder(3)
	elseif 	(option.eq.1) then
		x1=initialposition(1)
		y1=initialposition(3)
		x2=x1+velocity(1)
		y2=y1+velocity(3)
		x3= poscentrecylinder(1)
		y3= poscentrecylinder(3)
	else 
		write(*,*) 'error incorrect input for circle line intercept option should be 0 or 1 but is =', option 
	end if
	
	!Calcualtes coefficients of the quadratic equation which descrbives the intersect of a line and a circle
		a=((x2-x1)**2)+((y2-y1)**2)
		b=2*((x2-x1)*(x1-x3)+(y2-y1)*(y1-y3))
		c=(x3**2)+(y3**2)+(x1**2)+(y1**2)-2*(x3*x1+y3*y1)-(r**2)
		
		!Does the trajectory intersect the beam? If it does calcualte the routes
		If ((b**2)>(4*a*c)) then
		intersect= .true.
		u1=(-b+(sqrt(b**2-4*a*c)))/(2*a)
		u2=(-b-(sqrt(b**2-4*a*c)))/(2*a)
		
		!Calculates the entry and exit positions of the trajectory 
		
		posentry =initialposition+u1*velocity
		posexit=initialposition+u2*velocity
		else
		intersect=.false.
		end if
		
		!check that entry and exit positions the right way round
		!add this
		

		
	end subroutine circlelineintersect
	
subroutine TOFcalc(posintersect, t0, velocity, initialposition,tarr,TOF)
!Calcualtes the time of arrivals of each trajectory at the laser beam
	implicit none
	double precision, intent(in) :: initialposition, velocity,posintersect, t0
	Double precision, intent(out) :: tarr, TOF
		
	!calcualte the time of arrival at the laser beam
	TOF=(posintersect-initialposition)/velocity
	tarr= ((posintersect-initialposition)/velocity)+t0
		
end subroutine TOFcalc
	
subroutine TOF(tarr,tarrmin,tarrmax,initialposition, velocity,posintersect, t0)
!Calcualtes the time of arrivals of each trajectory at the laser beam
	implicit none
	double precision, intent(in) :: initialposition, velocity,posintersect, t0
	Double precision, intent(out) :: tarr, tarrmax, tarrmin

	!calcualte the time of arrival at the laser beam
	tarr= ((posintersect-initialposition)/velocity)+t0
		
				
	!If statement to find the largest and smallest time of arrival
		if (tarr .gt. tarrmax) then
			tarrmax=tarr
		elseif (tarrmin .gt. tarr) then
			tarrmin=tarr
		end if
			
end subroutine TOF

!Subroutine  TOFcount(nintersect,tarrmax,tarrmin)
!!Calcualtes whether each trajectory passes through the laser beam
	!implicit none
	!double precision :: initialposition, velocity,posintersect, t0
	!integer :: j, i, k, l, m
	!integer, intent(in) :: nintersect
	!Double precision :: tarrmax, tarrmin,tarr1
	!Double precision :: wbin !number and width of bin
	!Double precision, dimension(:),  allocatable :: tarr,valuebin,countbin,countbinfinal
	!Integer :: nbin

	!open(unit=13, file='TOF.txt')
	!open(unit=14, file='TOFcount.txt')
	!open(unit=15, file='tarrminmax.txt')	
	!open(unit=16, file='TOFcountotal.txt')
	
	!allocate(tarr(nintersect))

	!wbin=0.5
	!tarrmax=(int(tarrmax*1E6))+6
	!tarrmin=int(tarrmin*1E6)-6
	!nbin=(tarrmax-tarrmin)/wbin

	!write(15,*) tarrmax
	!write(15,*) tarrmin
	
	!allocate(valuebin(nbin))
	!allocate(countbin(nbin))	
	!allocate(countbinfinal(nbin))

	!reads in times of arrival converts from s to us and puts in an array
	!do j=1, nintersect
		!Read(13,'(i,1e20.6)') m, tarr1
		!tarr(j)=tarr1*1E6	
	!end do
	
	!!gives values to bins based on tarrmin and wbin
	!do i=1, nbin
	!valuebin(i)=tarrmin+i*wbin
	
	!!bins each arrival time in a bin with a max
		!do l=1, nintersect
		!If (tarr(l) .gt. valuebin(i)) then
		!countbin(i)=countbin(i)+1
		!end if
		!end do
	!end do
	
	
	!do k=1, nbin
	!write(16,'(2e20.6)') valuebin(k),countbin(k)
	!countbinfinal(k)=countbin(k)-countbin(k+1)
	!write(14,'(2e20.6)') valuebin(k),countbinfinal(k)
	!end do
	
!end subroutine TOFcount

Subroutine Dopplershift(velocity,Frequency)
!Calcualtes the dopplershift of trajectories
	implicit none
	Double precision, intent(in) :: velocity
	Double precision, intent(out) :: Frequency
	Double precision :: Frequency0

	Frequency0=2.998E8/790E-9
	frequency=velocity/2.998E8*Frequency0
	
	
end subroutine Dopplershift


