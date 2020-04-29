include "Maths/BoltzmannDistribution.f90"

!Integrating a TOF appearance profile (normalized to 0,1) gives a cumulative distribution function. Taking an inverse function of normalized (to 0,1) cumulative distribution
!function gives function ArrivalTime. By inputting a random number between 0 and 1, a corresponding arrival time is assigned. This program assignes arrival time from a real
!distribution taken for a single-point LIF measurement at valve-probe laser distance of 230 mm. As it's a real distribution, the spread due to differences in time of creation of each
!individual trajectory is already accounted for.
program getspeeds

	use class_BoltzmannDistribution

	implicit none
		
	integer :: n,ncyc
	double precision :: x, x0, h, s, Amax, Amin, dist, dischargepulse, t, mass, Temp
	double precision, dimension(:),  allocatable :: Speed, ArrivalTime, t0
	logical :: thermalDesorption, hit
	double precision :: probability, mostLikelyProbability, normalisedProbability
	double precision :: rand1, rand2
	double precision :: MBspeed
	double precision :: maxSpeed
	double precision :: mostProbableSpeed, BoltzmannConstant
	maxSpeed = 000
	BoltzmannConstant = 1.38064852E-23

	open(unit=10,file='getspeeds.txt')
	open(unit=11,file='getspeeds.inp')
	open(unit=12,file='getdirections.inp')
	open(unit=13,file='MBSpeeds.txt')
		
	read(11,*) x0
	read(11,*) Amax
	read(11,*) Amin
	read(11,*) h
	read(11,*) s
	read(11,*) dist
	read(11,*) dischargepulse
	read(11,*) thermalDesorption
	read(11,*) mass
	read(11,*) Temp
	read(12,*) ncyc
		
	allocate(Speed(ncyc))
	allocate(ArrivalTime(ncyc))
	allocate(t0(ncyc))

	if(thermalDesorption == .FALSE.) then

		write(10,'(3a)') "ncyc, Time of creation /s, Speed /ms-1"

		do n=1,ncyc
			
			!Call random number and use it to assign a random time of creation in the discharge pulse
			call random_number(t)
			t0(n)=t*dischargepulse
			
			!Call random number and use it to assign a time of arrival of a trajectory at a specified distance from valve using a real TOF appearance profile
			call random_number(x)
			ArrivalTime(n) = x0/(((Amax-Amin)/(x-Amin))**(1/s)-1)**(1/h)	!Inverse of a logistic curve (Logistic5 in Origin) resultant from integrating and normalising (0, 1) of the TOF profile
																			!This convolutes the spread in speeds due to imperfect cooling at expansion and due to the length of the discharge pulse
																			!The length of the discharge pulse is therefore artificially decreased from the nominal 10 microseconds
			Speed(n) = dist/(ArrivalTime(n)*1D-6)							!Tranforiming time of arrival into a speed
			
			write(10,'(i7,3e20.6)') n, t0(n), Speed(n)

		end do

	else

		write(13,'(3a)') "ncyc, Speed /ms-1"

		mostProbableSpeed = sqrt((2*BoltzmannConstant*Temp)/mass)
		!print *, mostProbableSpeed

		mostLikelyProbability = MBProbability(Temp, mostProbableSpeed, mass)

		do n=1, ncyc

			hit = .FALSE.

			do while (hit == .FALSE.)

				call random_number(rand1)
				MBspeed = rand1*maxSpeed

				probability = MBProbability(Temp, MBspeed, mass)
				!calculates the probability of the speed with respect to the most probable speed equalling 1. The Maxwell-Boltzmann distribution is already normalised to 1, meaning that the sum of all probabilities from zero to infinity will equal 1.
				normalisedProbability = probability/mostLikelyProbability

				call random_number(rand2)

				if (normalisedProbability .gt. rand2) then

					hit = .TRUE.

					write(13,'(i7,3e20.6)') n, MBspeed

				end if

			end do
		
		end do

	end if

end program getspeeds