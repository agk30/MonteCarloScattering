program mbtest
    implicit none

    double precision :: rand, rand2, maxSpeed1D, temp, mass, boltzmannConstant, pi, scatteredSpeed, &
        part1, part2, probability
    integer :: i, j
    logical :: accepted
    integer, dimension(50) :: array

    temp = 10D0
    mass = 2.8240519D-26
    boltzmannConstant = 1.38064852D-23
    pi = 3.141592653589793D0
    maxSpeed1D = 500D0

    open(unit=12,file= 'MBDistribution.txt')

    array = 0

    call random_seed

    accepted = .FALSE.

    do i = 1, 10000
        do
            call random_number(rand)

            scatteredSpeed = rand*maxSpeed1D

            !part1 = SQRT(mass/(2D0*pi*boltzmannConstant*temp))
            part2 = -(mass*scatteredSpeed*scatteredSpeed)/(2D0*boltzmannConstant*temp)
            probability = EXP(part2)

            !print *, probability, scatteredSpeed

            !print *, mass, scatteredSpeed, probability
            call random_number(rand2)

            if (probability .gt. rand2) then
                EXIT
            end if
        end do

        !print *, probability, scatteredSpeed

        do j = 1, 50
            if ((scatteredSpeed .gt. (j-1)*10) .and. (scatteredSpeed .lt. j*10)) then
                array(j) = array(j) + 1
                EXIT
            end if
        end do
    end do

    do i = 1, 50
        write(12,*) array(i)
    end do

end program mbtest