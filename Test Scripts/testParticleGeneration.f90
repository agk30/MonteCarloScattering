include "../Constants/mathConstants.f90"

program testParticle
    use mathConstants
    
    implicit none

    integer :: i
    real(kind=r14), dimension(3) :: particleVector, particleStartPos
    real(kind=r14) :: particleSpeed, particleTime, rand1, rand2, start

    open(unit=11,file='testParticlesUp.txt')

    write(11, '(7a)') "i, xstart, ystart, zstart, speed, xvector, yvector, zvector"

    do i = 1, 100000
    
        call random_number(rand1)
        call random_number(rand2)

        ! randomly distributes particles across a 10cmx10cm square centred at origin in x and z, but a fixed start point of 5cm in y axis
        particleStartPos(2) = -0.05D0
        particleStartPos(1) = rand1*0.1D0 - 0.05D0
        particleStartPos(3) = rand2*0.1D0 - 0.05D0

        particleSpeed = 800

        ! all particles travelling directly 'up' (in y direction)
        particleVector = 0D0
        particleVector(2) = 1D0

        write(11, '(i7,7e20.8)') i, particleStartPos(1), particleStartPos(2), particleStartPos(3), particleSpeed, particleVector(1), particleVector(2), particleVector(3)

    end do


end program testParticle