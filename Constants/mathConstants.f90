module mathConstants
    implicit none

    ! General use mathematical constants for use in all subroutines
    ! The r14 parameter is very important and the program will not compile without it
    ! It defines the kind for the real variables as double precision (14 decimal places, max exponent of 30)
    ! This method is more portable than the double precision variable
    integer, parameter :: r14 = selected_real_kind(14,30)
    real(kind=r14), parameter :: pi = 3.141592653589793D0
    real(kind=r14), parameter :: boltzmannConstant =1.38064852D-23

    public :: r14, pi, boltzmannConstant

end module mathConstants