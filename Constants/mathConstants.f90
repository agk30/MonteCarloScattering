module mathConstants
    implicit none

    integer, parameter :: r14 = selected_real_kind(14,30)
    real(kind=r14), parameter :: pi = 3.141592653589793D0
    real(kind=r14), parameter :: boltzmannConstant =1.38064852D-23

    public :: r14, pi, boltzmannConstant

end module mathConstants