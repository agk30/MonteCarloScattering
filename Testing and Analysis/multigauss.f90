include "../src/hwlib/speeds.f90"
include "../src/hwlib/mathConstants.f90"

program mg

    use speeds
    use mathConstants

	double precision, dimension(:), allocatable :: m_s, w_s, std_s
    double precision, dimension(:), allocatable :: m_t, w_t, std_t
    double precision :: dist, pulseLength, speed, t0
    logical :: gauss_time
    integer :: n_s, n_t

    call random_seed

    gauss_time = .False.
    n_s = 1
    n_t = 1
    pulseLength = 10E-6
    dist = 0.178

    allocate(m_s(1))
    allocate(w_s(1))
    allocate(std_s(1))
    allocate(m_t(1))
    allocate(W_t(1))
    allocate(std_t(1))

    m_s(1) = 87
    w_s(1) = 1
    std_s(1) = 10
    m_t(1) = 12
    w_t(1) = 1
    std_t(1) = 2

    call ingoing_speed_from_Gauss(w_s, m_s, std_s, w_t, m_t, std_t, n_s, n_t, gauss_time, dist, pulseLength, speed, t0)

    print *, speed, t0

end program mg