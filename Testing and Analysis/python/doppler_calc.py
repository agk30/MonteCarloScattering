def doppler_shift(frequency, velocity):
    return frequency * (velocity/2.99792458e8)

c = 2.99792458e8
v_t = 1000
wavelength = 283E-9

frequency = c/wavelength

change = doppler_shift(frequency, v_t)

shifted_frequency = (frequency + change)
shifted_waveltength = c/shifted_frequency

print (str(change/10E6)+" MHz")