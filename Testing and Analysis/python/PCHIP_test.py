import scipy.interpolate
import numpy

with open(r"C:\Users\adam\Documents\wedge 0.0.csv",'r') as f:
    data = numpy.loadtxt(f, delimiter=',')

#NB: Addresses are given as [row,column]
profile = data[:,[0,11]]

print (profile.shape)

fit = scipy.interpolate.InterpolatedUnivariateSpline(profile[:,0],profile[:,1], k=4)

roots = fit.derivative().roots()

fit_vals = fit(roots)

max_index = numpy.argmax(fit_vals)
max_val = roots[max_index]

print (max_val)