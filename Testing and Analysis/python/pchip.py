import scipy.interpolate
import filemanagement as fm
import numpy
import matplotlib.pyplot as plt

tof_profile_path = fm.get_input_file()
tof_profile = numpy.loadtxt(tof_profile_path, delimiter=",")

spline = scipy.interpolate.pchip(tof_profile[26:100,0], tof_profile[26:100,12])
spline_der = spline.derivative()
roots = spline.roots(discontinuity=False,extrapolate=False)

numpy.savetxt("roots.txt", roots)


#for root in roots:
#    if root < 177E-6 and root > 118E-6:
#        print(root)

#print (spline.derivative().roots())
#print (spline.derivative())

#print (spline_der.roots(discontinuity=False,extrapolate=False))
