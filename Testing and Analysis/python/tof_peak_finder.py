import scipy.interpolate
import numpy
import os
import filemanagement as fm
import matplotlib.pyplot as plt

folder_path = fm.get_input_folder()

#peak_array[wedges, arcs(+1 extra for wedge numbers)]
peak_array = numpy.zeros((12,8))

arc_list = [1, 2, 3, 4, 5, 6, 7]

arc_list_length = len(arc_list)

wedge_list = []

i = 0
for root, dirs, files in os.walk(folder_path):
    for name in files:
        if name.endswith(".csv"):
            with open(root+'/'+name,'r') as f:
                data = numpy.loadtxt(f, delimiter=',')

                numbers = fm.grab_number(name)
                wedge_number = int(numbers[0])
                wedge_list.append(wedge_number)
                print (wedge_number)

                for arc in arc_list:

                    #NB: Addresses are given as [row,column]
                    profile = data[:,[0,arc+arc_list_length]]

                    fit = scipy.interpolate.InterpolatedUnivariateSpline(profile[:,0],profile[:,1], k=4)
                    

                    roots = fit.derivative().roots()

                    fit_vals = fit(roots)



                    if fit_vals.any():
                        max_index = numpy.argmax(fit_vals)
                        max_val = roots[max_index]
                    else:
                        max_val = 0

                    if wedge_number == 0:
                        print(max_val, i, arc)
                    
                    peak_array[i,arc] = max_val

                    if wedge_number == 0 and arc == 2:
                        xs = numpy.linspace(0, 200E-6, 1000)
                        plt.ylim(ymax = 1, ymin = 0)
                        plt.plot(profile[:,0], profile[:,1], 'r')
                        plt.plot(xs, fit(xs), 'b')
                        fit.set_smoothing_factor(10)
                        plt.plot(xs, fit(xs), 'g')
                        plt.show()

                i += 1

for i in range(len(wedge_list)):
    peak_array[i,0] = wedge_list[i]

#output_list = numpy.stack((wedge_list,peak_array))
#output_list = numpy.swapaxes(peak_array, 0, 1)

with open ('peak_list.csv','w') as f:
    numpy.savetxt(f, peak_array, delimiter=',')

