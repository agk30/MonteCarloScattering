import roi
import filemanagement as fm
import numpy
import os
import sys
import matplotlib.pyplot as plt

# get folder path from command line argument
#folder_path = roi.get_args(sys.argv[1:])
folder_path = fm.get_input_folder()
print("Scanning "+ folder_path)

###############################################################
# parameters (these should really be given via command line but
# haven't gotten round to setting all that up yet ¯\_(ツ)_/¯ )
###############################################################

centre_point = [283, 209.5]
#centre_point = [294, 209.5]
num_arcs = 7
num_wedges = 12

xPx = 420
yPx = 420

max_radius = 140
#max_radius = 160

#startTime = 98
startTime = 74
#endTime = 208
endTime = 250
#endTime = 500
#endTime = 148
timeStep = 2

delimiter = "_"

output_directory = "Profiles"

bother_graphing = False

subtract_bg = False

################################################################
# End of input parameters
################################################################

temp_dir = os.path.split(folder_path)[0]
print (temp_dir)

output_directory = os.path.join(temp_dir, output_directory)

print (output_directory)

fm.dir_setup(output_directory)

# this next bunch of stuff sets up appropriate working variables based on input parameters
num_timepoints = int((endTime - startTime) / timeStep) + 1
half_wedge_step = 90/num_wedges
half_arc_step = (max_radius/num_arcs)/2

# arrays initialised
image = numpy.zeros((xPx,yPx))
bg_image = numpy.zeros((xPx,yPx))
outputArray = numpy.zeros((num_arcs,num_wedges,int((endTime-startTime)/timeStep)+1,2))

# generates lists of radius boundaries and wedge boundaries
radius = roi.generate_radii(num_arcs, max_radius)
wedge = roi.generate_wedges(num_wedges)

# takes a first look at the files in the given folder, finds out how many are in there
# then lists the number of files that will be processed given the input parameters
list = os.listdir(folder_path)
number_files = len(list)
print ('Number of files in directory = '+str(number_files))
print ('Number of files to process = '+str(num_timepoints)+', starting from '+str(startTime)+', ending at '+str(endTime))

# compiles list of surface out measurements
s_out_list = []
for file in list:
    #for part in file.split("_"):
    for part in file.split("Ch"):
        if part == "IB":
            s_out_list.append(file)

timepoint_list = []
#hacky way of making sure only images that fit the correct timestep are sampled
previous_value = startTime - timeStep
# Loops over every file in folder

first_run = False
bg_image = 0
for root, dirs, files in os.walk(folder_path):
    for name in files:
        file_path = root + "/" + name
        #surface, delay, transition = parse_file_name(file_path)
        delay = int(fm.simple_split(file_path, delimiter))
        if (int(delay) <= endTime) and (int(delay) >= startTime):
            if (int(delay) - previous_value) == timeStep:
                timepoint_list.append(delay)
                print ('Processing '+name,end='\r')
                # For each file, data are read into the image matrix
                image = fm.read_image(file_path)

                #if first_run:
                #    bg_image = image
                #    first_run = False

                #print(numpy.where(numpy.isclose(image, -1000)))
                image = image - bg_image
                # image goes to be processed, assigning the pixel intensity to the correct ROI
                outputArray[:,:,int((int(delay)-startTime)/timeStep),1] = roi.roi_assign(xPx, yPx, centre_point, radius, wedge, num_arcs, num_wedges, image)
                previous_value = delay

# comment here refers to comment that will go into output file.
# the comment just contains all the of the radii used in processing the image for reference
comment = str(radius-half_arc_step)

delay_list = numpy.zeros((num_timepoints,1))

# compiles a list of all the delays used in this process
for i in range(num_timepoints):
    delay_list[i] = (startTime + i*timeStep)*1E-6

# takes the ROI array and writes each wedge out into its own file
for j in range(num_wedges):
    with open(output_directory+'/wedge '+str(round(wedge[j]-90-half_wedge_step,2))+'.csv','wb') as f:
        write_array = outputArray[:,j,:,1]
        # axes are swapped because they just don't seem to be the right way round in the first place?
        write_array = numpy.swapaxes(write_array, 0, 1)
        # VERY important step. This prevents the write array from being modified when making the normalised array.
        # without specificially making a copy, the normalised_array object will just point to the write_array object and modify it unintentionally.
        normalised_array = write_array.copy()
        for i in range(num_arcs):
            max_value = max(write_array[:,i])
            if max_value > 0:
                normalised_array[:,i] = normalised_array[:,i]/max_value
            #else:
                #print ("Exiting: Image wasn't properly read, possible cause is the choice of delimiter for extracting delay from image.")
                #sys.exit()
        
        # staples the delay list to the ROI data
        write_array = numpy.hstack((delay_list,write_array))
        write_array = numpy.hstack((write_array,normalised_array))
        numpy.savetxt(f, write_array, fmt='%.5e', delimiter=',', header=comment)

# this step handles the writing of a 'parameters' file, contianing the radius and wedge boundaries used, once again for reference
wedge = numpy.array(wedge)
radius = numpy.array(radius)
# This part adjusts either the wedge or radius array to be the same size as the other, meaning that they can now be stacked (stapled) together in order to write them out to a file.
if num_arcs != num_wedges:
    if num_arcs > num_wedges:
        wedge = numpy.pad(wedge, (0,num_arcs-num_wedges))
        largest = num_arcs
    else:
        radius = numpy.pad(radius, (0,num_wedges-num_arcs))
        largest = num_wedges
else:
    largest = num_arcs

half_radius_increment = max_radius/num_arcs
for index, element in enumerate(radius):
    radius[index] -= half_radius_increment

array = numpy.zeros((largest,2))
array[:,0] = radius
array[:,1] = wedge
"""
print (" ")

for j in range(num_arcs):
    print("arc "+str(j))
    print(" ")
    for i in range(num_wedges):
        print(outputArray[j,i,0,1])
"""

# finally writes the parameters file
with open(output_directory+'/parameters.txt','wb') as f:
    header_str = 'start time = '+str(startTime)+', end time = '+str(endTime)+', num radii = '+str(num_arcs)+', num angles = '+str(num_wedges)
    numpy.savetxt(f, array , fmt='%s', delimiter=',',header=header_str)

if bother_graphing:
    real_data = numpy.loadtxt(r"C:\Users\adam\Desktop\Book2.csv", delimiter=',')

    probe_array = outputArray[3,6,:,1].copy()
    maximum = max(probe_array)

    plt.plot(delay_list, probe_array/maximum, delay_list, real_data[:,1])
    plt.show()

print ('')
print ('Done',end='\r')
print ('')