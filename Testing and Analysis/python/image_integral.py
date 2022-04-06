import filemanagement as fm
import os
import numpy
import re
import matplotlib.pyplot as plt

# Goal is to take integral intensity across some defined pixel range in images as functinon of time

# Inputs start here, these are things that will be changed by the user, should really be taken as cl arguments but we will get there one day
delay_start = 95
delay_end = 159
tstep = 4
reactive_surface = "LOA"
output_folder = "Output Data/Integral Intensities/"
# limits for selecting region in image to integrate over
x_lim = [100, 500]
y_lim = [100, 500]
surface_list = ["SQA","SQE","IB","PFPE","OA","Bkg","InstrumFunc","LOA"]
transition_list = ["Q11","Q12","Q13","Q14","Q15"]

# Sets up directory structure if it does not already exist
fm.dir_setup(output_folder)

# constructs a list of all the delays to be analysed based on user input start time, end time and timestep
delay_list = []
for i in range(int((delay_end-delay_start)/tstep)+1):
    delay_list.append(delay_start+(i*tstep))

# asks the user for the folder containing the profiles we wish to analyse. Profiles must be contained within directories inside this main directory
input_folder = fm.get_input_folder()

# compiles list of directories (read: profiles) contained in this main directory
dir_list = []
for root, dirs, files in os.walk(input_folder):
    for dir in dirs:
        dir_list.append(root+"/"+dir)

# initialises the image array. Array is quite large and can potentially use a lot of memory so there might be a better way of doing this if it becomes a probelem
# indices are as follows: ([surface out = 0, surface in = 1], [number of delays in sequence], [number of transitions], [number of image x-pixels], [number of image y-pixels]
image_array = numpy.ones((2,int((delay_end - delay_start)/tstep)+1,4,677,630))

# begins process of analysing images
successful_images = []
# loops over all folders (profiles or image sequences) found before
for dir in dir_list:
    # loops over every file in the directory
    for file in os.listdir(dir):
        # funky regex to extract meta data from file name. It shouldn't be like this, we need to get meta data some other way ideally
        meta_data = re.findall("(LOA|IB|PFPE|OA)_(Q12|Q13|Q14|Q15)_(v1|v2|v3|v4|v5|v6|v7|v8|v9|v10)_ChC([0-9]{3})",file)
        # unpacks the tuple from regex function to get the right variables
        surface, transition, run, delay = meta_data[0]
        # turns the delay string into an integer
        delay = int(delay)

        # uses the meta data to find the right index in the image array this sequence should modify
        surface_index, transition_index = fm.get_indices(surface, transition, surface_list, transition_list)

        # since the data directories we typically create in the lab contain other stuff than just TOF profiles, gotta check that we only try to analyse the right things
        # check the get_indices function for how indices are assigned, but essentially -1 is assigned if a proper surface or transition is not found
        if (surface_index != -1) and (transition_index != -1):
            if (surface == reactive_surface) or (surface == "IB"):
                print ('Processing '+file,end='\r')
                # working_image is used to hold the image data before adding to the main array
                working_image = fm.read_image(dir+"/"+file)
                # all profiles of same surface and transition are summed together
                image_array[surface_index, delay_list.index(delay), transition_index, :, :] += working_image[:,:]
                # keeps track of all files that were successfully processed to help identify problems if things don't work how they are supposed to
                successful_images.append(dir+"/"+file)

# throws the successful_images file into the output folder
with open(output_folder+'successful_images.txt','w') as f:
    for file in successful_images:
        f.write(file+"\n")

# now for the maths
# pulls both the surface in and out images from the main array
# normally this should just point to them rather than creating new objects which would end up using a ton of memory
surface_in = image_array[1,:,:,:,:]
surface_out = image_array[0,:,:,:,:]

# subtraction leaves (hopefully) only scattered signal behind
subtracted_images = surface_in - surface_out

# loop over transitions, keep them separate though
for i in range(4):
    integral_intensities = []
    # loop over delays
    for j in range(len(delay_list)):
        integral_intensities.append(numpy.sum(subtracted_images[j, i, x_lim[0]:x_lim[1], y_lim[0]:y_lim[1]]))
    
    # output_file is reused for each iteration of the loop to hold the 2 dimensional array which will be written to individual files
    output_file = numpy.stack((delay_list, integral_intensities), axis=-1)

    # if file naming is not working for you, feel free to change the string into whatever you need, just make sure to change variables into strings
    with open(output_folder+reactive_surface+'_Q1'+str(i+2)+'_integral_intensities.txt','w') as f:
        numpy.savetxt(f, output_file, delimiter=",")    

#plt.plot(delay_list, integral_intensities)
#plt.show()

print ('')
print ('Done',end='\r')
print ('')