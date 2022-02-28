import roi
import numpy
import os
import shutil

output_dir = "Output Images/Summed Images/Generic/"
shutil.rmtree(output_dir)
os.mkdir(output_dir)

folder_list = []

folder_list.append(r"D:\Valve Distance Variation Experiment\80mm Back\05072021_1_Q11_IB_TOF Profile")
folder_list.append(r"D:\Valve Distance Variation Experiment\80mm Back\05072021_2_Q11_IB_TOF Profile")
folder_list.append(r"D:\Valve Distance Variation Experiment\80mm Back\05072021_3_Q11_IB_TOF Profile")

summed_images, delay_list = roi.sum_tofs(folder_list, "ChC")

background = summed_images[:,:,0].copy()

for i in range(len(delay_list)):
    summed_images[:,:,i] = summed_images[:,:,i] - background


for i in range(len(delay_list)):
    file_name = "summed_image_"+str(delay_list[i])
    with open(output_dir+file_name, 'wb') as f:
        numpy.savetxt(f, summed_images[:,:,i], fmt='%.5e')

#print (summed_images[0,0,0])