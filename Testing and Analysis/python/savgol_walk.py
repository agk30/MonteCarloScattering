import sys
import os
import numpy

sys.path.insert(1, 'C:/Users/adam/Documents/Code/MontePython/src')

import convolution as convolution

folder_path = "E:/270921/Q13IB/Q13IB"

xPx = 677
yPx = 630
kSize = 27

kernel = convolution.savgol_coef_array(kSize)

for root, dirs, files in os.walk(folder_path):
    for name in files:
        file_path = root + "/" + name
        print (name)
        with open(file_path, "r") as f:
            matrix = convolution.multi_sav_gol(xPx,yPx,kSize,f,kernel)
        with open('Output Images/'+'sv_'+name+'.txt','wb') as f:
            numpy.savetxt(f, matrix, fmt='%.5e')