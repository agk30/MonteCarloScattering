import filemanagement as fm
import numpy

file = fm.get_input_file()
output_file = "slice.txt"

row_range = [200, 205]
#row_range = [210, 215]
#row_range = [180, 185]
column_range = [0, 419]

with open(file,'r') as f:
    # remember: address for image is image(row,column) also can be seen as image(y,x)
    image = numpy.loadtxt(f)
    #bg = numpy.loadtxt("F:/Analysis Data Holding Cell/80mm Back/Summed Images/summed_image_100")
    #print (sum(bg))
    #image = image - bg
    slice = image[row_range[0]:row_range[1], column_range[0]:column_range[1]]
    print(image[200,220])

summed_slice = numpy.zeros((column_range[1] - column_range[0] + 1))

for i in range(column_range[1] - column_range[0]):
    summed_slice[i] = numpy.sum(slice[:,i])


with open(output_file, 'w') as f:
    numpy.savetxt(f, summed_slice)

