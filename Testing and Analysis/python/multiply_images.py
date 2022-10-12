import numpy
import filemanagement as fm

print ("Select Image File")
image_path = fm.get_input_file()
image = fm.read_image(image_path)

print ("Select Instrument Functinon File")
if_image_path = fm.get_input_file()
if_image = fm.read_image(if_image_path)

if_image = if_image/numpy.amax(if_image)

image = image*if_image
numpy.savetxt("Output Data/if_adjusted_image.txt", image)