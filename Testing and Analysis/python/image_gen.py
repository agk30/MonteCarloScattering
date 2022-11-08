import numpy

image = numpy.ones((420,420))

with open('test_image.txt','wb') as f:
    numpy.savetxt(f, image , fmt='%s', delimiter='  ')