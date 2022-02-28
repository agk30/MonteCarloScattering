import numpy
import math

probe_height = 0
centre_point = []
num_arcs = 0
num_wedges = 0

xPx = 0
yPx = 0

max_num_rads = 0
max_num_wedges = 0

image = numpy.zeros((xPx,yPx))
radius = numpy.zeros((max_num_rads))
wedge = numpy.zeros((max_num_wedges))

# for every pixel, decides what sector it belongs to
for row in range(xPx):
    for column in range(yPx):
        # decides if pixel lies within largest semi-circle at all
        dist_from_centre = math.sqrt((row-centre_point[1])**2 + (column-centre_point[2])**2)
        if dist_from_centre < radius[max_num_rads-1]:

            # finds the arc in which the pixel lies
            for i in range(max_num_rads):
                if dist_from_centre < radius[i]:
                    selected_arc = i
                    # Arc found!

            # finds the wedge in which the pixel lies
            angle = math.acos(abs(centre_point[2]-column)/dist_from_centre)
            for j in range(max_num_wedges):
                if angle < wedge[j]:
                    selected_wedge = j