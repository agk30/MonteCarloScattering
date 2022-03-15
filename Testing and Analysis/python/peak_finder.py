import numpy
import scipy.signal
import matplotlib.pyplot as plt
from tkinter import Tk     # from Tkinter import Tk for earlier than Python 3.x
from tkinter.filedialog import askopenfilename

# we don't want a full GUI, so keep the root window from appearing
Tk().withdraw() 
# show an "Open" dialog box and return the path to the selected file
input_path = askopenfilename() 

min_peak_height = 0.0015

# try to use 'with open()' as much as you can for handling files, it essentially makes sure the file is closed properly once you are done using it
with open(input_path,'r') as file:

    # take spectrum file and take only the list of wavelengths and the background subtracted column
    array = numpy.loadtxt(file)
    array_shape = numpy.shape(array)
    spectrum = numpy.zeros((array_shape[0],2))
    spectrum[:,0] = array[:,0]
    spectrum[:,1] = -array[:,3]

    # extracts peaks. note that the resulting array is the indices in the spectrum array which contain the peaks, not the peaks themselves
    peaks, _ = scipy.signal.find_peaks(spectrum[:,1], min_peak_height)

    # set up plots for spectrum and the peaks
    plt.plot(spectrum[:,0],spectrum[:,1])
    plt.plot(spectrum[peaks,0], spectrum[peaks,1], "x")

    # assign annotations to the peaks
    for i, txt in enumerate(peaks):
        plt.annotate(spectrum[txt,0], (spectrum[peaks[i],0], spectrum[peaks[i],1]))
        print(spectrum[txt,0], spectrum[i,1])
    
    # finally plot the graph
    plt.show()