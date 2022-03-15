import os
import shutil
import numpy
import sys
from tkinter import Tk     # from Tkinter import Tk for earlier than Python 3.x
from tkinter.filedialog import askdirectory, askopenfilename

def read_image(image_path):

    with open(image_path, "r") as f:
        image = numpy.loadtxt(f)

    return image

def parse_file_name(file_path):

    name_list = file_path.split("_")
    surface_list = ["SQA","SQE","IB","PFPE","OA","Bkg","InstrumFunc","LOA"]
    transition_list = ["Q11","Q12","Q13","Q14","Q15"]

    junk = []

    for name in name_list:
        if name in surface_list:
            surface = name
        elif name.isnumeric():
            delay = name
        elif name in transition_list:
            transition = name
        else:
            junk.append(name)

    delay = int(simple_split(file_path,"ChC"))
    
    if 'surface' not in locals():
        surface = -1

    if 'delay' not in locals():
        delay = -1

    if 'transition' not in locals():
        transition = -1

    return surface, delay, transition

def simple_split(file_path, delimiter):

    file_path = file_path.split(".")
    stem = file_path[0]
    name_list = stem.split(delimiter)

    found_delay = False

    for name in name_list:
        if name.isnumeric():
            delay =  name
            found_delay = True
            
    if found_delay:
        return delay
    else:
        print ("No delay found in image file name, did you use the right delimiter?")
        sys.exit()

# creates output folders if they do not already exist
def dir_setup(directory):

    if not os.path.isdir(directory):
        os.makedirs(directory)
    else:
        shutil.rmtree(directory)
        os.mkdir(directory)

def get_input_folder():

    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw() 
    # show an "Open" dialog box and return the path to the selected file
    input_path = askdirectory() 
    
    return input_path

def get_input_file():

    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw() 
    # show an "Open" dialog box and return the path to the selected file
    input_path = askopenfilename() 

    return input_path

def get_indices(surface, transition):

    surface_list = ["SQA","SQE","PFPE","OA","Bkg","InstrumFunc","LOA"]
    transition_list = ["Q12","Q13","Q14","Q15"]

    if surface == "IB":
        surface_index = 0
    elif surface in surface_list:
        surface_index = 1
    else:
        surface_index = -1

    if transition in transition_list:
        transition_index = transition_list.index(transition)
    else:
        transition_index = -1

    return surface_index, transition_index