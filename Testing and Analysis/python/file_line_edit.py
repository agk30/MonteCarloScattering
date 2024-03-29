import filemanagement as fm
import fileinput
import sys
import os

def replacement(file, previousw, nextw):
   for line in fileinput.input(file, inplace=1):
       line = line.replace(previousw, nextw)
       sys.stdout.write(line)

old_line = "fullSim = T"
new_line = "fullSim = F"

# opening the file in read mode
folder = fm.get_input_folder()

for root, dirs, files in os.walk(folder):
    for file in files:
        replacement(root+"/"+file, old_line, new_line)