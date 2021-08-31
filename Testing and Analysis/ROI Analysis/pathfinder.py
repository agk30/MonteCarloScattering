import os

# Path to data must be provided
path="D:/Dev/Data/0 degrees/sqe/Averages/"
datasets = []
# Lists are initialised here so that values can be assigned to their indices later
surface_out = ["","",""]
pfpe = ["","",""]
sqa = ["","",""]
sqe = ["","",""]

counter = 1

# Finds all files and folders within given path.
# Takes a sample file name and stores it in the directories list
for (dirpath, dirnames, filenames) in os.walk(path):
    for i in filenames: 
        file_path = dirpath + "/" + i
        datasets.insert(counter,file_path)
        counter = counter + 1
        break

# Parses the appropriate transition and surface information from directory name
for set in datasets:
    prestr, transition, surface, endstr = set.split("_",3)
    prestr, date = prestr.rsplit("/",1)
    # surface, endstr = endstr.split(" ",1)

    # Sorts datasets into their respective surface lists.
    # To create a standard, Q12 data goes to the 0th index, Q13 to 1st, and Q14 to 2nd
    if surface == "Ingoing Beam":
        if transition == "Q12":
            surface_out[0] = set
        elif transition == "Q13":
            surface_out[1] = set
        elif transition == "Q14":
            surface_out[2] = set

    if surface == "PFPE":
        if transition == "Q12":
            pfpe[0] = set
        elif transition == "Q13":
            pfpe[1] = set
        elif transition == "Q14":
            pfpe[2] = set

    if surface == "SQA" or "Squalane":
        if transition == "Q12":
            sqa[0] = set
        elif transition == "Q13":
            sqa[1] = set
        elif transition == "Q14":
            sqa[2] = set

    if surface == "SQE" or "Squalene":
        if transition == "Q12":
            sqe[0] = set
        elif transition == "Q13":
            sqe[1] = set
        elif transition == "Q14":
            sqe[2] = set

# Creates system call to fortran program which processes the data
# Call must be in format of "<program name> <surface in data> <surface out data>"
for i in range(len(surface_out)):
    if pfpe[i] != "":
        command =  ("a.exe" + ' "' + pfpe[i] + '" "' + surface_out[i] + '"')
        print (command)
        os.system(command)
    if sqa[i] != "":
        command =  ("a.exe" + ' "' + sqa[i] + '" "' + surface_out[i] + '"')
        print (command)
        os.system(command)
    if sqe[i] != "":
        command =  ("a.exe" + ' "' + sqe[i] + '" "' + surface_out[i] + '"')
        print (command)
        os.system(command)
    print ("done Q1" + str(i+2) + " set")