import os

path="D:/Dev/Data/0 degrees/sqa/Sequences/"
datasets = []
surface_out = ["","",""]
pfpe = ["","",""]
sqa = ["","",""]
sqe = ["","",""]

counter = 1

for (dirpath, dirnames, filenames) in os.walk(path):
    for i in filenames: 
        file_path = dirpath + "/" + i
        datasets.insert(counter,file_path)
        counter = counter + 1
        break

for set in datasets:
    prestr, run, transition, endstr = set.split("_",3)
    prestr, date = prestr.rsplit("/",1)
    surface, endstr = endstr.split(" ",1)

    if surface == "IB":
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

    if surface == "SQA":
        if transition == "Q12":
            sqa[0] = set
        elif transition == "Q13":
            sqa[1] = set
        elif transition == "Q14":
            sqa[2] = set

    if surface == "SQE":
        if transition == "Q12":
            sqe[0] = set
        elif transition == "Q13":
            sqe[1] = set
        elif transition == "Q14":
            sqe[2] = set

for set in surface_out:
    for i in range(len(surface_out)):
        print (pfpe[i] + surface_out[i])
        break