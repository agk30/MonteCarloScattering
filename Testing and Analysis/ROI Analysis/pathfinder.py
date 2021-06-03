import os
path = "D:/Dev/Data/"
dirs = []
q12 = []
q13 = []
q14 = []

for (dirpath, dirnames, filenames) in os.walk(path):
    dirs.extend(dirnames)

for dir in dirs:
    print (dir)

for dir in dirs:
    date, run, transition, endstr = dir.split("_",3)
    surface, endstr = endstr.split(" ",1)

    if transition == "Q12":
        if surface == "IB":
            q12[1] = dir
        elif surface == "PFPE":
            q12[2] = dir
        elif surface == "SQA":
            q12[3] = dir
        elif surface == "SQE":
            q12[4] = dir
    elif transition == "Q13":
        if surface == "IB":
            q13.insert(1,dir)
        elif surface == "PFPE":
            q13.insert(2,dir)
        elif surface == "SQA":
            q13.insert(3,dir)
        elif surface == "SQE":
            q13.insert(4,dir)
    elif transition == "Q14":
        if surface == "IB":
            q14[1] = dir
        elif surface == "PFPE":
            q14[2] = dir
        elif surface == "SQA":
            q14[3] = dir
        elif surface == "SQE":
            q14[4] = dir

print (q13)



