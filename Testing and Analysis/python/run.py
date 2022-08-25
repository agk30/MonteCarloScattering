import os

input_dir = "inputs/"

for root, dirs, files in os.walk(input_dir):
    for name in files:
        path = input_dir+name
        print("")
        print("Starting new run: "+path)
        os.system("./MCScattering "+'"'+path+'"')