import glob
import re

input_folder="mapping/logs/*"
files = glob.glob(input_folder)
ERRORCHECK = False

for file in files:
    # ERRORCHECK = False
    valsfile = open(file)
    for line in valsfile:  
        if "error" in line:
            print (file+" Alignment didn't work")
            