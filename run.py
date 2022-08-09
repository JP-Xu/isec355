import string
from main import *
import os
import MDAnalysis as mda

# change working path
os.chdir("/Users/jiamingxu/Downloads/65DYDM")

files = os.listdir()
dcdfiles = [ file for file in files if file.endswith('dcd') ]
psffiles = [ file for file in files if file.endswith('psf') ]
gammas = [-7, 0, 7, 15]
ts = int(input("Enter the timestep of simulation in fs/step: "))

for gamma in gammas:

    ## Check if there is one or more than one psf file in the path.

    if len(psffiles) == 1:
        dcdfiles_in_gamma = [ file for file in dcdfiles if 'gamma' + str(gamma) in file]
        u = mda.Universe(psffiles[0], ",".join(dcdfiles_in_gamma) )

    # Let user select which psf file that user wants to process.
    elif len(psffiles) > 1:
        string_for_query = "Enter which psf file you want to process: "
        for i, file in enumerate(psffiles):
            string_for_query += str(i) + " " + file + "\n"
        input_index = input(string_for_query)
        while True:
            try:
                input_index = int(input_index)
                assert input_index in range(len(psffiles)), "Input invalid, index must be the one in above."
                break
            except:
                print("Input invalid, try again.")
            
        u = mda.Universe(psffiles[input_index])

    # Raise error if no psf file exists
    else:
        raise  FileNotFoundError("a psf file must be placed in the same path.")

    ## using isec355
    lobby = isec355(u, ts)
    print( lobby.get_apl())

#u = mda.Universe()