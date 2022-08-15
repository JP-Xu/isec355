from main import *
import os
import MDAnalysis as mda
import natsort
import csv

# change working path
# automatically and save to csv file with header:
# apl  se  composition  gamma
csv_header = ['APL', 'Error', 'Composition', 'Gamma'] 
csv_data = []

# Folder names
download_folder = "/Users/jiamingxu/Downloads/"
folders = get_bilipid_pathes(download_folder)
apl_results = []
## Looping all folders just scaned

for folder in folders:
    print(folder)
    os.chdir(folder)
    
    files = os.listdir()
    dcdfiles = get_files('dcd')
    psffiles = get_files('psf')
    gammas = [-7, 0, 7, 15]
    ts = 2 #int(input("Enter the timestep of simulation in fs/step: "))

    for gamma in gammas:

        dcdfiles_in_gamma = [ file for file in dcdfiles if file.startswith('gamma' + str(gamma)) ]
        ## Check if there is one or more than one psf file in the path.

        #if len(psffiles) == 1:
        #    
        #    ## if dcd files exist.
        #    if dcdfiles_in_gamma != []:
        #        u = mda.Universe(psffiles[0], dcdfiles_in_gamma)
        #    else:
        #        print("DCD files for gamma {} doesn't exist. Composition: {}".format(gamma, folder))
        #
        ## Let user select which psf file that user wants to process.
        #elif len(psffiles) > 1:
        #    string_for_query = "Enter which psf file you want to process: \n\n"
        #    for i, file in enumerate(psffiles):
        #        string_for_query += str(i) + " " + file + "\n"
        #    input_index = input(string_for_query)
        #    while True:
        #        try:
        #            input_index = int(input_index)
        #            assert input_index in range(len(psffiles)), "Input invalid, index must be the one in above."
        #            break
        #        except:
        #            print("Input invalid, try again.")
        #        
        #    u = mda.Universe(psffiles[input_index], dcdfiles_in_gamma)

        # Raise error if no psf file exists
        #else:
        #    raise  FileNotFoundError("a psf file must be placed in the same path.")
        #u = mda.Universe("step5_input.psf", dcdfiles_in_gamma)
        ## using isec355 to get apl and standard error.
        lobby = Isec355("step5_input.psf", dcdfiles_in_gamma)
        apl_result = lobby.get_apl()
        apl_results.append(list(apl_result) + [folder] + [gamma])

    os.chdir("../")

with open('apl.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(csv_header)

    # write multiple rows
    writer.writerows(apl_results)