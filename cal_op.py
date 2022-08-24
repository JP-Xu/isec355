from tqdm import tqdm
from main import *
import os
import MDAnalysis as mda
import natsort
import csv



def append_list_as_row(file_name, list_of_elem):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)


# change working path
# automatically and save to csv file with header:
# apl  se  composition  gamma

# Folder names
target_folder = input("Enter the target folder: ")
num_acyl_catoms = input("What's the length of the longest acyl chain? ")
num_acyl_catoms = int(num_acyl_catoms)
folders = get_bilipid_pathes(target_folder)
acyl_cnames = [ 'C'+str(x+2) for x in range(num_acyl_catoms-2)]
csv_header = ['composition', 'gamma', 'lipid', 'SN'] +  acyl_cnames
csv_data = []
## Looping all folders just scaned
#folders=['100DOPC','100DYPC']


for folder in tqdm(folders):
    print("Processing compisition {}, {} outof {} folders.".format(folder, folders.index(folder), len(folders)), end='\r')
    os.chdir(folder + '/namd')
    
    files = os.listdir()
    dcdfiles = get_files('dcd')
    psffiles = get_files('psf')
    gammas = [-7, 0, 7, 15]
    ts = 2 #int(input("Enter the timestep of simulation in fs/step: "))
    
    for gamma in gammas:

        dcdfiles_in_gamma = [ file for file in dcdfiles if file.startswith('gamma' + str(gamma)) ]
        try:
            lobby = Isec355("step5_input.psf", dcdfiles_in_gamma)
        except:
            print("Unable to analyze for system {}, gamma {}".format(folder, gamma))
            continue
        names_d, result_d = lobby.order_param(-500, -1) # last 500 frames, 100 ns.
        for l in names_d.keys():
            for i in range(2):
                csv_data.append([folder, gamma, l, 'SN'+str(i+1)] + list(result_d[l][i]))
        #print(names_d)
        #print(result_d)
    os.chdir("../../")
try:
    os.chdir('Results')
except:
    print("Failed to enter the results folder, save file to {}.".format(os.getcwd()))


with open('op_.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(csv_header)

    # write multiple rows
    writer.writerows(csv_data)
