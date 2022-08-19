from main import *
import os
import MDAnalysis as mda
import natsort
import csv

# change working path
# automatically and save to csv file with header:
# apl  se  composition  gamma
csv_header = ['std', 'Composition', 'Gamma'] 
csv_data = []

# Folder names
target_folder = input("Enter the target folder: ")
folders = get_bilipid_pathes(target_folder)
## Looping all folders just scaned

for folder in folders:
    print(folder)
    os.chdir(folder + '/namd/')
    
    files = os.listdir()
    dcdfiles = get_files('dcd')
    psffiles = get_files('psf')
    gammas = [-7, 0, 7, 15]
    ts = 2 #int(input("Enter the timestep of simulation in fs/step: "))
    
    for gamma in gammas:

        dcdfiles_in_gamma = [ file for file in dcdfiles if file.startswith('gamma' + str(gamma)) ]
        lobby = Isec355("step5_input.psf", dcdfiles_in_gamma)
        
        p_l1_result = np.array([])
        p_l2_result = np.array([])
        p_l1_std    = np.array([])
        p_l2_std    = np.array([])

        for frame in lobby.trajectory[:]:
            P_select_string = "(resname " + " or resname ".join([ x for x in set(lobby.residues.resnames) if 'PC' in x ]) + ') and name P' 
            try:
                Isec355.repair_bilayer(lobby, frame)
            
            except FloatingPointError:
                print("Cannot reshape bilayer, skip frame {} for gamma {}, composition {}.".format(frame.frame, gamma, folder) )
                continue
            
            P_coors_l1 = lobby.select_atoms(P_select_string + ' and resid 1:100').positions
            P_coors_l2 = lobby.select_atoms(P_select_string + ' and resid 101:200').positions
            p_l1_result = np.append(p_l1_result, np.average(P_coors_l1[:,2]))
            p_l2_result = np.append(p_l2_result, np.average(P_coors_l2[:,2]))
            p_l1_std = np.append(p_l1_std, np.std(P_coors_l1[:,2]))
            p_l2_std = np.append(p_l2_std, np.std(P_coors_l2[:,2]))
        
        # error propagation, sigma = |mean|*sqrt(sum( (sigma_i/mean_i)^2 )))
        p_l1_std = np.average(p_l1_result) * np.sqrt(sum((p_l1_std / p_l1_result) ** 2))
        p_l2_std = np.average(p_l2_result) * np.sqrt(sum((p_l2_std / p_l2_result) ** 2))


        csv_data.append( [p_l1_std, p_l2_std] + [folder] + [gamma])
    os.chdir("../../")
try:
    os.chdir('Results')
except:
    print("Failed to enter the results folder, save file to {}.".format(os.getcwd()))
with open('surface_std.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(csv_header)

    # write multiple rows
    writer.writerows(csv_data)