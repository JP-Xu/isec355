import os

def get_files(path='./', pre = 'gamma', ext = 'dcd'):
    ''' This function returns files under current folder with defined prefix and extension.
    input pre: string, prefix of wanted files
    input ext: string, extension of wanted files.
    return: a list of file names.
    '''
    import regex as re
    import natsort 
    
    result = []
    for file_ in os.scandir(path):
        if file_.is_file() and re.search(r'^(%s)'%pre, file_.name) and re.search(r'(%s)$'%ext, file_.name):
            result.append(file_.name)
    result = natsort.natsorted(result, alg=natsort.FLOAT)        
    return result

def get_dirs(path='./', pre = '65', end = 'input'):
    ''' This function returns files under current folder with defined prefix and extension.
    input pre: string, prefix of wanted files
    input ext: string, extension of wanted files.
    return: a list of file names.
    '''
    import regex as re
    import natsort 
    result = []
    for file_ in os.scandir(path):
        if file_.is_dir() and re.search(r'^(%s)'%pre, file_.name) and re.search(r'(%s)$'%end, file_.name):
            result.append(file_.name)
    result = natsort.natsorted(result, alg=natsort.FLOAT)        
    return result

all_folders = get_dirs('/work/hung_group/xu.jiam/liposomes/namd_sims/DOPC_DS', pre='DS', end='.')
for fd in all_folders:
    os.chdir('/work/hung_group/xu.jiam/liposomes/namd_sims/DOPC_DS/'+fd+'/namd/')
    # copy inp, slurm and py files to current folder.
    os.system("cp /work/hung_group/xu.jiam/liposomes/namd_sims/Results/gamma*inp ./")
    os.system("cp /work/hung_group/xu.jiam/liposomes/namd_sims/Results/pair_cal.py ./")
    os.system("cp /work/hung_group/xu.jiam/liposomes/namd_sims/Results/pair_cal.slurm ./")
    os.system("cp /work/hung_group/xu.jiam/liposomes/namd_sims/Results/toppar/* ./toppar/") 
    # change job name based on current dir
    data = fd.split('-')
    job_name = fd #data[0]+data[1][:2]+data[2][:2]
    
    os.system(f"sed -i '/#SBATCH -J/s/.*/#SBATCH -J {job_name}/g' pair_cal.slurm")
    os.system("sbatch pair_cal.slurm")
    #print("System" + fd + "has been finished.", end='')
