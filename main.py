from ast import Is
import numpy as np
import os
import natsort
import MDAnalysis
import regex as re
from scipy import stats
import csv


class Isec355(MDAnalysis.core.universe.Universe):
    """
    An extension class for MDAnalysis.
    This class contains several functions:
    length of side of x-y plane.
    """
    
    def __init__(self, struc_file, traj_file, ts=2):
        """ 
        input Universe: MDAnalysis Universe.
        attributes:
        2. lipnames: lipid names in file.
        1. dt: time between frames.
        """
        self.struc_file = struc_file
        self.traj_file = traj_file
        np.seterr(all='raise') # raise any warning as an error in numpy.
        super().__init__(struc_file, traj_file)
        self.lipnames = [ x for x in set(self.residues.resnames) if "PC" in x ]
        self.dt = ts * self.trajectory.dt

    @property
    def uniq_resname(self):
        """ returns a list of resnames. """
        return list(set(self.residues.resnames))

    def atoms_in_residues(self, resname):
        """ returns a list containing atom names in its input residue.
        """
        assert resname in self.uniq_resname , "Input resname isn't found in this trajectory."
        return np.unique(self.select_atoms('resname %s'%resname).atoms.names)
        

    def __repr__(self):
        resname_string = "<Residues are " + ", ".join(set(self.residues.resnames)) + ">\n"
        sim_info = "<contains {} atoms, {} frames with {} ps timesteps>".format(self.atoms.n_atoms,
                                                                               len(self.trajectory),
                                                                                self.dt)
        return resname_string + sim_info
    
    def get_box_x(self):
        """returns a list contains cell length in x direction"""
        result = []
        for frame in self.trajectory:
            result += [frame.dimensions[0]]
        return result
    
    @property
    def num_lipids(self):
        """ return number of lipids (assuming PC) in one leaflet, assuming bilayer."""
        result = 0
        
        for resname in self.residues.resnames:
            if resname in self.lipnames:
                result += 1
        return round(result/2)

    def bilayer_reshape(self):

        return 
    
    def get_apl(self):
        """ returns the mean apl and stadard error of simulation.
        """
        boxx = self.get_box_x()
        _, n_samples = self.corr_time(boxx)
        apl = np.average(boxx) ** 2 / 100 # in angstrom^2
        apl_se = np.std(boxx)/np.sqrt(n_samples)
        
        return apl, apl_se

    @staticmethod
    def aveNearTwo(array: np.array) -> np.array:
        a = array[:-1]
        b = array[1:]
        return ( a + b ) /2

    @staticmethod
    def repair_bilayer(u, frame):
        """ This function repairs bilayer based on the status of water. If there's only one layer of water,
        bilayer must be broken.
        input:
            1. u: universe of specific frame.
            2. frame: frame of current frame stores dimension.
        no returns, modify coordinates of current frame to make bilayer whole.
        """
        water_coors = u.select_atoms('resname TIP3').positions
        dens_, coors_ = np.histogram(water_coors[:,2], bins=30)

        # if there's only one layer of water bulk, create a blank space
        if 0 not in dens_:
            water_coors[:,2][water_coors[:,2] < 0 ] += frame.dimensions[2]  # move water lower than 0 plane one z pbc up.
            dens_, coors_ = np.histogram(water_coors[:,2], bins=30)
            coors_ = u.aveNearTwo(coors_)
            blank_z_coor = np.average(coors_[ dens_ == 0 ])
            u.atoms.translate([0, 0, -blank_z_coor]).wrap()
            
        elif dens_[0] != 0 and dens_[-1] != 0:
            coors_ = u.aveNearTwo(coors_)
            blank_z_coor = np.average(coors_[ dens_ == 0 ])
            u.atoms.translate([0, 0, -blank_z_coor]).wrap()
            # If there are two layers of water bulk, find the center z coor of blank bulk.
            
        elif dens_[0] == 0 or dens_[-1] == 0:
            coors_ = u.aveNearTwo(coors_)
            blank_z_coor = np.average(coors_[ dens_ == 0 ])
            u.atoms.translate([0, 0, -blank_z_coor]).wrap()
                
            
        
    def get_thickness(self):
        """ returns a list of thickness and standard error of thickness. Calculated in terms of phosphorous atoms.
        """
        l1, l1_std, l2, l2_std = self.get_P_z_average()
        thickness = abs(l1 - l2)
        std = thickness * np.sqrt(l1_std ** 2 + l2_std ** 2)
        return [thickness, std]



    #####
    ##        
    ##   Calculate autocorrelation function.
    ##
    #####


    @staticmethod
    def corr_time(array):
        """ 
        This function returns the autocorrelation function of its input.
        input array: a numeric array.
        return: index of 0 point in its input array, numeric.
        
        Algorithms from https://doi.org/10.33011/livecoms.1.1.5067.

        Returns:
        c_f: aotucorrelation function. Can be used for plotting.
        n_ind: number of independent samples. 

        variables:
        _mean, mean of its input 
        _variance, variance of its input
        c_f, average amount of autocorrelation function between snapshots separated by t'.
        _length, length of its input array which is the length of simulation.

        """

        _mean = np.average(array)
        _variance = np.var(array)
        c_f = []
        _length = len(array)

        # looping through its input array for different t'.
        for t_prime in range(1, _length):
            c_f += [ np.average((array[:-t_prime]-_mean)*(array[t_prime:]-_mean)) / _variance]
        
        x_zero, y_zero = Isec355.first_nearest_zero(c_f)
        n_ind = _length//(1+2*np.sum(c_f[:x_zero]))
        
        return c_f, n_ind
     
    @staticmethod
    def first_nearest_zero(array):
        """ returns the index and value of the first point whose y value changes sign in the array (+ -> - ).
        """
    
        for i, value in enumerate(array):
            if array[i] > 0 and array[i+1] < 0:
                return [i, value]
            


    ###
    #
    #  Calculate order parameter of acyl chains.
    #
    ###

    def order_param(self, init_frame = 0, final_frame = -1):
        """ returns the order parameter of two acyl chains.
        inputs init_frame, final_frame: two integers define start and final frame to analyze

        """
        def Sc(v1, v2):
            norm = np.array([0, 0, 1])
            """ returns order parameter of two numpy arrays contain x-y-z coordinates.
            inputs"""
            cos_theta = (np.dot(v1 - v2, norm) / (np.linalg.norm(v1 - v2, axis=2)))
            return 0.5 * (3 * cos_theta **2 - 1)

        c_names_dict = {}
        op_results = {}
        ## Get lipid resnames 
        lip_names = [ x for x in self.uniq_resname if 'PC' in x ]
        #print(lip_names)
        ## Get acyl chain carbon atom names
        for lip_name in lip_names:
            chain1_atom_names = natsort.natsorted(re.findall( r'C2[0-9]+', ','.join(self.atoms_in_residues(lip_name))))
            chain2_atom_names = natsort.natsorted(re.findall( r'C3[0-9]+', ','.join(self.atoms_in_residues(lip_name))))
            c_names_dict[lip_name] = {'c1': chain1_atom_names, 'c2': chain2_atom_names,
                                      'c1_string': "resname %s and "%lip_name + "(name " + " or name ".join(chain1_atom_names) + ")",
                                      'c2_string': "resname %s and "%lip_name + "(name " + " or name ".join(chain2_atom_names) + ")"}
            op_results[lip_name] = [[], []]
        ## Get repaired bilayer.
        for frame in self.trajectory[init_frame:final_frame]:
            try:
                Isec355.repair_bilayer(self, frame) # try to make bilayer whole
            
            except FloatingPointError:
                print("Cannot reshape bilayer, skip frame {} for trajectory {}.".format(frame.frame, self.traj_file))
                continue
            
            ## Loop through lipids
            for lip_name in lip_names:
                ## Get positions of carbons in tails.
                c1_pos = self.select_atoms(c_names_dict[lip_name]['c1_string']).positions
                c2_pos = self.select_atoms(c_names_dict[lip_name]['c2_string']).positions
                ## Handle these positions. there are n - 2 bonds among n carbon atoms. therefore, the size of
                ## order parameter would be n - 2.

                ## Loop through the first acyl
                y1 = self.select_atoms(c_names_dict[lip_name]['c1_string']).positions.reshape(( -1, len(c_names_dict[lip_name]['c1']), 3))
                y2 = self.select_atoms(c_names_dict[lip_name]['c2_string']).positions.reshape(( -1, len(c_names_dict[lip_name]['c2']), 3))

                op_results[lip_name][0].append(np.average(Sc(y1[:,:-2], y1[:,2:]), axis=0))
                op_results[lip_name][1].append(np.average(Sc(y2[:,:-2], y2[:,2:]), axis=0))
                
        for lip_name in lip_names:    
            op_results[lip_name] = np.average(op_results[lip_name], axis=1)

        return c_names_dict, op_results



        

def get_bilipid_pathes(path):
    """ returns folder names contains results of bi-lipid bilayers under its input directory.
    input paht: directory contains result folders.

    This function works as following:
    1. Navigates to its input folders
    2. Searches for folders whose naming creteria is [0-9]{2-3}D[OY]D[DHLM]
    3. Returns a list of folders
    """
    import regex as re
    import natsort

    os.chdir(path)
    result = []
    for path in os.scandir():
        if path.is_dir() and re.match("[0-9]+D[YO][PD][CHDLM]", path.name):
            result += [path.name]

    pathes_for_sorting = [ re.sub('[0-9]+', "", x) + re.match('[0-9]+', x).group(0) for x in result ]
    sorted_index = natsort.index_natsorted(pathes_for_sorting)
    sorted_result = natsort.order_by_index( result, sorted_index)
    return sorted_result

def get_files(ext):
    """
    returns a list of files ends with its input in current folder with an accending order.
    input path: path to target directory. default: current folder
    """
    files = os.listdir()
    extfiles = [ file for file in files if file.endswith(ext) ]
    extfiles = natsort.natsorted(extfiles, alg=natsort.ns.REAL)
    return extfiles

def chunk(array, arr_size=1):
    """ returns a list contains arr_size of sublists of its input array
    """
    result = []
    arr_length = len(array)
    sub_array = range(0, arr_length, int( (arr_length - 1) / arr_size))
    for i in range(len(sub_array) - 1):
        result += [array[ sub_array[i] : sub_array[i+1]]]

    return result

def get_ka(autocorr = True):
    """ returns the compressibility modulus calculated based on apl.
    autocorrelation is used for splitting the apl of whole simulation into several
    parts, then calculate ka of each part.
    """
    ## Create some variables
    all_apl_list = []
    all_n_samples = []
    kas = [] # stores ka of each chunk
    ## Searching for dcd and psf files.
    dcdfiles = get_files('dcd')
    psffiles = get_files('psf')
    if dcdfiles == []:
        raise FileNotFoundError("No dcd files start with gamma found in current folder.")
    elif psffiles == []:
        raise FileNotFoundError("No psf file found in current folder.")
    elif len(psffiles) > 1:
        raise FileExistsError("There are multiple psf files in current folder, be sure to have only one.")
    ## dcd files must starts with gamma
    dcdfiles = [ x for x in dcdfiles if x.startswith('gamma')]
    gammas = list(set([ int(re.findall('-?[0-9]+', x)[0]) for x in dcdfiles ]))
    gammas = natsort.natsorted(gammas, alg=natsort.ns.REAL)
    print("Current composition: ".format())
    print("Surface tensions considered: ".format(", ".join([ str(x) for x in gammas])))
    print("DCD files used: , PSF file used: .".format(", ".join(dcdfiles), psffiles[0]))
    ## loop through results of each surface tension
    for gamma in gammas:
        r = re.compile('gamma'+str(gamma))
        dcdfiles_in_gamma = list( filter( r.match, dcdfiles))
        lobby = Isec355(psffiles[0], dcdfiles_in_gamma)
        # get a list of apl
        apl_list = [ x ** 2 / 100 for x in lobby.get_box_x() ] ## a list of apl in Angstrom^2
        all_apl_list += [apl_list]

        if autocorr:
        # get number of independent frames
            _, n_samples = lobby.corr_time(apl_list)
            all_n_samples += [n_samples]

    if autocorr:
        # number of chunks to split:
        min_n_samples = min(all_n_samples)
        all_apl_list = [ chunk(x, min_n_samples) for x in all_apl_list ]
        ## Calculate mean apl for each chunk.
        all_apl_list = [ np.average(x) for y in all_apl_list for x in y ]
        all_apl_list = np.array(all_apl_list).reshape((len(gammas), -1))
        # Loop each column of different gamma apls and apply linear regression
        for i in range(min_n_samples):
            area_strain = [ y[i] / all_apl_list[0][i] - 1  for y in all_apl_list ]
            ka = stats.linregress(area_strain, gammas)
            kas.append(ka.slope)
        return [np.average(kas), np.std(kas)/np.sqrt(min_n_samples)]
        
    else:
        all_apl_list = [ np.average(x) for x in all_apl_list ]
        area_strain = [ x / all_apl_list[0] - 1 for x in all_apl_list ]
        ka = stats.linregress(area_strain, gammas)
        return [ka.slope, ka.stderr]




if __name__ == "__main__":
    download_folder = "/Users/jiamingxu/Downloads"
    folders = get_bilipid_pathes(download_folder)
    print(folders)
    for folder in folders:
        os.chdir(folder)
        ka = get_ka(True)
        print(folder, ka)
        os.chdir('../')
