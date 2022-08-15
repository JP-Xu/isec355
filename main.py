import numpy as np
import os
import natsort
import MDAnalysis

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
        super().__init__(struc_file, traj_file)
        self.lipnames = [ x for x in set(self.residues.resnames) if "PC" in x ]
        self.dt = ts * self.trajectory.dt
        
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
        """ returns a list of apl and stadard error of apl.
        """
        boxx = self.get_box_x()
        _, n_samples = Isec355.corr_time(boxx)
        apl = np.average(boxx) ** 2 / 100 # in angstrom^2
        apl_se = np.std(boxx)/np.sqrt(n_samples)
        
        return apl, apl_se
        
    def get_thickness(self):
        """ returns a list of thickness and standard error of thickness.
        """

        
    ## Calculate autocorrelation function.
    @staticmethod
    def corr_time(array):
        """ 
        This function returns the autocorrelation function of its input.
        input array: a numeric array.
        return: index of 0 point in its input array, numeric.
        
        Returns:
        c_f: aotucorrelation function. Can be used for plotting.
        n_samples: number of samples of simulation trajectory calculated by autocorrelation.

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
        n_samples = round(_length/x_zero)
        
        return c_f, n_samples
     
    @staticmethod
    def first_nearest_zero(array):
        """ returns the index and value of the first point whose y value changes sign in the array (+ -> - ).
        """
    
        for i, value in enumerate(array):
            if array[i] > 0 and array[i+1] < 0:
                return [i, value]
            


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
        if path.is_dir() and re.match("[0-9]+D[YO]D[HDLM]", path.name):
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

    




if __name__ == "__main__":
    download_folder = "/Users/jiamingxu/Downloads"
    folders = get_bilipid_pathes(download_folder)
    print(folders)
    os.chdir('65DODH')
    dcdfiles = get_files('dcd')
    print(dcdfiles)
    os.chdir('../')