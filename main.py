import numpy as np


class isec355:
    """
    An extension class for MDAnalysis.
    This class contains several functions:
    length of side of x-y plane.
    """
    
    def __init__(self, Universe, ts=2):
        """ 
        input Universe: MDAnalysis Universe.
        """
        self.u = Universe
        self.lipnames = [ x for x in set(self.u.residues.resnames) if "PC" in x ]
        self.dt = ts * Universe.trajectory.dt
        
    def __repr__(self):
        resname_string = "<Residues are " + ", ".join(set(self.u.residues.resnames)) + ">\n"
        sim_info = "<contains {} atoms, {} frames with {} ps timesteps>".format(self.u.atoms.n_atoms,
                                                                               len(self.u.trajectory),
                                                                                self.dt)
        return resname_string + sim_info
    
    def get_box_x(self):
        """returns a list contains cell length in x direction"""
        result = []
        for frame in self.u.trajectory:
            result += [frame.dimensions[0]]
        return result
    
    def get_num_lipids(self):
        """ return number of lipids (assuming PC) in one leaflet, assuming bilayer."""
        result = 0
        
        for resname in self.u.residues.resnames:
            if resname in self.lipnames:
                result += 1
        return round(result/2)
    
    def get_apl(self):
        """ returns apl and stadard error of apl.
        """
        boxx = self.get_box_x()
        _, n_samples = isec355.corr_time(boxx)
        apl = np.average(boxx) ** 2 / 100 # in angstrom^2
        apl_se = np.std(boxx)/np.sqrt(n_samples)
        
        return apl, apl_se
        
        
        
    ## Calculate autocorrelation function.
    
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
        
        x_zero, y_zero = isec355.first_nearest_zero(c_f)
        n_samples = round(_length/x_zero)
        
        return c_f, n_samples
                          
    def first_nearest_zero(array):
        """ returns the index and value of the first point whose y value changes sign in the array (+ -> - ).
        """
    
        for i, value in enumerate(array):
            if array[i] > 0 and array[i+1] < 0:
                return [i, value]
            


