
# occiput.io 
# Stefano Pedemonte
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Aalto University, Finland 
# 2013-2015


# Code to generate subsets for OSEM (ordered subsets expectation maximization) and algorithms of that sort, such 
# as stochastic optimization. 


__all__ = ["SubsetGenerator"]

from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 
from numpy.random import randint 



class SubsetGenerator():  
    """This object has the machinery to generate subsets of the axial and longitudinal directions for PET reconstruction (e.g. OSEM algorithm). 
    It implements various algorithms to create the subsets. """
    def __init__(self,N_azimuthal,N_axial):
        self._N_axial     = N_axial
        self._N_azimu     = N_azimuthal
        self._index = 0

    def new_subset(self, mode, subset_size, azimuthal_range=None): 
        """Returns a new subset. """
        if mode=='random': 
            return self._subsets_random_no_replacement(subset_size,azimuthal_range) 
        elif mode=='random_axial': 
            return self._subsets_random_axial(subset_size)
        elif mode=='ordered_axial': 
            return self._subsets_ordered_axial(subset_size)
        else: 
            raise UnexpectedParameter("'mode' parameter %s not recognised."%str(mode))
        
    def all_active(self): 
        """Returns full set.  """
        return ones((self._N_azimu, self._N_axial),dtype=uint32,order="C") 

    def _subsets_random_no_replacement(self,subset_size,azimuthal_range=None): 
        """Generates subsets randomly - no replacement. """
        if subset_size is None: 
            return self.all_active() 
        if subset_size>=self._N_axial*self._N_azimu: 
            return self.all_active() 
        M = zeros((self._N_azimu, self._N_axial),dtype=uint32,order="C") 
        n = 0
        while n<subset_size:
            active_axial = randint(self._N_axial)
            if azimuthal_range is None: 
                active_azimu = randint(self._N_azimu) 
            else: 
                active_azimu = azimuthal_range[ randint(len(azimuthal_range)) ]
            if M[active_azimu,active_axial] == 0: 
                M[active_azimu,active_axial] = 1
                n+=1
        return M

    def _subsets_ordered_axial(self, subset_size): 
        """Generates ordered subsets; use all azimuthal angles, subsample axially. """
        if subset_size is None: 
            return self.all_active() 
        if subset_size>=self._N_axial*self._N_azimu: 
            return self.all_active() 
        M = zeros((self._N_azimu, self._N_axial),dtype=uint32,order="C") 
        for i in range(self._index, self._N_axial, self._N_axial/subset_size): 
            M[:,i] = 1
        self._index += 1 
        if self._index == self._N_axial: 
            self._index = 0 
        return M
        
    def _subsets_random_axial(self, subset_size): 
        """Generates random subsets; use all azimuthal angles, subsample axially. """
        if subset_size is None: 
            return self.all_active() 
        if subset_size>=self._N_axial*self._N_azimu: 
            return self.all_active() 
        M = zeros((self._N_azimu, self._N_axial),dtype=uint32,order="C") 
        n = 0
        while n < subset_size:
            active_axial = randint(self._N_axial)
            if M[0,active_axial] == 0: 
                M[:,active_axial] = 1
                n+=1
        return M

    
