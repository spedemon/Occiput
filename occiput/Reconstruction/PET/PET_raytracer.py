
# occiput.io 
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# 2013 - 2015
# Boston, MA


# Interface to the ray-tracer for projection and back-projection. 
# The interface to the ray-tracer for scatter simulation is in another file. 


__all__ = ["DEFAULT_PROJECTION_PARAMETERS","DEFAULT_BACKPROJECTION_PARAMETERS","ProjectionParameters","BackprojectionParameters", 
           "PET_project_compressed","PET_backproject_compressed","has_NiftyPy"]

from occiput.Core.NiftyPy_wrap import PET_project_compressed, PET_backproject_compressed, has_NiftyPy


DEFAULT_PROJECTION_PARAMETERS  =  {
                         "N_samples":                256, 
                         "sample_step":              0.8, 
                         "background_activity":      0.0, 
                         "background_attenuation":   0.0, 
                         "truncate_negative_values": 0,  
                         "gpu_acceleration":         1,      
                         "direction":                7,      # affects performance only, between 1 and 6; 2 and 5 are normally the best value
                         "block_size":               512 }   # affects performance only

DEFAULT_BACKPROJECTION_PARAMETERS  =  {
                         "N_samples":                256, 
                         "sample_step":              0.8, 
                         "background_activity":      0.0, 
                         "background_attenuation":   0.0, 
                         "truncate_negative_values": 0,  
                         "gpu_acceleration":         1,      
                         "direction":                7,      # affects performance only, between 1 and 6; 2 and 5 are normally the best value
                         "block_size":               512 }   # affects performance only



class ProjectionParameters(): 
    """Data structure containing the parameters of a projector for PET. """
    default_parameters = DEFAULT_PROJECTION_PARAMETERS
    def __init__(self,parameters=None): 
        self._initialised = 0
        self.name = "Unknown binning name"
        if parameters is None: 
            self.load_from_dictionary(self.default_parameters)  
        elif type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        elif type(parameters) in [list, tuple]: 
            if len(parameters) == len(self.default_parameters.keys()): 
                self.N_samples                = parameters[0] 
                self.sample_step              = parameters[1]
                self.background_activity      = parameters[2]
                self.background_attenuation   = parameters[3] 
                self.truncate_negative_values = parameters[4]
                self.gpu_acceleration         = parameters[5]
                self.direction                = parameters[6]
                self.block_size               = parameters[7]
            else: 
                raise UnknownParameter('Parameter %s specified for ProjectionParameters is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for ProjectionParameters is not compatible. '%str(parameters)) 
            
    def load_from_dictionary(self,dictionary):
        self.N_samples                = dictionary['N_samples']                 # Number of samples along a line when computing line integrals
        self.sample_step              = dictionary['sample_step']               # distance between consecutive points along a line when computing line integrals (this is in the same unit measure as the size of the imaging volume (activity_size and attenuation_size))
        self.background_activity      = dictionary['background_activity']       # Activity in voxels outside of the imaging volume
        self.background_attenuation   = dictionary['background_attenuation']    # Attenuation in voxels outside of the imaging volume 
        self.truncate_negative_values = dictionary['truncate_negative_values']  # If set to 1, eventual negative values obtained when projecting are set to 0
                                                                                # (This is meant to remove eventual unwanted small negative values due to FFT-based smoothing within the projection algorithm 
                                                                                # - note that numerical errors may produce small negative numbers when doing FFT-IFFT even if the function is all positive )
        self.gpu_acceleration         = dictionary['gpu_acceleration']          # Whether to use GPU acceleration (if available) or not 
        self.direction                = dictionary['direction']                 # Direction parameter for the projector and b-proj., it affects performance only; 2 and 5 are normally the best values. 
        self.block_size               = dictionary['block_size']                # Number of blocks per thread for the GPU accelerated projector and b-proj. It affects performance only. Typical values: 256, 512, 768
        


class BackprojectionParameters(ProjectionParameters): 
    """Data structure containing the parameters of a projector for PET. """
    default_parameters = DEFAULT_BACKPROJECTION_PARAMETERS
    # Note: this seems unfinished, but it inherits all that is necessary from ProjectionParameters


      
