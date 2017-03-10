
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Jan. 2015, Boston
# Feb. 2015, Helsinki
# Nov. 2015, Boston 


# Here is the library of PET scanners! Access to listmode data is provided by an external Package. 
# Don't be put off, the mechanism is quite simple. 

__all__ = ["Generic","Brain_PET","Biograph_mMR","Discovery_RX","get_scanner_by_name"]


from occiput.Reconstruction.PET.PET_meshing import Michelogram
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 


# Import scanner definitions 
try: 
    from Occiput_Interface_Brain_PET import Brain_PET
except: 
    Brain_PET = None
try: 
    from Occiput_Interface_Biograph_mMR import Biograph_mMR
except: 
    Biograph_mMR = None
try: 
    from Occiput_Interface_Discovery_RX import Discovery_RX
except: 
    Discovery_RX = None



class Generic(): 
    def __init__(self):
        self.model        = "Generic PET Scanner" 
        self.manufacturer = "Occiput's immagination"
        self.version      = "1.0"
        self.supports_listmode = False
        self.uses_meshing      = False

        self.N_u    = 128
        self.N_v    = 128
        self.size_u = 2.0*128 
        self.size_v = 2.0*128
        self.N_azimuthal                        = 5
        self.N_axial                            = 120 
        self.angles_azimuthal                   = float32([-0.5, -0.25, 0.0, 0.25, 0.5])
        self.angles_axial                       = float32( linspace(0,pi-pi/self.N_axial,self.N_axial) )

        self.scale_activity                     = 1.0

        self.activity_N_samples_projection_DEFAULT       = 150 
        self.activity_N_samples_backprojection_DEFAULT   = 150 
        self.activity_sample_step_projection_DEFAULT     = 2.0 
        self.activity_sample_step_backprojection_DEFAULT = 2.0 
        self.activity_shape_DEFAULT                      = [128,128,128] 
        self.activity_size_DEFAULT                       = float32([2.0, 2.0, 2.0])*float32(self.activity_shape_DEFAULT)

        self.attenuation_N_samples_projection_DEFAULT       = 150 
        self.attenuation_N_samples_backprojection_DEFAULT   = 150 
        self.attenuation_sample_step_projection_DEFAULT     = 2.0 
        self.attenuation_sample_step_backprojection_DEFAULT = 2.0 
        self.attenuation_shape_DEFAULT                      = [128,128,128] 
        self.attenuation_size_DEFAULT                       = float32([2.0, 2.0, 2.0])*float32(self.attenuation_shape_DEFAULT)

        self.listmode = None
        self.physiology = None 


def get_scanner_by_name(name): 
    if name == "Generic": 
        return Generic
    elif name=="Brain_PET" or name=="BrainPET": 
        return Brain_PET
    elif name=="Biograph_mMR" or name=="Siemens_Biograph_mMR" or name=="mMR" or name=="Siemens_mMR": 
        return Biograph_mMR
    elif name=="Discovery_RX" or name=="GE_Discovery_RX" or name=="GE_RX" or name =="RX": 
        return Discovery_RX
    else: 
        return None 


    

    
