
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Apr. 2014, Boston, MA

import numpy
from occiput.Core import Image3D, Grid3D, Transform_6DOF
from ilang.Models import Model 




class SSD_ilang( Model ):
    variables = {'source':'continuous','target':'continuous','transformation':'continuous','sigma':'continuous'} 
    dependencies = [['source','target','directed'],['transformation','target','directed'],['sigma','target','directed']]

    def __init__( self, name=None ):
        if name == None:  
            name = self.__class__.__name__
        Model.__init__(self, name) 
    
    def log_conditional_probability_transformation(self,transformation): 
        source = self.get_value('source') 
        target = self.get_value('target') 
        sigma  = self.get_value('sigma')
        log_p = 0.0
        return log_p

    def log_conditional_probability_gradient_transformation(self,transformation): 
        source = self.get_value('source') 
        target = self.get_value('target') 
        sigma  = self.get_value('sigma')
        gradient = numpy.asarray([0.0,0.0,0.0,0.0,0.0,0.0]) 
        return gradient 
        
    def sample_conditional_probability_target(self): 
        return 0    

    def init(self): 
        pass 











