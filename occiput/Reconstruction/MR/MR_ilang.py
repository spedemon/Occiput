
# occiput 
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# March 2015, Boston, MA



import ilang 
import ilang.Models 
from ilang.Models import Model 
from ilang.Graphs import ProbabilisticGraphicalModel 

__all__ = ['MR_Static_Gaussian','MR_Dynamic_Gaussian','ProbabilisticGraphicalModel']



class MR_Static_Gaussian(Model): 
    variables = {'x':'continuous','k':'continuous','sigma':'continuous'} 
    dependencies = [['x','k','directed'],['sigma','k','directed']]

    def __init__(self, MR_scan, name=None):
        if name == None:  
            name = self.__class__.__name__
        Model.__init__(self, name) 
        # PET scan
        self.MR_scan = MR_scan    
        # small number
        self.EPS = 1e9

    def set_MR_scan(self, MR_scan): 
        self.MR_scan = MR_scan 

    def init(self): 
        pass 

    def log_conditional_probability_x(self,x_): 
        sigma = self.get_value('sigma')
        return 0 

    def log_conditional_probability_gradient_x(self,x_): 
        sigma = self.get_value('sigma')
        #projection_data = self.PET_scan.project(lambda_, alpha) 
        #Nx = lambda_.shape[0]; Ny = lambda_.shape[1]; Nz = lambda_.shape[2]
        #gradient        = self.PET_scan.backproject(self.PET_scan.get_measurement()[0]/(projection_data[0]+self.EPS), Nx,Ny,Nz, alpha) 
        gradient = 0
        return gradient 

    def log_conditional_probability_sigma(self,sigma): 
        x_ = self.get_value('x') 
        return 0 

    def log_conditional_probability_gradient_sigma(self,sigma): 
        x_ = self.get_value('x')  
        return 0
        
    def sample_conditional_probability_k(self): 
        return 0    



class MR_Dynamic_Gaussian(Model): 
    variables = {'x':'continuous','sigma':'continuous','roi_1':'continuous','roi_1':'continuous','k_1':'continuous','k_2':'continuous'} 
    dependencies = [['x','k_1','directed'],['x','k_2','directed'],['sigma','k_1','directed'],['sigma','k_2','directed'],['roi_1','k_1','directed'],['roi_2','k_2','directed']] 

    def __init__(self, MR_scan, name=None): 
        if name == None:  
            name = self.__class__.__name__
        Model.__init__(self, name) 
        # PET scan object: 
        self.MR_scan = MR_scan 
        self.N_time_bins = len(self.MR_scan.time_bins)
        # Variables and dependencies: 
        self._x      = None 
        self._sigma  = None 
        self._rois   = None 
        self.variables = {'x':'continuous', 'sigma':'continuous'} 
        self.dependencies = []
        for t in range(self.N_time_bins): 
            var_name_counts = 'k_'+str(t)
            var_name_roi    = 'roi_'+str(t)
            self.variables[var_name_counts]='continuous'
            self.variables[var_name_roi]='continuous'
            self.dependencies.append(['x',var_name_counts,'directed'])
            self.dependencies.append(['sigma',var_name_counts,'directed']) 
            self.dependencies.append([var_name_roi,var_name_counts,'directed']) 

    def set_MR_scan(self, MR_scan): 
        self.MR_scan = MR_scan 

    def init(self): 
        pass 

    def log_conditional_probability_lambda(self): 
        return 0 

    def log_conditional_probability_gradient_lambda(self): 
        return numpy.zeros([100,1])

    def log_conditional_probability_alpha(self): 
        return 0 

    def log_conditional_probability_gradient_alpha(self): 
        return numpy.zeros([100,1])
        
    def sample_conditional_probability_z(self): 
        return 0    





        