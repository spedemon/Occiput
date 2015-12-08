
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Apr. 2014, Boston, MA


import numpy
from .ilang_models import SSD_ilang
from ilang.Graphs import ProbabilisticGraphicalModel 
from ilang.Samplers import Sampler 
from occiput.Visualization import MultipleVolumesNiftyPy
from occiput.Core import Image3D, Transform_Translation 


DEFAULT_N_POINTS = [100,100,100]
DEFAULT_N_ITER   = 30
DEFAULT_SIGMA    = 3000000


class Registration_Two_Images(object): 
    def __init__( self, source=None, target=None, degrees_of_freedom=3, sigma=DEFAULT_SIGMA, initial_transformation=None ): 
        self.set_source(source)
        self.set_target(target)  
        self.set_sigma(sigma)
        self.__make_graph() 
        if initial_transformation == None: 
            self.set_transformation( numpy.zeros(degrees_of_freedom) )
        else: 
            self.set_transformation( initial_transformation ) 

    def __make_graph(self): 
        self.ilang_model = SSD_ilang() 
        # The following is an alternative use of ilang models: define the log_p and gradient functions here in order to have access to 
        # the local variables. 
        self.ilang_model.log_conditional_probability_transformation = self.__P
        self.ilang_model.log_conditional_probability_gradient_transformation = self.__G

        self.ilang_graph = ProbabilisticGraphicalModel(['source','target','transformation','sigma'])
        self.ilang_graph.set_nodes_given(['sigma','source','target'],True) 
        self.ilang_graph.add_dependence(self.ilang_model,{'source':'source','target':'target','sigma':'sigma','transformation':'transformation'}) 
        self.ilang_sampler = Sampler(self.ilang_graph)

    def __set_space(self): 
        # make sure that the images are in the same space
        if self.source.space != self.target.space: 
            raise "Source [%s] and Target [%s] must be in the same space"%(source.space,target.space)
        self.space = self.source.space 

    def set_transformation(self,tr):  
        self.__transformation = tr

    def __get_transformation(self): 
        return self.__transformation 

    def set_cost_function(self,cost): 
        pass 

    def __initialize_registration(self, optimization_method='QuasiNewton_L_BFGS_B', n_points=DEFAULT_N_POINTS): 
        self.ilang_graph.set_node_value('source',self.source.data)
        self.ilang_graph.set_node_value('target',self.target.data)
        self.ilang_graph.set_node_value('sigma',self.sigma)

        if n_points == None: 
            n_points = self.target.shape
        self.grid = self.target.get_world_grid(n_points) 
        self.resampled_target = self.target.compute_resample_on_grid(self.grid) 

        self.ilang_sampler.set_node_sampling_method_manual('transformation',optimization_method)

    def register(self, optimization_method='GradientAscent',iterations=DEFAULT_N_ITER, n_points=DEFAULT_N_POINTS):  
        self.__initialize_registration(optimization_method, n_points)
        self.ilang_graph.set_node_value('transformation',self.transformation)
        self.ilang_sampler.sample(iterations) 
        self.transformation = self.ilang_graph.get_node_value('transformation')

    def set_source(self,source): 
        self.source = source
        if hasattr(self,'target'): 
            self.__set_space()
        
    def set_target(self,target): 
        self.target = target
        if hasattr(self,'source'): 
            self.__set_space()

    def set_sigma(self,sigma):
        self.sigma = sigma

    def get_result(self): 
        result = self.source.copy()
        result.transform( self.get_transformation_matrix() )
        return result
    
    def get_transformation_matrix(self):
        #rot = self.__transformation[0:3]
        #tra = self.__transformation[3:6]    
        tra = self.__transformation
        return Transform_Translation(tra, self.space, self.space)    #FIXME: rotation and translation 

    def display(self): 
        D = MultipleVolumesNiftyPy([self.source, self.target, self.get_result()])
        #D = MultipleVolumes([self.source, self.target])
        return D

    def __G(self,transformation): 
#        print "Transformation: ", transformation 
        #rot = transformation[0:3]
        #tra = transformation[3:6]
        tra = transformation
        T = Transform_Translation(tra, self.space, self.space)    #FIXME: rotation and translation 
        re_s = self.source.compute_resample_on_grid(self.grid, affine_grid_to_world=T)  #FIXME: memoize
        gr_s = self.source.compute_gradient_on_grid(self.grid, affine_grid_to_world=T)
        re_t = self.resampled_target 
        G_tra0 = (re_t.data-re_t.data*gr_s[0]).sum()
        G_tra1 = (re_t.data-re_t.data*gr_s[1]).sum()
        G_tra2 = (re_t.data-re_t.data*gr_s[2]).sum()
        G_tra = [G_tra0, G_tra1, G_tra2]
        G = numpy.asarray([G_tra[0], G_tra[1], G_tra[2]]) / self.sigma
#        print "gradient:       ", G
#        print "log_likelihood: ", self.__P(transformation)
        return G

    def __P(self,transformation): 
        # FIXME: memoize (now image is resampled to compute log_p and the gradient) 
        #print transformation 
        #rot = transformation[0:3]
        #tra = transformation[3:6]
        tra  = transformation
        T = Transform_Translation(tra, self.space, self.space)    #FIXME: rotation and translation 
        resampled_source = self.source.compute_resample_on_grid(self.grid, affine_grid_to_world=T) 
        P = -numpy.linalg.norm(resampled_source.data - self.resampled_target.data) / self.sigma  #FIXME: verify
#        print "log_likelihood: ", P
        return P

    def _repr_html_(self):
        #return self.ilang_graph._repr_html_()
        return self.display()._repr_html_()

    def __repr__(self): 
        return "Registration of 2 images."

    transformation = property(__get_transformation, set_transformation)



class Registration_Longitudinal(): 
    def __init__(self, images ): 
        self.__images = images 
        self.__make_graph() 
        
    def set_transformation(self,tr): 
        pass 

    def set_cost_function(self,cost): 
        pass 
        
    def register(self): 
        pass 

    def __make_graph(self): 
        self.ilang_graph = ProbabilisticGraphicalModel()
        for i in range(len(self.__images)): 
            self.ilang_graph.add_node('im_%d'%i) 
            self.ilang_graph.set_nodes_given(['im_%d'%i,],True)
        for i in range(len(self.__images)-1): 
            self.ilang_graph.add_nodes(['sigma_%d'%i,'T_%d'%i])
            self.ilang_graph.set_nodes_given(['sigma_%d'%i,],True) 
            model = SSD_ilang('SSD_%d'%i) 
            self.ilang_graph.add_dependence(model,{'sigma':'sigma_%d'%i,'transformation':'T_%d'%i,'source':'im_%d'%i,'target':'im_%d'%(i+1)}) 
        self.ilang_sampler = Sampler(self.ilang_graph)

    def _repr_html_(self):
        return self.ilang_graph._repr_html_()
    
    def __repr__(self): 
        return "Longitudinal registation."






class Registration_N_Images(): 
    def __init__(self, images ): 
        self.__images = images 
        self.__make_graph() 
        
    def set_transformation(self,tr): 
        pass 

    def set_cost_function(self,cost): 
        pass 
        
    def register(self): 
        pass 

    def __make_graph(self): 
        self.ilang_graph = ProbabilisticGraphicalModel()
        self.ilang_graph.add_node('template') 
        for i in range(len(self.__images)): 
            self.ilang_graph.add_node('im_%d'%i) 
            self.ilang_graph.set_nodes_given(['im_%d'%i,],True)
        for i in range(len(self.__images)): 
            self.ilang_graph.add_nodes(['sigma_%d'%i,'T_%d'%i])
            self.ilang_graph.set_nodes_given(['sigma_%d'%i,],True) 
            model = SSD_ilang('SSD_%d'%i) 
            self.ilang_graph.add_dependence(model,{'sigma':'sigma_%d'%i,'transformation':'T_%d'%i,'source':'im_%d'%i,'target':'template'}) 
        self.ilang_sampler = Sampler(self.ilang_graph)

    def _repr_html_(self):
        return self.ilang_graph._repr_html_()
    
    def __repr__(self): 
        return "Registration of N images. "



