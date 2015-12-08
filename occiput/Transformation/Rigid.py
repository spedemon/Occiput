
# occiput   
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging
# Jan 2014, Boston, MA, USA 


__all__ = ['RigidTransformationSSD']

from ilang.Models import Model 
from ilang.Graphs import ProbabilisticGraphicalModel 
from ilang.Samplers import Sampler

from occiput.Core import Image3D 
from occiput.Visualization import MultipleVolumes

try: 
    from NiftyPy.NiftyReg import resample_image_rigid
    from NiftyPy.NiftyReg import deriv_intensity_wrt_space_rigid
    from NiftyPy.NiftyReg import deriv_intensity_wrt_transformation_rigid
    from NiftyPy.NiftyReg import deriv_ssd_wrt_transformation_rigid
    from NiftyPy.NiftyReg import gaussian_smoothing
except: 
    has_NiftyPy = False
    print "Please install NiftyPy"
else: 
    has_NiftyPy = True 

import numpy





class ModelRigidSSD(Model): 
    variables          = {'source':'continuous','target':'continuous','transformation':'continuous'} 
    dependencies       = [['source','target','directed'],['transformation','target','directed']] 
    preferred_samplers = {'transformation':['QuasiNewton_L_BFGS_B']}

    # init
    def init(self): 
        pass 
         
    # expose to sampler: 
    def log_conditional_probability_transformation(self,T): 
        source    = self.get_value('source')
        target    = self.get_value('target')
        #return -.5*numpy.dot(numpy.dot((x-mu),hessian),(x-mu).T)
        FIXME
        return 0

    def log_conditional_probability_gradient_transformation(self,T): 
        source    = self.get_value('source')
        target    = self.get_value('target')      
        #return -.5*numpy.dot((x-mu),hessian+hessian.T) 
        # FIXME
        return numpy.random.rand(6)*4/30000






class RigidTransformationSSD(): 
    def __init__(self,target,source): 
        self.target = target 
        self.source = source 
        self._resampled_source_cache = None 
        self._resampled_source_need_update = True
        self._contruct_ilang_model()
        self._set_transformation(numpy.zeros((1,6))) 
        
    def _contruct_ilang_model(self): 
        self.ilang_model   = ModelRigidSSD('RigidSSD')
        self.graph         = ProbabilisticGraphicalModel(['target','source','transformation']) 
        self.graph.set_nodes_given(['target','source'],True) 
        self.graph.set_node_value('source',self.source.get_data())
        self.graph.set_node_value('target',self.target.get_data())
        self.graph.set_node_value('transformation',numpy.zeros((1,6)))
        self.graph.add_dependence(self.ilang_model,{'target':'target','source':'source','transformation':'transformation'}) 
        self.sampler       = Sampler(self.graph) 

    def _set_transformation(self,transformation):
        self.graph.set_node_value('transformation',transformation) 
        self._resampled_source_need_update = True 
    
    def _get_transformation(self): 
        return self.graph.get_node_value('transformation')
      
    def estimate_transformation(self,method=None,iterations=30000,trace=True,parameters=None, display_progress=True): 
        if method!=None: 
            self.sampler.set_node_sampling_method_manual('transformation',method)
        else: 
            self.sampler.set_node_sampling_method_auto('transformation')
        last_sample = self.sampler.sample_node('transformation',nsamples=iterations,trace=trace) 
        # This is not necessary, but it triggers the cache of the resampled source image: 
        self._set_transformation(last_sample) 
        return last_sample 

    def resample_in_target_space(self,image): 
        # FIXME 
        #print "Resampling .."
        transformation = self._get_transformation() 
        return image

    def _get_resampled_source(self): 
        #FIXME: remove the next few lines
        if self.transformation.sum() != 0: 
            return self.target
        
        if self._resampled_source_need_update: 
            resampled_source = self.resample_in_target_space(self.source)
            self._resampled_source_cache = resampled_source
            self._resampled_source_need_update = False 
        else: 
            resampled_source = self._resampled_source_cache
        return resampled_source 

    def display_in_browser(self,axis=0,shrink=256,rotate=90,subsample_slices=4,scale_factors=None):
        self.display(axis,shrink,rotate,scale_factors,open_browser=True)
        
    def display(self,axis=0,shrink=256,rotate=90,subsample_slices=4,scale_factors=None,open_browser=False): 
        resampled = self.resampled_source
        D = MultipleVolumes([resampled,self.target], axis, shrink, rotate, subsample_slices, scale_factors, open_browser)
        return D.display() 

    def _repr_html_(self): 
        return self.display()._repr_html_()

    transformation   = property(_get_transformation, _set_transformation, None)
    resampled_source = property(_get_resampled_source, None, None)

