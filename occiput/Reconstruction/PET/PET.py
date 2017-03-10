
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Jan. 2015, Boston
# June 2015, Helsinki 
# Nov. 2015 - Mar. 2017, Boston 


# If you are looking for PET reconstruction, this is where to start. 
# The objects defined here provide a abstractions for Static and Dynamic PET reconstruction, 
# abstracting the scanner geometries and vendor models and providing an interface to the 
# software tools for projection, backprojection and reconstruction. 


__all__ = ['PET_Static_Scan', 'PET_Multi2D_Scan', 'PET_Dynamic_Scan', 'PET_Cyclic_Scan', 'Binning', 'PET_Projection_Sparsity', 'PET_Projection', 'RigidTransform'] 


# Set verbose level 
# This is a global setting for occiput. There are 3 levels of verbose: high, low, no_printing 
from occiput.global_settings import * 
#set_verbose_high() 
#set_verbose_low() 
set_verbose_no_printing() 


# Import occiput: 
from occiput.Core import Image3D
from occiput.Core.Errors import FileNotFound, UnknownParameter, UnexpectedParameter
from occiput.Core.Print import millisec_to_min_sec, pretty_print_large_number, print_percentage
from occiput.Core.Print import rad_to_deg, deg_to_rad, array_to_string
from occiput.Core import Transform_Identity, Transform_6DOF, Transform_Affine, Transform_Scale

from occiput.Reconstruction.PET.PET_projection import PET_Projection, Binning, PET_Projection_Sparsity, PET_compress_projection
from occiput.Reconstruction.PET.PET_projection import PET_uncompress_projection, PET_initialize_compression_structure
from occiput.Reconstruction.PET.PET_projection import display_PET_Projection_geometry
from occiput.Reconstruction.PET.PET_raytracer import ProjectionParameters, BackprojectionParameters
from occiput.Reconstruction.PET.PET_raytracer import PET_project_compressed, PET_backproject_compressed
from occiput.Reconstruction.PET.PET_subsets import SubsetGenerator 

from occiput.Reconstruction.PET.PET_scanners import Generic, get_scanner_by_name
from occiput.DataSources.Synthetic.Shapes import uniform_cylinder
from occiput.DataSources.FileSources.vNAV import load_vnav_mprage
from occiput.DataSources.FileSources import import_nifti
from occiput.DataSources.FileSources import import_interfile_volume, export_interfile_volume 
from occiput.DataSources.FileSources import import_interfile_projection, export_interfile_projection
from occiput.DataSources.FileSources import guess_file_type_by_name 
from occiput.DataSources.FileSources import import_PET_Projection 

from occiput.Visualization import *
from occiput.Visualization.Colors import *
from occiput.Visualization import ipy_table, has_ipy_table, svgwrite, has_svgwrite 


# Import ilang (inference language; optimisation) 
from PET_ilang import PET_Static_Poisson, PET_Dynamic_Poisson, ProbabilisticGraphicalModel
from ilang.Samplers import Sampler 


# Import DisplayNode to produce ipython notebook visualisations
from DisplayNode import DisplayNode


# Import interfile data handling module 
from interfile import Interfile


# Import other modules
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, where, ndarray, nan, tile
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 
import os
import h5py 
import time
import copy
try: 
    import pylab
except: 
    has_pylab = False
else: 
    has_pylab = True


# Default parameters 
DEFAULT_SUBSET_SIZE       = 24
DEFAULT_RECON_ITERATIONS  = 10
DEFAULT_N_TIME_BINS       = 15
EPS = 1e-6



def f_continuous(var): 
    """Makes an nd_array Fortran-contiguous. """
    if isinstance(var,ndarray): 
        if not var.flags.f_contiguous: 
            var = asarray(var,order='F')
    else: 
        if hasattr(var,'data'): 
            if isinstance(var.data,ndarray): 
                if not var.data.flags.f_contiguous: 
                    var.data = asarray(var.data,order='F')           
    return var 

        

# FIXME: eliminate the class RigidTransform; use transformation matrices in Image3D instead, for activity and attenuation volumes. # Use Core.Transform_6DOF or Core.Transform_Affine if required, to parameterize the projector and back_projector. 

class RigidTransform(): 
    """Region of Interest. Legacy! """
    def __init__(self,parameters=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):  
        if type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        elif type(parameters) in [list, tuple]: 
            if len(parameters) == 6: 
                self.x       = parameters[0]
                self.y       = parameters[1]
                self.z       = parameters[2]
                self.theta_x = parameters[3]
                self.theta_y = parameters[4]
                self.theta_z = parameters[5]
            else: 
                raise UnknownParameter('Parameter %s specified for the construction of RigidTransform is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for the construction of RigidTransform is not compatible. '%str(parameters)) 
    
    def load_from_dictionary(self,dictionary):
        self.x       = dictionary['x']                           # Translation along x
        self.y       = dictionary['y']                           # Translation along y
        self.z       = dictionary['z']                           # Translation along z
        self.theta_x = dictionary['theta_x']                     # Rotation around x
        self.theta_y = dictionary['theta_y']                     # Rotation around y
        self.theta_z = dictionary['theta_z']                     # Rotation around z

    def __repr__(self):
        s = "PET volume location (RigidTransform): \n"
        s = s+" - x:             %f \n"%self.x
        s = s+" - y:             %f \n"%self.y
        s = s+" - z:             %f \n"%self.z
        s = s+" - theta_x:       %f \n"%self.theta_x
        s = s+" - theta_y:       %f \n"%self.theta_y
        s = s+" - theta_z:       %f \n"%self.theta_z
        return s

    def _repr_html_(self): 
        if not has_ipy_table: 
            return "Please install ipy_table."
        table_data = [['x',self.x],['y',self.y],['z',self.z],['theta_x',self.theta_x],['theta_y',self.theta_y],['theta_z',self.theta_z]] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()









class PET_Static_Scan(): 
    """PET Static Scan. """
    def __init__(self): 
        self.use_gpu(True)                               # by default, use GPU.   
        self.set_scanner(Generic)                        # set scanner geometry and load interface 
        self.activity      = None                        # memoization of activity.  
        self.attenuation   = None                        # memoization of attenuation. 
        self.attenuation_projection = None               # memoization of attenuation projection. 
        self.sensitivity   = None                        # sensitivity is a permanent parameter. 
        self.prompts       = None                        # measurement: prompts. Initialized as empty data structure.
        self.randoms       = None
        self.scatter       = None
        self._normalization = None                        # normalization volume - for all projections - memoize
        self._need_normalization_update = True            # If True, the normalization volume needs to be recomputed 
        self.use_compression(False)

        self.set_transform_scanner_to_world(Transform_Identity(map_from='scanner',map_to='world'))
        
        self._time_profiling_init()
        #self._construct_ilang_model() 
        

    def set_transform_scanner_to_world(self, transform): 
        #FIXME: verify that the transform maps from 'scanner' to 'world' 
        self.transform_scanner_to_world = transform 

    def _make_Image3D_activity(self, data=None): 
        shape = float32(self.activity_shape)
        size  = float32(self.activity_size)
        T_scanner_to_world = self.transform_scanner_to_world 
        T_pix_to_scanner = Transform_Scale(size/shape, map_from='pixels_PET_Static', map_to='scanner') 
        T_pix_to_world = T_scanner_to_world.left_multiply(T_pix_to_scanner)
        image = Image3D(data=data, affine=T_pix_to_world, space='world')
        return image

    def _make_Image3D_attenuation(self, data=None): 
        shape = float32(self.attenuation_shape)
        size  = float32(self.attenuation_size)
        T_scanner_to_world = self.transform_scanner_to_world 
        T_pix_to_scanner = Transform_Scale(size/shape, map_from='pixels_PET_Static', map_to='scanner') 
        T_pix_to_world = T_scanner_to_world.left_multiply(T_pix_to_scanner)
        image = Image3D(data=data, affine=T_pix_to_world, space='world')
        return image

    def time_profile_projection(self):
        return self._timing_projection
    
    def time_profile_backprojection(self):
        return self._timing_backprojection
    
    def time_profile_reconstruction(self): 
        return self._timing_reconstruction

    def _time_profiling_init(self):
        self._timing_reconstruction = {}
        self._timing_projection = {
             'T0_transfer_to_gpu': 0.0,
             'T1_alloc': 0.0,
             'T2_rotate': 0.0,
             'T3_resample': 0.0,
             'T4_integral': 0.0,
             'T5_transfer_to_host': 0.0,
             'T6_free': 0.0,
             't1_project_total': 0.0,
             't2_make_continuous': 0.0,
             't3_prepare_compressed': 0.0,
             't4_scale': 0.0,
             't5_make_projection': 0.0}
        self._timing_backprojection = {
             'T0_transfer_to_gpu': 0.0,
             'T1_alloc': 0.0,
             'T2_rotate': 0.0,
             'T3_resample': 0.0,
             'T4_integral': 0.0,
             'T5_transfer_to_host': 0.0,
             'T6_free': 0.0,
             'T7_accumulate': 0.0,
             'T8_clear_memory': 0.0,
             'T9_copy_texture': 0.0,
             't1_backproject_total': 0.0,
             't2_prepare_compressed': 0.0}
    
    def _time_profiling_reset(self):
        self._timing_reconstruction["T1_project_transfer"] = 0.0
        self._timing_reconstruction["T2_backproject_transfer"] = 0.0
        self._timing_reconstruction["T3_project_alloc"] = 0.0
        self._timing_reconstruction["T4_backproject_alloc"] = 0.0       
        self._timing_reconstruction["T5_project_integral"] = 0.0
        self._timing_reconstruction["T6_backproject_integral"] = 0.0
        self._timing_reconstruction["T7_project_rotate"] = 0.0
        self._timing_reconstruction["T8_backproject_rotate"] = 0.0
        self._timing_reconstruction["T9_project_resample"] = 0.0
        self._timing_reconstruction["T10_backproject_resample"] = 0.0
        self._timing_reconstruction["T11_project_subsets"] = 0.0
        self._timing_reconstruction["T12_backproject_subsets"] = 0.0
        self._timing_reconstruction["T13_compose_projections"] = 0.0
        self._timing_reconstruction["T14_compute_update"] = 0.0

    def _time_profiling_record_projection(self):
        p = self._timing_projection
        self._timing_reconstruction["T1_project_transfer"] += p['T0_transfer_to_gpu']+p['T5_transfer_to_host']
        self._timing_reconstruction["T3_project_alloc"] += p['T1_alloc']+p['T6_free']
        self._timing_reconstruction["T5_project_integral"] += p['T4_integral']
        self._timing_reconstruction["T7_project_rotate"] += p['T2_rotate']
        self._timing_reconstruction["T9_project_resample"] += p['T3_resample']
        self._timing_reconstruction["T11_project_subsets"] += p['t3_prepare_compressed']

    def _time_profiling_record_backprojection(self): 
        b = self._timing_backprojection 
        self._timing_reconstruction["T2_backproject_transfer"] += b['T0_transfer_to_gpu']+ \
            b['T9_copy_texture']+ b['T5_transfer_to_host']
        self._timing_reconstruction["T4_backproject_alloc"] += b['T1_alloc'] + b['T6_free'] 
        self._timing_reconstruction["T6_backproject_integral"] += b['T4_integral']
        self._timing_reconstruction["T8_backproject_rotate"] += b['T2_rotate']
        self._timing_reconstruction["T10_backproject_resample"] += b['T3_resample']
        self._timing_reconstruction["T12_backproject_subsets"] += b['t2_prepare_compressed']
    
    def _time_profiling_record_norm(self):
        pass 
        
    def set_activity_shape(self, activity_shape): 
        if not len(activity_shape) == 3: 
            print "Invalid activity shape"  #FIXME: raise invalid input error
        else: 
            self.activity_shape = activity_shape 

    def set_activity_size(self, activity_size): 
        if not len(activity_size) == 3: 
            print "Invalid activity size"  #FIXME: raise invalid input error
        else: 
            self.activity_size = activity_size 
        self._adapt_line_step_size_activity() 

    def set_attenuation_shape(self, attenuation_shape): 
        if not len(attenuation_shape) == 3: 
            print "Invalid attenuation shape"  #FIXME: raise invalid input error
        else: 
            self.attenuation_shape = attenuation_shape 

    def set_attenuation_size(self, attenuation_size): 
        if not len(attenuation_size) == 3: 
            print "Invalid attenuation size"  #FIXME: raise invalid input error
        else: 
            self.attenuation_size = attenuation_size 
        self._adapt_line_step_size_attenuation() 

    def _adapt_line_step_size_activity(self):   # FIXME: move this calculation in the raytracer 
        if not hasattr(self,'activity_size'): 
            activity_size = float32([0,0,0])
        elif self.activity_size is None: 
            activity_size = float32([0,0,0])
        else: 
            activity_size = float32(self.activity_size)
        diagonal = sqrt((activity_size**2).sum())
        self.activity_projection_parameters.sample_step = diagonal / self.activity_projection_parameters.N_samples 
        self.activity_backprojection_parameters.sample_step = diagonal / self.activity_backprojection_parameters.N_samples 

    def _adapt_line_step_size_attenuation(self):   # FIXME: move this calculation in the raytracer 
        if not hasattr(self,'attenuation_size'): 
            attenuation_size = float32([0,0,0])
        elif self.attenuation_size is None: 
            attenuation_size = float32([0,0,0])
        else: 
            attenuation_size = float32(self.attenuation_size)
        diagonal = sqrt((attenuation_size**2).sum())
        self.attenuation_projection_parameters.sample_step = diagonal / self.attenuation_projection_parameters.N_samples 
        self.attenuation_backprojection_parameters.sample_step = diagonal / self.attenuation_backprojection_parameters.N_samples 

    def _construct_ilang_model(self): 
        # define the ilang probabilistic model 
        self.ilang_model = PET_Static_Poisson(self) 
        # construct a basic Directed Acyclical Graph
        self.ilang_graph = ProbabilisticGraphicalModel(['lambda','alpha','counts']) 
        self.ilang_graph.set_nodes_given(['counts','alpha'],True) 
        self.ilang_graph.add_dependence(self.ilang_model,{'lambda':'lambda','alpha':'alpha','counts':'counts'}) 
        # construct a basic sampler object
        #self.sampler     = Sampler(self.ilang_graph)
 
    def set_binning(self, binning): 
        if isinstance(binning,Binning): 
            self.binning = binning
        else:
            self.binning = Binning(binning)
        self._subsets_generator = SubsetGenerator(self.binning.N_azimuthal,self.binning.N_axial) 
        return self.binning 

    def set_scanner(self,scanner): 
        if type(scanner) is type(""): 
            scanner = get_scanner_by_name(scanner) 
        self.scanner = scanner() 

        self.activity_projection_parameters     = ProjectionParameters()       
        self.activity_backprojection_parameters = BackprojectionParameters() 
        self.activity_projection_parameters.N_samples     = self.scanner.activity_N_samples_projection_DEFAULT
        self.activity_projection_parameters.sample_step   = self.scanner.activity_sample_step_projection_DEFAULT
        self.activity_backprojection_parameters.N_samples = self.scanner.activity_N_samples_backprojection_DEFAULT 
        self.activity_backprojection_parameters.sample_step = self.scanner.activity_sample_step_backprojection_DEFAULT

        self.set_activity_shape(self.scanner.activity_shape_DEFAULT)
        self.set_activity_size(self.scanner.activity_size_DEFAULT)
        
        self.activity_projection_parameters.gpu_acceleration = self._use_gpu
        self.activity_backprojection_parameters.gpu_acceleration = self._use_gpu 

        self.attenuation_projection_parameters     = ProjectionParameters()       
        self.attenuation_backprojection_parameters = BackprojectionParameters() 
        self.attenuation_projection_parameters.N_samples     = self.scanner.attenuation_N_samples_projection_DEFAULT
        self.attenuation_projection_parameters.sample_step   = self.scanner.attenuation_sample_step_projection_DEFAULT
        self.attenuation_backprojection_parameters.N_samples = self.scanner.attenuation_N_samples_backprojection_DEFAULT 
        self.attenuation_backprojection_parameters.sample_step = self.scanner.attenuation_sample_step_backprojection_DEFAULT

        self.set_attenuation_shape(self.scanner.attenuation_shape_DEFAULT)
        self.set_attenuation_size(self.scanner.attenuation_size_DEFAULT)

        self.attenuation_projection_parameters.gpu_acceleration = self._use_gpu
        self.attenuation_backprojection_parameters.gpu_acceleration = self._use_gpu 

        binning = Binning()
        binning.size_u = self.scanner.size_u
        binning.size_v = self.scanner.size_v
        binning.N_u    = self.scanner.N_u
        binning.N_v    = self.scanner.N_v
        binning.N_axial          = self.scanner.N_axial
        binning.N_azimuthal      = self.scanner.N_azimuthal 
        binning.angles_axial     = self.scanner.angles_axial
        binning.angles_azimuthal = self.scanner.angles_azimuthal
        self.binning = binning 
        
        self.set_scale_activity(self.scanner.scale_activity) 
        
        self._subsets_generator = SubsetGenerator(self.binning.N_azimuthal,self.binning.N_axial) 
   

    def use_gpu(self, use_it): 
        self._use_gpu = use_it

    def use_compression(self, use_it): 
        self._use_compression = use_it 
        if not use_it: 
            if self.prompts is not None: 
                if self.prompts.is_compressed(): 
                    self.set_prompts( self.prompts.uncompress_self() )
            if self.randoms is not None: 
                if self.randoms.is_compressed(): 
                    self.set_randoms( self.randoms.uncompress_self() )   
            if self.sensitivity is not None: 
                if self.sensitivity.is_compressed(): 
                    self.set_sensitivity( self.sensitivity.uncompress_self() )
            if self.scatter is not None: 
                if self.scatter.is_compressed(): 
                    self.set_scatter( self.scatter.uncompress_self() ) 
        else: 
            if hasattr(self,"_use_compression"): 
                if self._use_compression is False and use_it is True: 
                    # FIXME 
                    #print "Not able to compress once uncompressed. Please implement PET_Projection.uncompress_self() to "
                    #print "enable this functionality. "
                    return 
            if self.prompts is not None: 
                if not self.prompts.is_compressed(): 
                    self.set_prompts( self.prompts.compress_self() )
            if self.randoms is not None: 
                if not self.randoms.is_compressed(): 
                    self.set_randoms( self.randoms.compress_self() )   
            if self.sensitivity is not None: 
                if not self.sensitivity.is_compressed(): 
                    self.set_sensitivity( self.sensitivity.compress_self() )
            if self.scatter is not None: 
                if not self.scatter.is_compressed(): 
                    self.set_scatter( self.scatter.compress_self() ) 

    def get_scatter(self): 
        return self.scatter 

    def set_scatter(self, scatter, duration_ms=None): 
        self.scatter = scatter  
        if duration_ms is not None: 
            self.scatter.time_bins = int32([0,duration_ms])

    def get_prompts(self): 
        return self.prompts
               
    def set_prompts(self,prompts): 
        if isinstance(prompts, PET_Projection): 
            self.prompts = prompts 
            self.sparsity = self.prompts.sparsity     #update self.sparsity (self.sparsity exists to store sparsity
                                                      #information in case there is no prompts data)
#            self.set_binning(prompts.get_binning()) #FIXME: check if it is compatible with the scanner
        elif self.prompts is not None:  
                prompts = PET_Projection( self.prompts.get_binning(), prompts, self.prompts.sparsity.offsets,
                                          self.prompts.sparsity.locations, self.prompts.get_time_bins()) 
                self.prompts = prompts 
        else: 
            print "Prompts data should be an instance of PET_Projection or an array whose dimension"
            print "matches the sparsity pattern of the current projection data. " 
            # FIXME: raise input error and to a try-except when creating the instance of PET_Projection 

    def get_randoms(self): 
        return self.randoms

    def set_randoms(self,randoms): 
        if isinstance(randoms, PET_Projection): 
            self.randoms = randoms
            self.sparsity_delay = self.randoms.sparsity    #update self.sparsity (self.sparsity exists to store 
                                                           #sparsity information in case there is not randoms data)
            #self.set_binning(randoms.get_binning())   #FIXME: make sure binning is consistent with randoms 
        elif self.randoms is not None:  
                randoms = PET_Projection( self.randoms.get_binning(), randoms, self.randoms.sparsity.offsets, self.randoms.sparsity.locations, self.randoms.get_time_bins()) 
                self.randoms = randoms 
        else: 
            print "Delay randoms data should be an instance of PET_Projection or an array whose dimension"
            print "matches the sparsity pattern of the current projection data. " 
            # FIXME: raise input error and to a try-except when creating the instance of PET_Projection 

    def get_sensitivity(self): 
        return self.sensitivity 

    def set_sensitivity(self, sensitivity): 
        #FIXME: verify type: PET_projection or nd_array (the latter only in full sampling mode)
        self.sensitivity = sensitivity 


    def export_sensitivity(self, filename): 
        if self.sensitivity is None: 
            print "Sensitivity has not been loaded"
        else: 
            self.get_sensitivity().save_to_file(filename)

    def export_prompts(self,filename):
        self.get_prompts().save_to_file(filename) 

    def export_randoms(self,filename):
        self.get_randoms().save_to_file(filename) 

    def export_scatter(self,filename): 
        self.get_randoms().save_to_file(filename)


    # FIXME: when importing, compress if compression is enabled 
    def import_sensitivity(self, filename, datafile='', vmin = 0.00, vmax = 1e10):  
        filetype = guess_file_type_by_name(filename)
        if filetype is "h5": 
            sensitivity = import_PET_Projection(filename)
        elif filetype is "interfile_projection_header": 
            sensitivity = import_interfile_projection(filename,self.binning,self.scanner.michelogram,datafile,True,vmin,vmax)
#            if self.prompts is not None:  # FIXME: sensitivity loaded from interfile with some manufacturers has non-zero value 
#                                          # where there are no detectors - set to zero where data is zero 
#                                          #(good approx only for long acquisitions). See if there is a better 
#                                          # way to handle this. 
#                sensitivity.data[self.prompts.data==0]=0
#            else: 
#                print "Warning: If loading real scanner data, please load prompts before loading the sensitivity. Ignore this message if this is a simulation. See the source code for more info. " # FIXME: see comment two lines up 
        elif filetype is "mat": 
            print "Sensitivity from Matlab not yet implemented. All is ready, please spend 15 minutes and implement. "
            return 
        else: 
            print "File type unknown. "
            return 
        sensitivity.data = float32(sensitivity.data)
        if self._use_compression is False: 
            sensitivity = sensitivity.uncompress_self() 
        self.set_sensitivity(sensitivity)

    def import_scatter(self, filename, datafile='', duration_ms=None): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
            projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_scatter: file type unknown. "
            return 
        projection.data = float32(projection.data)
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_scatter(projection, duration_ms) 

    def import_randoms(self, filename, datafile=''): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
             projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_randoms: file type unknown. "
            return
        projection.data = float32(projection.data)
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_randoms(projection)

    def import_prompts(self, filename, datafile=''):
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
            projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile, load_time=True)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_prompts: file type unknown. "
            return
        projection.data = float32(projection.data)
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_prompts(projection)

    def import_attenuation(self, filename, datafile='', filename_hardware='', datafile_hardware=''): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_volume_header": 
            volume = import_interfile_volume(filename, datafile)
        elif filetype is "nifti": 
            print "Nifti attenuation file not supported. Everything is ready to implement this, please implement it. "
            # FIXME: if nifti files are used, sum the hardware image using resampling in the common space
        elif filetype is "h5": 
            print "H5 attenuation file not supported. Everything is ready to implement this, please implement it. "
            # FIXME: if h5 files are used, sum the hardware image using resampling in the common space
        elif filetype is "mat": 
            print "Matlab attenuation file not supported. Everything is ready to implement this, please implement it. "
        else: 
            print "PET.import_attenuation: file type of %s unknown. Unable to load attenuation tomogram. "%filename
            return 
        if filename_hardware is not '': 
            filetype = guess_file_type_by_name(filename_hardware) 
            if filetype is "interfile_volume_header": 
                volume_hardware = import_interfile_volume(filename_hardware, datafile_hardware)
            else: 
                print "File type of %s unknown. Unable to load hardware attenuation tomogram. "%filename_hardware
            volume.data = volume.data + volume_hardware.data  
        volume.data = float32(volume.data)
        self.set_attenuation(volume)

    def import_attenuation_projection(self, filename, datafile=''): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
            projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile, load_time=True)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_attenuation_projection: file type unknown. "
            return
        projection.data = float32(projection.data)
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_attenuation_projection(projection)

    def export_attenuation_projection(self,filename): 
        self.get_attenuation_projection().save_to_file(filename)

    def get_attenuation_projection(self): 
        return self.attenuation_projection

    def set_attenuation_projection(self, attenuation_projection): 
        self.attenuation_projection = attenuation_projection
    
    def import_listmode(self, filename, datafile=None, time_range_ms=[0,None], display_progress=False ): 
        """Load measurement data from a listmode file. """
        print_debug("- Loading static PET data from listmode file "+str(filename) )
        hdr = Interfile.load(filename) 

        # 2) Guess the path of the listmode data file, if not specified or mis-specified; 
        #  1 - see if the specified listmode data file exists 
        if datafile is not None: 
            datafile = datafile.replace("/",os.path.sep).replace("\\",os.path.sep)          # cross platform compatibility 
            if not os.path.exists(datafile): 
                raise FileNotFound("listmode data",datafile)  
        #  2 - if the listmode data file is not specified, try with the name (and full path) contained in the listmode header
        datafile      = hdr['name of data file']['value']
        datafile = datafile.replace("/",os.path.sep).replace("\\",os.path.sep)              # cross platform compatibility
        if not os.path.exists(datafile): 
        #  3 - if it doesn't exist, look in the same path as the header file for the listmode 
        #      data file with name specified in the listmode header file 
            datafile = os.path.split(filename)[0]+os.path.sep+os.path.split(datafile)[-1]  
            if not os.path.exists(datafile): 
        #  4 - if it doesn't exist, look in the same path as the header file for the listmode data 
        #      file with same name as the listmode header file, replacing the extension: ".l.hdr -> .l" 
                if filename.endswith(".l.hdr"): 
                    datafile = filename.replace(".l.hdr",".l") 
                    if not os.path.exists(datafile):     
                        raise FileNotFound("listmode data",datafile)  
        #  5 - if it doesn't exist, look in the same path as the header file for the listmode data 
        #      file with same name as the listmode header file, replacing the extension: ".hdr -> .l" 
                elif filename.endswith(".hdr"): 
                    datafile = filename.replace(".hdr",".l") 
                    if not os.path.exists(datafile):     
                        raise FileNotFound("listmode data",datafile)  
           
        # 3) Determine duration of the acquisition 
        n_packets              = hdr['total listmode word counts']['value'] 
        scan_duration          = hdr['image duration']['value']*1000            # milliseconds
        
        # 4) determine scanner parameters
        n_radial_bins          = hdr['number of projections']['value'] 
        n_angles               = hdr['number of views']['value'] 
        n_rings                = hdr['number of rings']['value'] 
        max_ring_diff          = hdr['maximum ring difference']['value']
        n_sinograms            = n_rings+2*n_rings*max_ring_diff-max_ring_diff**2-max_ring_diff

        # Determine the time binning  
        time_range_0 = time_range_ms[0]
        if time_range_ms[1] is not None: 
            time_range_1 = time_range_ms[1]
        else: 
            time_range_1 = scan_duration  
        time_bins = int32(linspace(time_range_0,time_range_1,2))

        # Display information 
        print_debug(" - Number of packets:    %d       "%n_packets )
        print_debug(" - Scan duration:        %d [sec] "%(scan_duration/1000.0) )
        print_debug(" - Listmode data file:   %s       "%datafile )
        print_debug(" - Listmode header file: %s       "%filename)
        print_debug(" - Number of time bins:  %d       "%(len(time_bins)-1) )
        print_debug(" - Time start:           %f [sec] "%(time_bins[ 0]/1000.0) )
        print_debug(" - Time end:             %f [sec] "%(time_bins[-1]/1000.0) )
        print_debug(" - time_bins:            %s       "%str(time_bins))
        print_debug(" - n_radial_bins:        %d       "%n_radial_bins)
        print_debug(" - n_angles:             %d       "%n_angles)
        print_debug(" - n_angles:             %d       "%n_sinograms)

        if display_progress: 
            progress_bar = ProgressBar()
            progress_callback = progress_bar.set_percentage
        else: 
            def progress_callback(value): 
                pass 
        
        # Load the listmode data 
        M = self.scanner.michelogram

        R = self.scanner.listmode.load_listmode(datafile, n_packets, time_bins, self.binning,n_radial_bins, n_angles, 
                                       n_sinograms, M.span, M.segments_sizes, M.michelogram_sinogram, 
                                       M.michelogram_plane, progress_callback) 
        progress_callback(100) 

        # Load static measurement data 
        self._load_static_measurement() 
        
        # Free structures listmode data
        self.scanner.listmode.free_memory()

        # Construct ilang model 
        self._construct_ilang_model()
    
    def quick_inspect(self, index_axial=0, index_azimuthal=5, index_bin=60): 
        if self.randoms is not None and not isscalar(self.randoms): 
            randoms = self.randoms.to_nd_array()[index_axial,index_azimuthal,:,index_bin]
        else: 
            randoms = 0.0 
        if self.prompts is not None and not isscalar(self.prompts): 
            prompts = self.prompts.to_nd_array()[index_axial,index_azimuthal,:,index_bin]
        else: 
            prompts = 0.0 
        if self.sensitivity is not None and not isscalar(self.sensitivity): 
            sensitivity = self.sensitivity.to_nd_array()[index_axial,index_azimuthal,:,index_bin]
        else: 
            sensitivity = 1.0 
        if self.scatter is not None and not isscalar(self.scatter): 
            scatter = self.scatter.to_nd_array()[index_axial,index_azimuthal,:,index_bin]
            if self.scatter.get_duration() is not None: 
                if self.scatter.get_duration() > 1e-6: 
                    if self.prompts.get_duration() is not None: 
                        if self.prompts.get_duration() > 1e-6: 
                            scatter = scatter * self.prompts.get_duration() / self.scatter.get_duration()
        else: 
            scatter = 0.0 
        if has_pylab: 
            pylab.plot(prompts - randoms)
            pylab.hold(1)
            pylab.plot(sensitivity*scatter,'g') 
        else: 
            print "quick_inspect uses Pylab to display imaging data. Please install Pylab. " 

    def _get_sparsity(self): 
        #if self.prompts == None: 
        # in this case returns sparsity pattern for uncompressed projection 
        sparsity = PET_Projection_Sparsity(self.binning.N_axial, self.binning.N_azimuthal, self.binning.N_u, self.binning.N_v)
        #else: 
        #    sparsity = self.prompts.sparsity
        return sparsity 

    def set_scale_activity(self,scale): 
        self.scale_activity = scale 


    def project_attenuation(self, attenuation=None, unit='inv_cm', transformation=None, sparsity=None, \
                            subsets_matrix=None, exponentiate=True): 
        if attenuation is None: 
            attenuation = self.attenuation 
        if isinstance(attenuation,ndarray): 
            attenuation_data = f_continuous(float32(attenuation)) 
        else: 
            attenuation_data = f_continuous(float32(attenuation.data)) 

        if not list(attenuation_data.shape) == list(self.attenuation_shape): 
            raise UnexpectedParameter("Attenuation must have the same shape as self.attenuation_shape") 

        # By default, the center of the imaging volume is at the center of the scanner 
        tx = 0.5*(self.attenuation_size[0] - self.attenuation_size[0] / self.attenuation_shape[0])
        ty = 0.5*(self.attenuation_size[1] - self.attenuation_size[1] / self.attenuation_shape[1])
        tz = 0.5*(self.attenuation_size[2] - self.attenuation_size[2] / self.attenuation_shape[2])
        if transformation is None:  
            transformation = RigidTransform((tx,ty,tz,0,0,0)) 
        else: 
            transformation = copy.copy(transformation)
            transformation.x = transformation.x + tx
            transformation.y = transformation.y + ty
            transformation.z = transformation.z + tzy

        # Scale according to the unit measure of the specified attenuation. It is assumed that the attenuation map 
        # is constant in a voxel, with the value specified in 'attenuation', of unit measure 'unit'. 
        if unit == 'inv_mm': 
            invert = False 
            scale  = 1.0
        elif unit == 'inv_cm': 
            invert = False 
            scale  = 10.0 
        elif unit == 'mm': 
            invert = True 
            scale  = 1.0 
        elif unit == 'cm': 
            invert = True 
            scale  = 10.0 
        else: 
            print "Unit measure unknown. Assuming inv_cm. Keep track of the unit measures! "
            invert = False
            scale  = 10.0
            
        if invert: 
            attenuation_data = 1.0 / (attenuation_data + eps) 
        step_size_mm = self.attenuation_projection_parameters.sample_step 
        step_size = step_size_mm / scale

        # Optionally project with a sparsity pattern not equal to sparsity associated to the loaded prompts data 
        # Note: if prompts have not been loaded, self._get_sparsity() assumes no compression. 
        
        if sparsity is None: 
            sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            sparsity_subset = sparsity
            angles = self.binning.get_angles()
        else: 
            sparsity_subset = sparsity.get_subset(subsets_matrix)
            angles = self.binning.get_angles(subsets_matrix)
 
        offsets = sparsity_subset.offsets
        locations = sparsity_subset.locations 
        activations = ones([angles.shape[1], angles.shape[2]],dtype="uint32")
        
        # Call the raytracer 
        projection_data, timing = PET_project_compressed(attenuation_data,None,offsets, locations, activations, 
            angles.shape[2], angles.shape[1], angles, 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            0.0,0.0,0.0, 
            transformation.x, transformation.y, transformation.z, transformation.theta_x, transformation.theta_y, transformation.theta_z, 
            0.0,0.0,0.0,0.0,0.0,0.0, 
            self.attenuation_projection_parameters.gpu_acceleration, self.attenuation_projection_parameters.N_samples,
            self.attenuation_projection_parameters.sample_step, self.attenuation_projection_parameters.background_attenuation, 
            0.0, self.attenuation_projection_parameters.truncate_negative_values, 
            self.attenuation_projection_parameters.direction, self.attenuation_projection_parameters.block_size) 
        
        self._timing_projection = timing
        
        # Fix scale and exponentiate 
        if exponentiate: 
            projection_data = exp(-projection_data * step_size) 
        else: 
            projection_data = projection_data * step_size

        # Create object PET_Projection: it contains the raw projection data and the description of the projection geometry 
        # and sparsity pattern. 
        time_bins = int32([0,0])  # Projection of the attenuation does not have timing information 
        projection = PET_Projection( self.binning, projection_data, sparsity.offsets, sparsity.locations, time_bins, subsets_matrix)  
        self.set_attenuation_projection(projection)
        return projection 


    def backproject_attenuation(self, projection, unit="inv_cm", transformation=None, sparsity=None, subsets_matrix=None): 
        if isinstance(projection,ndarray): 
            projection_data = float32(projection)
        else: 
            projection_data = float32(projection.data)

        # By default, the center of the imaging volume is at the center of the scanner 
        tx = 0.5*(self.attenuation_size[0] - self.attenuation_size[0] / self.attenuation_shape[0])
        ty = 0.5*(self.attenuation_size[1] - self.attenuation_size[1] / self.attenuation_shape[1])
        tz = 0.5*(self.attenuation_size[2] - self.attenuation_size[2] / self.attenuation_shape[2])
        if transformation is None:  
            transformation = RigidTransform((tx,ty,tz,0,0,0)) 
        else: 
            transformation = copy.copy(transformation)
            transformation.x = transformation.x + tx
            transformation.y = transformation.y + ty
            transformation.z = transformation.z + tz

        # Scale according to the unit measure of the specified attenuation. It is assumed that the attenuation map 
        # is constant in a voxel, with the value specified in 'attenuation', of unit measure 'unit'. 
        if unit == 'inv_mm': 
            invert = False 
            scale  = 1.0
        elif unit == 'inv_cm': 
            invert = False 
            scale  = 10.0 
        elif unit == 'mm': 
            invert = True 
            scale  = 1.0 
        elif unit == 'cm': 
            invert = True 
            scale  = 10.0 
        else: 
            print "Unit measure unknown. Assuming inv_cm. Keep track of the unit measures! "
            invert = False
            scale  = 10.0
            
        if invert: 
            attenuation_data = 1.0 / (attenuation_data + eps) 
        step_size_mm = self.attenuation_projection_parameters.sample_step 
        step_size = step_size_mm / scale

        if sparsity is None: 
                sparsity = self._get_sparsity() 
        
        if isinstance(projection,ndarray): 
            projection_data = float32(projection)
            offsets = sparsity.offsets
            locations = sparsity.locations
            angles = self.binning.get_angles(subsets_matrix)
            activations = ones([self.binning.N_azimuthal,self.binning.N_axial],dtype="uint32")
        else: 
            projection_data = float32(projection.data)
            offsets = projection.sparsity.offsets
            locations = projection.sparsity.locations 
            angles = projection.get_angles()
            activations = ones([projection.sparsity.N_azimuthal,projection.sparsity.N_axial],dtype="uint32")        
        
        #print "backproject attenuation" 
        #print "projection_data",projection_data.shape
        #print "offsets",offsets.shape
        #print "locations",locations.shape
        #print "activations",activations.shape
        #print "angles",angles.shape
        
        # Call ray-tracer 
        backprojection_data, timing = PET_backproject_compressed(projection_data,None, offsets, locations, activations, 
            angles.shape[2], angles.shape[1], angles, 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.attenuation_shape[0], self.attenuation_shape[1], self.attenuation_shape[2], 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            0.0, 0.0, 0.0, 
            transformation.x, transformation.y, transformation.z, transformation.theta_x, transformation.theta_y, transformation.theta_z, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            self.attenuation_backprojection_parameters.gpu_acceleration, self.attenuation_backprojection_parameters.N_samples,
            self.attenuation_backprojection_parameters.sample_step, 
            self.attenuation_backprojection_parameters.background_attenuation, 0.0, 
            self.attenuation_backprojection_parameters.direction, self.attenuation_backprojection_parameters.block_size) 

        self._timing_backprojection = timing
        backprojection_data = backprojection_data * step_size 
        
        # Set the correct scale - unit measure and return Image3D - FIXME: set scale for requested unit measure 
        return self._make_Image3D_attenuation(backprojection_data)  


    def project_activity(self, activity, unit="Bq/mm3", transformation=None, sparsity=None, subsets_matrix=None):
        t0 = time.time()
        if isinstance(activity,ndarray): 
            activity_data = f_continuous(float32(activity))
        else: 
            activity_data = f_continuous(float32(activity.data))
        
        t1_make_continuous = (time.time()-t0)*1000
        t0 = time.time()
        
        # By default, the center of the imaging volume is at the center of the scanner; no rotation 
        tx = 0.5*(self.activity_size[0] - self.activity_size[0] / self.activity_shape[0])
        ty = 0.5*(self.activity_size[1] - self.activity_size[1] / self.activity_shape[1])
        tz = 0.5*(self.activity_size[2] - self.activity_size[2] / self.activity_shape[2])
        if transformation is None:  
            transformation = RigidTransform((tx,ty,tz,0,0,0)) 
        else: 
            transformation = copy.copy(transformation)
            transformation.x = transformation.x + tx
            transformation.y = transformation.y + ty
            transformation.z = transformation.z + tz
            
        # Optionally project with a sparsity pattern not equal to sparsity associated to the loaded prompts data 
        # Note: if prompts have not been loaded, self._get_sparsity() assumes no compression. 
        if sparsity is None: 
            sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            sparsity_subset = sparsity
            angles = self.binning.get_angles()
        else: 
            sparsity_subset = sparsity.get_subset(subsets_matrix)
            angles = self.binning.get_angles(subsets_matrix)
        
        scale = 1.0  #FIXME: change this according to the input unit measure - check how this is done in project_attenuation
        step_size_mm = self.activity_projection_parameters.sample_step 
        step_size = step_size_mm / scale

        offsets = sparsity_subset.offsets
        locations = sparsity_subset.locations 
        activations = ones([angles.shape[1],angles.shape[2]],dtype="uint32")
        
        #print locations[:,0:20]
        #print locations.flags
        #print sparsity.locations[:,0:20]
        #print sparsity.locations.flags
        
        #print "project activity"
        #print "activity",activity_data.shape
        #print "offsets",offsets.shape
        #print "locations",locations.shape
        #print "activations",activations.shape
        #print "angles.shape",angles.shape

        t2_prepare_compressed = (time.time()-t0)*1000
        t0 = time.time()
        
        # Call the raytracer 
        projection_data, timing = PET_project_compressed(activity_data, None, offsets, locations, activations, 
            angles.shape[2], angles.shape[1], angles,  
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            0.0,0.0,0.0, 
            transformation.x, transformation.y, transformation.z, transformation.theta_x, transformation.theta_y, transformation.theta_z, 
            0.0,0.0,0.0,0.0,0.0,0.0, 
            self.activity_projection_parameters.gpu_acceleration, self.activity_projection_parameters.N_samples,
            self.activity_projection_parameters.sample_step, self.activity_projection_parameters.background_activity, 
            0.0, self.activity_projection_parameters.truncate_negative_values, 
            self.activity_projection_parameters.direction, self.activity_projection_parameters.block_size) 

        t3_project = (time.time()-t0)*1000
        t0 = time.time()
        
        self._timing_projection = timing
        # Create object PET_Projection: it contains the raw projection data and the description of the projection geometry 
        # and sparsity pattern. 
        time_bins = int32([0,1000.0])     # 1 second - projection returns a rate - by design 
        projection_data = projection_data * step_size
        
        t4_scale = (time.time()-t0)*1000
        t0 = time.time()
        
        projection = PET_Projection( self.binning, projection_data, sparsity.offsets, sparsity.locations, time_bins, subsets_matrix)  

        t5_make_projection = (time.time()-t0)*1000

        # Optionally scale by sensitivity, attenuation, time, global sensitivity 
        #if attenuation is not None: 
        #    projection.data = projection.data * attenuation.data #attenuation.compress_as(projection).data
        #if self.sensitivity is not None: 
        #    projection.data = projection.data * self.sensitivity.data.reshape(projection_data.shape)
        #self.sensitivity.compress_as(projection).data

        self._timing_projection["t1_project_total"] = t3_project
        self._timing_projection["t2_make_continuous"] = t1_make_continuous
        self._timing_projection["t3_prepare_compressed"] = t2_prepare_compressed
        self._timing_projection["t4_scale"] = t4_scale
        self._timing_projection["t5_make_projection"] = t5_make_projection
        
        return projection 


    def backproject_activity(self, projection, transformation=None, subsets_matrix=None): 
        # By default, the center of the imaging volume is at the center of the scanner 
        t0 = time.time()
        tx = 0.5*(self.activity_size[0] - self.activity_size[0] / self.activity_shape[0])
        ty = 0.5*(self.activity_size[1] - self.activity_size[1] / self.activity_shape[1])
        tz = 0.5*(self.activity_size[2] - self.activity_size[2] / self.activity_shape[2])
        if transformation is None:  
            transformation = RigidTransform((tx,ty,tz,0,0,0)) 
        else: 
            transformation = copy.copy(transformation)
            transformation.x = transformation.x + tx
            transformation.y = transformation.y + ty
            transformation.z = transformation.z + tz
            
        if not isinstance(projection,ndarray): 
            if not subsets_matrix is None: 
                projection_subset = projection.get_subset(subsets_matrix)
                sparsity_subset = projection_subset.sparsity
                angles = projection_subset.get_angles()  
                projection_data = float32(projection_subset.data)
            else: 
                projection_subset = projection
                sparsity_subset = projection.sparsity
                angles = projection_subset.get_angles()  
                projection_data = float32(projection_subset.data)
        else: 
            sparsity = self._get_sparsity() 
            if not subsets_matrix is None: 
                # doesn't look right 
                indexes = subsets_matrix.flatten() == 1
                projection_data = float32(projection.swapaxes(0,1). \
                                          reshape((sparsity.N_axial*sparsity.N_azimuthal, \
                                                   self.binning.N_u,self.binning.N_v))[indexes,:,:])
                sparsity_subset = sparsity.get_subset(subsets_matrix)
                angles = self.binning.get_angles(subsets_matrix)  
            else: 
                sparsity_subset = sparsity
                angles = self.binning.get_angles()  
                projection_data = float32(projection)
                
        offsets = sparsity_subset.offsets
        locations = sparsity_subset.locations
        activations = ones([sparsity_subset.N_azimuthal,sparsity_subset.N_axial],dtype="uint32")         
            
        scale = 1.0  #FIXME: change this according to the input unit measure - check how this is done in project_attenuation
        step_size_mm = self.activity_projection_parameters.sample_step 
        step_size = step_size_mm / scale  

        t2_prepare_compressed = (time.time()-t0)*1000
        t0 = time.time()
        
        #print "backproject activity" 
        #print "projection_data",projection_data.shape
        #print "offsets",offsets.shape
        #print "locations",locations.shape
        #print "activations",activations.shape
        #print "angles",angles.shape
        #print angles[:,0,0:5]
        #print offsets[0:3,0:5]
        #print locations[:,0:5]
        #time.sleep(0.2)
        
        # Call ray-tracer 
        backprojection_data, timing = PET_backproject_compressed(projection_data, None, offsets, locations, activations, 
            angles.shape[2], angles.shape[1], angles,
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_shape[0], self.activity_shape[1], self.activity_shape[2], 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            0.0, 0.0, 0.0, 
            transformation.x, transformation.y, transformation.z, transformation.theta_x, transformation.theta_y, transformation.theta_z, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            self.activity_backprojection_parameters.gpu_acceleration, self.activity_backprojection_parameters.N_samples,
            self.activity_backprojection_parameters.sample_step, 
            self.activity_backprojection_parameters.background_activity, 0.0, 
            self.activity_backprojection_parameters.direction, self.activity_backprojection_parameters.block_size) 

        t1_backproject_total = (time.time()-t0)*1000
        self._timing_backprojection = timing
        self._timing_backprojection["t1_backproject_total"] = t1_backproject_total
        self._timing_backprojection["t2_prepare_compressed"] = t2_prepare_compressed 
        
        backprojection_data = backprojection_data * step_size
        return self._make_Image3D_activity(backprojection_data)  


    def _load_static_measurement(self, time_bin=None):  
        if time_bin is None: 
            Rp = self.scanner.listmode.get_measurement_static_prompt() 
            Rd = self.scanner.listmode.get_measurement_static_delay() 
        else: 
            Rp = self.scanner.listmode.get_measurement_prompt(time_bin) 
            Rd = self.scanner.listmode.get_measurement_delay(time_bin) 
        time_start        = Rp['time_start'] 
        time_end          = Rp['time_end'] 

        time_bins = int32(linspace(time_start, time_end, 2))
        prompts = PET_Projection( self.binning, Rp['counts'], Rp['offsets'], Rp['locations'], time_bins) 
        randoms = PET_Projection( self.binning, Rd['counts'], Rd['offsets'], Rd['locations'], time_bins) 
        if self._use_compression: 
            self.set_prompts( prompts )
            self.set_randoms( randoms )
        else:
            #print "Uncompressing"
            self.set_prompts( prompts.uncompress_self() )
            self.set_randoms( randoms.uncompress_self() )
        self._construct_ilang_model() 
        

#    def get_activity(self):
#        return self.activity 

#    def set_activity(self,activity): 
#        self.ilang_graph.set_node_value('lambda',activity) 
#        print "PET_Static_Scan.set_activity(): This is for integration with iLang - please implement"
        # FIXME: how about the transformation ? 

    def get_attenuation(self):
        return self.attenuation 

    def set_attenuation(self,attenuation): 
        #self.ilang_graph.set_node_value('alpha',attenuation) 
        # FIXME: how about the transformation ? 
        # FIXME: setting activity and attenuation as members here is not in the spirit of iLang - memoization 
        self.attenuation = attenuation 

    def get_normalization(self, attenuation_times_sensitivity=None, transformation=None, sparsity=None, duration_ms=None, subsets_matrix=None, epsilon=None): 
        # FIXME: memoization 
        if attenuation_times_sensitivity is None: 
            attenuation_times_sensitivity = self.sensitivity  #FIXME: include attenuation here - memoization mumap proj
        if isscalar(attenuation_times_sensitivity) or attenuation_times_sensitivity is None: 
            attenuation_times_sensitivity = ones(self.prompts.data.shape)
        if duration_ms is None: 
            duration_ms = self.prompts.get_duration() 
        duration_sec = duration_ms / 1000.0 
        alpha = self.scale_activity 
        normalization = self.backproject_activity(attenuation_times_sensitivity*duration_sec*alpha,transformation,subsets_matrix) 
        return normalization 

    def get_gradient_activity(self, activity, attenuation=None, unit_activity="Bq/mm3", transformation_activity=None, 
                              sparsity=None, duration_ms=None, subset_size=None, subset_mode='random', subsets_matrix=None,
                              azimuthal_range=None, separate_additive_terms=False, epsilon=None): 
        # Optionally use only a subset of the projections - use, in order: subsets_matrix; subset_size, subset_mode and az_range
        if subsets_matrix is None: 
            if subset_size is not None: 
                if subset_size >= 0: 
                    subsets_matrix = self._subsets_generator.new_subset(subset_mode,subset_size,azimuthal_range) 

        # Optionally use the specified value of epsilon (small number added to the denominator is divisions)
        if epsilon is None: 
            epsilon = EPS 

        if attenuation is None: 
            attenuation = 1.0 
        
        if self.prompts is None: 
            print "self.prompts is None, please set prompts. " 
            return 
            # FIXME : throw an error 

        # By default use the timing information stored in the prompts, however optionally enable overriding 
        if duration_ms is None: 
                duration_ms = self.prompts.get_duration() 
        duration_sec = duration_ms / 1000.0

        # Precompute attenuation*sensitivity - FIXME: do it only is the subset, same for other calculation in proj space
        alpha = self.scale_activity
        prompts = self.prompts 
        randoms = self.randoms
        scatter = self.scatter 
        sensitivity = self.sensitivity 
        if randoms is None: 
            randoms = 0.0
        if scatter is None: 
            scatter = 0.0
        if sensitivity is None or sensitivity is 1.0: 
            att_sens = attenuation
        else: 
            att_sens = sensitivity * attenuation
            
            
        att_sens = att_sens.get_subset(subsets_matrix)
        
        
        # Compute the firt term of the gradient: backprojection of the sensitivity of the scanner
        # If it is requested that the gradient is computed using all the projection measurements, use the 
        #memoized normalization. self.get_normalization() takes care of memoization.  
#        gradient_term1 = self.get_normalization(att_sens, transformation_activity, sparsity, duration_ms, subsets_matrix, epsilon=epsilon)
        norm = PET_Projection(self.binning, data=1.0, subsets_matrix=subsets_matrix)
        gradient_term1 = self.backproject_activity(norm)
        print "the two lines above (1033-1034 PET.py) are temporary"
        self._time_profiling_record_norm()
        
        # Compute the second term of the gradient: backprojection of the ratio between the measurement and the projection of 
        # current activity estimate... Ordinary Poisson to include scatter and randoms.
        projection = self.project_activity(activity, unit=unit_activity, transformation=transformation_activity, sparsity=sparsity,
                                           subsets_matrix=subsets_matrix) 
        self._time_profiling_record_projection()
        prompts_subset = prompts.get_subset(subsets_matrix)
        gradient_term2 = self.backproject_activity( prompts_subset / (projection+randoms/(att_sens*alpha*duration_sec+epsilon)+scatter/(attenuation*alpha*duration_sec+epsilon)+epsilon),
                              transformation=transformation_activity)
        self._time_profiling_record_backprojection()
        if separate_additive_terms: 
            return (gradient_term1, gradient_term2, subsets_matrix) 
        else: 
            gradient = gradient_term1+gradient_term2
            return (gradient, subsets_matrix) 

    def get_gradient_attenuation(self, attenuation, activity, 
                                 sparsity=None, duration_ms=None, subset_size=None, subset_mode='random', subsets_matrix=None,
                                 azimuthal_range=None, epsilon=None): 
        # Optionally use only a subset of the projections - use, in order: subsets_matrix; subset_size, subset_mode and az_range
        if subsets_matrix is None: 
            if subset_size is not None: 
                if subset_size >= 0: 
                    subsets_matrix = self._subsets_generator.new_subset(subset_mode,subset_size,azimuthal_range) 

        # Optionally use the specified value of epsilon (small number added to the denominator is divisions)
        if epsilon is None: 
            epsilon = EPS 

        if attenuation is None: 
            attenuation = 1.0 
        
        if self.prompts is None: 
            print "self.prompts is None, please set prompts. " 
            return 
            # FIXME : throw an error 

        # By default use the timing information stored in the prompts, however optionally enable overriding 
        if duration_ms is None: 
                duration_ms = self.prompts.get_duration() 
        duration_sec = duration_ms / 1000.0

        # Precompute attenuation*sensitivity - FIXME: do it only is the subset, same for other calculation in proj space
        alpha = self.scale_activity
        prompts = self.prompts 
        randoms = self.randoms
        scatter = self.scatter 
        sensitivity = self.sensitivity 
        if randoms is None: 
            randoms = 0.0
        if scatter is None: 
            scatter = 0.0
        if sensitivity is None: 
            sensitivity = 1.0 

        attenuation_projection = self.project_attenuation(attenuation, unit='inv_cm', transformation=None, sparsity=sparsity, subsets_matrix=subsets_matrix, exponentiate=True) 
        
        #FIXME: transformation = None 
        pr_activity = self.project_activity(activity, transformation=None, sparsity=sparsity,
                                           subsets_matrix=subsets_matrix)*sensitivity*attenuation_projection*duration_sec*alpha
        gradient = self.backproject_attenuation(pr_activity - prompts/ (randoms/(pr_activity+epsilon) +  scatter/(pr_activity/(sensitivity+epsilon)+epsilon) + 1 ), unit="inv_cm", transformation=None, sparsity=sparsity, 
                                                subsets_matrix=subsets_matrix)
        return gradient 

    def estimate_activity_and_attenuation(self, activity=None, attenuation=None, iterations=DEFAULT_RECON_ITERATIONS,
                             sparsity=None,
                             subset_size=DEFAULT_SUBSET_SIZE,
                             subset_mode='random', epsilon=None, subsets_matrix=None, azimuthal_range=None, show_progressbar=True): 
        # FIXME: save time: don't compute twice the proj of the attenuation
        activity    = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
        attenuation = self._make_Image3D_attenuation(zeros(self.attenuation_shape,dtype=float32,order="F"))
        if show_progressbar: 
            progress_bar = ProgressBar() 
            progress_bar.set_percentage(0.1) 
        for iteration in range(iterations): 
            activity = self.estimate_activity(activity,attenuation,1,sparsity,subset_size,subset_mode,epsilon,subsets_matrix,azimuthal_range, show_progressbar=False) 
            attenuation = self.estimate_attenuation(activity,attenuation,1,sparsity,subset_size,subset_mode,epsilon,subsets_matrix,azimuthal_range, show_progressbar=False) 
            if show_progressbar: 
                progress_bar.set_percentage((iteration+1)*100.0/iterations) 
        if show_progressbar: 
            progress_bar.set_percentage(100.0) 
        return (activity, attenuation) 

    def estimate_attenuation(self, activity=None, attenuation=None, iterations=DEFAULT_RECON_ITERATIONS, sparsity=None,
                             subset_size=DEFAULT_SUBSET_SIZE,
                             subset_mode='random', epsilon=None, subsets_matrix=None,
                             azimuthal_range=None,show_progressbar=True):
        if show_progressbar: 
            progress_bar = ProgressBar() 
            progress_bar.set_percentage(0.1) 

        if attenuation is None: 
            attenuation = self._make_Image3D_attenuation(zeros(self.attenuation_shape,dtype=float32,order="F"))
        for iteration in range(iterations): 
            attenuation = attenuation + self.get_gradient_attenuation(attenuation, activity, 
                                 sparsity, duration_ms=None, subset_size=subset_size, subset_mode=subset_mode,
                                 subsets_matrix=subsets_matrix,
                                 azimuthal_range=azimuthal_range, epsilon=epsilon) 
            if show_progressbar: 
                progress_bar.set_percentage((iteration+1)*100.0/iterations) 
        if show_progressbar: 
            progress_bar.set_percentage(100.0) 
        return attenuation
    
    def estimate_activity(self, activity=None, attenuation=None, iterations=DEFAULT_RECON_ITERATIONS, sparsity=None,
                          subset_size=DEFAULT_SUBSET_SIZE, subset_mode='random', epsilon=None, subsets_matrix=None,
                          azimuthal_range=None, show_progressbar=True): 
        # Optionally use the specified value of epsilon (small number added to the denominator is divisions)
        if epsilon is None: 
            epsilon = EPS 

        if self.prompts is None: 
            print "self.prompts is None, please set prompts. " 
            return 
            # FIXME : throw an error 

        duration_ms = self.prompts.get_duration() 
        
        if show_progressbar: 
            progress_bar = ProgressBar() 
            progress_bar.set_percentage(0.1) 
        
        #print "Projection of the attenuation. "
        if attenuation is None:
            attenuation = self.attenuation 
        if attenuation is not None: 
            self.attenuation_projection = self.project_attenuation(attenuation) #FIXME: now it's only here that this is defined
        else: 
            self.attenuation_projection = 1.0
        
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
        # FIXME: use transformation - also notice that transformation_activity is always set to None here
    
        self._time_profiling_reset()
        for iteration in range(iterations): 
            [gradient1, gradient2, subsets_matrix] = self.get_gradient_activity(activity,self.attenuation_projection,
                                                                         transformation_activity=None, 
                                                                         sparsity=sparsity,
                                                                         duration_ms=duration_ms, subset_size=subset_size,
                                                                         subset_mode=subset_mode, subsets_matrix=subsets_matrix,
                                                                         azimuthal_range=azimuthal_range,
                                                                         separate_additive_terms=True, epsilon=epsilon) 
            activity = activity * gradient2 / (gradient1+epsilon) 
            if show_progressbar: 
                progress_bar.set_percentage((iteration+1)*100.0/iterations) 
        if show_progressbar: 
            progress_bar.set_percentage(100.0) 
        return activity 

    
    def osem_reconstruction(self, iterations=10, activity=None, attenuation_projection=None, subset_mode="random", \
                            subset_size=64, transformation=None):
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
            subsets_generator = SubsetGenerator(self.binning.N_azimuthal, self.binning.N_axial)
        
        if self.sensitivity is None: 
            sensitivity = self.prompts.copy()
            sensitivity.data = 0.0*sensitivity.data+1
            self.set_sensitivity(sensitivity)
        
        for i in range(iterations):
            print i
            subsets_matrix = subsets_generator.new_subset(subset_mode, subset_size)
            activity = self.osem_step(activity, subsets_matrix, attenuation_projection, transformation)
        return activity
    
    def mlem_reconstruction(self, iterations=10, activity=None, attenuation_projection=None, transformation=None):
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
        
        if self.sensitivity is None: 
            sensitivity = self.prompts.copy()
            sensitivity.data = 0.0*sensitivity.data+1
            self.set_sensitivity(sensitivity)
        
        for i in range(iterations):
            print i
            subsets_matrix = None
            activity = self.osem_step(activity, subsets_matrix, attenuation_projection, transformation)
        return activity
    
    def osem_step(self, activity, subsets_matrix=None, attenuation_projection=None, transformation=None): 
        epsilon = 1e-08
        
        prompts = self.prompts 
        if self._use_compression: 
            prompts = prompts.uncompress_self()
        
        duration_ms = prompts.get_duration() 
        if duration_ms is None: 
            print "Acquisition duration unknown (self.prompts.time_bins undefined); assuming 60 minutes. "
            duration_ms = 1000*60*60 
        duration = duration_ms / 1000.0 
        alpha = self.scale_activity
        
        
        if attenuation_projection is not None: 
            attenuation_projection = attenuation_projection.get_subset(subsets_matrix)
        elif self.attenuation_projection is not None: 
            attenuation_projection = self.attenuation_projection.get_subset(subsets_matrix)
        elif self.attenuation is not None:
            print "Projecting attenuation"
            self.attenuation_projection = self.project_attenuation(self.attenuation)
            attenuation_projection = self.attenuation_projection.get_subset(subsets_matrix)
            print "Done"
        else: 
            attenuation_projection = 1.0
        
        if self.sensitivity is not None: 
            sens_x_att  = self.sensitivity.get_subset(subsets_matrix) * attenuation_projection
        else: 
            sens_x_att   = attenuation_projection
    
        if self.randoms is not None: 
            randoms = self.randoms 
            if self._use_compression: 
                randoms = randoms.uncompress_self()
            randoms = (randoms.get_subset(subsets_matrix) + epsilon) / (sens_x_att * alpha * duration + epsilon) 
        
        if self.scatter is not None: 
            mscatter     = (self.scatter.get_subset(subsets_matrix) + epsilon) / (attenuation_projection * alpha * duration + epsilon)
            # Scale scatter: this is used in dynamic and kinetic imaging, when scatter is calculated using the ativity for a time period longer than the current frame: 
            if self.scatter.get_duration() is not None: 
                if self.scatter.get_duration() > 1e-6: 
                    mscatter = mscatter * duration / self.scatter.get_duration()

        norm = self.backproject_activity(sens_x_att * alpha * duration, transformation=transformation)
    
        projection = self.project_activity(activity, subsets_matrix = subsets_matrix, transformation=transformation)
    
        if self.randoms is not None: 
            if self.scatter is not None:
                update1 = self.backproject_activity(prompts.get_subset(subsets_matrix)/(projection + randoms + mscatter + epsilon), transformation=transformation)
            else: 
                update1 = self.backproject_activity(prompts.get_subset(subsets_matrix)/(projection + randoms + epsilon), transformation=transformation)
    
        else:
            if self.scatter is not None: 
                update1 = self.backproject_activity(prompts.get_subset(subsets_matrix)/(projection + mscatter + epsilon), transformation=transformation)
            else:
                update1 = self.backproject_activity(prompts.get_subset(subsets_matrix)/(projection + epsilon), transformation=transformation)
    
        #activity = activity - gradient1
        activity = (activity /(norm+epsilon)) * update1

        return activity
    

    def brain_crop(self, bin_range=(100,240)):
        if self._use_compression is True: 
            print "Projection cropping currently only works with uncompressed data. "
            print "In order to enable cropping, please complete the implementation of PET_Projection.get_subset()"
            print "Now PET_Projection.get_subset() only works with uncompressed data. "
            return 
        if hasattr(self, "_cropped"): 
            return 
        A = bin_range[0]
        B = bin_range[1]
        self.binning.size_u = (1.0*self.binning.size_u) / self.binning.N_u * (B-A)
        self.binning.N_u = B-A
        if self.prompts is not None: 
            self.prompts = self.prompts.crop((A,B))
        if self.randoms is not None: 
            self.randoms = self.randoms.crop((A,B))          
        if self.scatter is not None: 
            self.scatter = self.scatter.crop((A,B))     
        if self.sensitivity is not None: 
            self.sensitivity = self.sensitivity.crop((A,B))
        self._cropped = True
 
    
    def volume_render(self,volume,scale=1.0): 
        # FIXME: use the VolumeRender object in occiput.Visualization (improve it), the following is a quick fix: 
        [offsets,locations] = PET_initialize_compression_structure(180,1,256,256)
        if isinstance(volume,ndarray): 
            volume = float32(volume)
        else: 
            volume = float32(volume.data)
        subsets_generator = SubsetGenerator(1,180) 
        subsets_matrix=subsets_generator.all_active() 
        mask = uniform_cylinder(volume.shape,volume.shape,
                                [0.5*volume.shape[0],0.5*volume.shape[1],0.5*volume.shape[2]],
                                0.5*min(volume.shape[0]-1,volume.shape[1]),volume.shape[2],2,1,0)
        volume[where(mask.data==0)]=0.0
        direction=7
        block_size=512
        proj, timing = PET_project_compressed(volume, None, offsets,locations, subsets_matrix, 
            180, 1, pi/180, 0.0, 256, 256, 256.0, 256.0, 
            256.0, 256.0, 256.0,  
            256.0, 256.0, 256.0, 
            128.0, 128.0, 128.0, 0.0, 0.0, 0.0,  
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            1, 256, 1.5, 
            0.0, 0.0, 0,
            direction, block_size) 
        proj[where(proj>proj.max()/scale )]=proj.max()/scale
        binning = Binning() 
        binning.N_axial = 180
        binning.N_azimuthal = 1
        binning.angles_axial =     float32( linspace(0,pi-pi/180.0,180) )
        binning.angles_azimuthal = float32( linspace(0,0,1) )
        binning.size_u = 256.0
        binning.size_v = 256.0 
        binning.N_u = 256
        binning.N_v = 256
        projection = PET_Projection(binning, proj, offsets, locations)
        return projection.uncompress_self() 

    def display_geometry(self):
        return display_PET_Projection_geometry()  

    def __repr__(self): 
        s = "Static PET acquisition:  \n" 
        s = s+" - Time_start:                   %s \n"%millisec_to_min_sec(self.prompts.get_time_start() )
        s = s+" - Time_end:                     %s \n"%millisec_to_min_sec(self.prompts.get_time_end() )
        s = s+" - Duration:                     %s \n"%millisec_to_min_sec(self.prompts.get_time_end()-
                                                                           self.prompts.get_time_start())
        s = s+" - N_counts:                     %d \n"%self.prompts.get_integral() 
        s = s+" - N_locations:                  %d \n"%self.prompts.sparsity.get_N_locations()
        #s = s+" - compression_ratio:            %d \n"%self.prompts.sparsity.compression_ratio
        #s = s+" - listmode_loss:                %d \n"%self.prompts.sparsity.listmode_loss
        s = s+" = Scanner: \n"
        s = s+"     - Name:                     %s \n"%self.scanner.model 
        s = s+"     - Manufacturer:             %s \n"%self.scanner.manufacturer 
        s = s+"     - Version:                  %s \n"%self.scanner.version 
        s = s+" * Binning: \n"        
        s = s+"     - N_axial bins:             %d \n"%self.binning.N_axial 
        s = s+"     - N_azimuthal bins:         %d \n"%self.binning.N_azimuthal 
        s = s+"     - Angles axial:             %s \n"%array_to_string(self.binning.angles_axial) 
        s = s+"     - Angles azimuthal:         %s \n"%array_to_string(self.binning.angles_azimuthal)
        s = s+"     - Size_u:                   %f \n"%self.binning.size_u
        s = s+"     - Size_v:                   %f \n"%self.binning.size_v        
        s = s+"     - N_u:                      %s \n"%self.binning.N_u
        s = s+"     - N_v:                      %s \n"%self.binning.N_v
        return s

    def _repr_html_(self):
        if not has_ipy_table: 
            return "Please install ipy_table."
        if self.scanner is not None: 
          table_data = [
          ['Time_start',millisec_to_min_sec(self.prompts.get_time_start())],
          ['Time_end',millisec_to_min_sec(self.prompts.get_time_end())],
          ['Duration',millisec_to_min_sec(self.prompts.get_time_end()-self.prompts.get_time_start())],
          ['N_counts',pretty_print_large_number(self.prompts.get_integral() )],
          ['N_locations',pretty_print_large_number(self.prompts.sparsity.get_N_locations)],
          #['compression_ratio',print_percentage(self.compression_ratio)],
          #['listmode_loss',self.listmode_loss], 
          ['Scanner Name',self.scanner.model],['Scanner Manufacturer',self.scanner.manufacturer],
                      ['Scanner Version',self.scanner.version], ] 
        else:
          table_data = [
          ['Time_start',millisec_to_min_sec(self.prompts.get_time_start())],
          ['Time_end',millisec_to_min_sec(self.prompts.get_time_end())],
          ['Duration',millisec_to_min_sec(self.prompts.get_time_end()-self.prompts.get_time_start())],
          ['N_counts',pretty_print_large_number(self.prompts.get_integral() )],
          ['N_locations',pretty_print_large_number(self.prompts.sparsity.get_N_locations())],     ]
          #['compression_ratio',print_percentage(self.compression_ratio)],
          #['listmode_loss',self.listmode_loss], ] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()
        
        


        

class PET_Multi2D_Scan(PET_Static_Scan):
    def __init__(self):
        self.n_slices = 0 
        self.scatter_duration = None
        self.prompts_duration = None
        self.attenuation_projection = None
        PET_Static_Scan.__init__(self)
        
    def set_activity_size(self, activity_size):
            self.activity_size = ( activity_size[0], activity_size[1], self.n_slices )

    def set_activity_shape(self, activity_shape):
            self.activity_shape = ( activity_shape[0], activity_shape[1], self.n_slices )

    def set_attenuation_projection(self, attenuation_projection):
        self.attenuation_projection = attenuation_projection
            
    def set_prompts_duration(self, duration_array):
        self.prompts_duration = duration_array
        
    def set_scatter_duration(self, duration_array):
        self.scatter_duration = duration_array
            
    def set_number_of_slices(self, n_slices):
        self.n_slices = n_slices
        self.binning.N_azimuthal = 1 
        self.binning.angles_azimuthal = int32([0,])
        self.binning.N_v = self.n_slices + 1
        self.binning.size_v = self.n_slices
        self.set_binning(self.binning)
        if self.scatter_duration is None: 
            self.scatter_duration = ones([self.n_slices+1,])
        if self.prompts_duration is None: 
            self.prompts_duration = ones([self.n_slices+1,])

    def osem_reconstruction(self, iterations=10, activity=None, subset_mode="random", subset_size=64, transformation=None):
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
            subsets_generator = SubsetGenerator(self.binning.N_azimuthal, self.binning.N_axial)
        
        if self.sensitivity is None: 
            sensitivity = self.prompts.copy()
            sensitivity.data = 0.0*sensitivity.data+1
            self.set_sensitivity(sensitivity)
        
        for i in range(iterations):
            print i
            subsets_matrix = subsets_generator.new_subset(subset_mode, subset_size)
            activity = self.osem_step(activity, subsets_matrix, transformation)
        return activity
    
    def mlem_reconstruction(self, iterations=10, activity=None, transformation=None):
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
        
        if self.sensitivity is None: 
            sensitivity = self.prompts.copy()
            sensitivity.data = 0.0*sensitivity.data+1
            self.set_sensitivity(sensitivity)
        
        for i in range(iterations):
            print i
            subsets_matrix = None
            activity = self.osem_step(activity, subsets_matrix, transformation)
        return activity
    
    def osem_step(self, activity, subsets_matrix, transformation): 
        epsilon = 1e-08
        
        prompts = self.prompts 
        if self._use_compression: 
            prompts = prompts.uncompress_self()
        prompts = prompts.get_subset(subsets_matrix)

        # make matrices of duration arrays 
        prompts_duration = tile( self.prompts_duration / 1000.0, (prompts.data.shape[0], prompts.data.shape[1], \
                                                                  prompts.data.shape[2], 1) )
        scatter_duration = tile( self.scatter_duration / 1000.0, (prompts.data.shape[0], prompts.data.shape[1], \
                                                                  prompts.data.shape[2], 1) )
        alpha = self.scale_activity
        
        
        if self.attenuation_projection is not None: 
            attenuation_projection = self.attenuation_projection.get_subset(subsets_matrix)
        else: 
            attenuation_projection = 1.0
        
        if self.sensitivity is not None: 
            sens_x_att  = self.sensitivity.get_subset(subsets_matrix) * attenuation_projection
        else: 
            sens_x_att   = attenuation_projection
    
        if self.randoms is not None: 
            randoms = self.randoms 
            if self._use_compression: 
                randoms = randoms.uncompress_self()
            randoms = (randoms.get_subset(subsets_matrix) + epsilon) / (sens_x_att * alpha * prompts_duration + epsilon) 
        
        if self.scatter is not None: 
            mscatter     = (self.scatter.get_subset(subsets_matrix) + epsilon) / (attenuation_projection * alpha * prompts_duration + epsilon)
            # Scale scatter: this is used in dynamic and kinetic imaging, when scatter is calculated using the ativity for a time period longer than the current frame: 
            mscatter = mscatter * prompts_duration / scatter_duration

        norm = self.backproject_activity(sens_x_att * alpha * prompts_duration, transformation=transformation)
    
        projection = self.project_activity(activity, subsets_matrix = subsets_matrix, transformation=transformation)
    
        if self.randoms is not None: 
            if self.scatter is not None:
                update1 = self.backproject_activity(prompts/(projection + randoms + mscatter + epsilon), transformation=transformation)
            else: 
                update1 = self.backproject_activity(prompts/(projection + randoms + epsilon), transformation=transformation)
    
        else:
            if self.scatter is not None: 
                update1 = self.backproject_activity(prompts/(projection + mscatter + epsilon), transformation=transformation)
            else:
                update1 = self.backproject_activity(prompts/(projection + epsilon), transformation=transformation)
    
        #activity = activity - gradient1
        activity = (activity /(norm+epsilon)) * update1

        return activity    
        
    def quick_inspect(self, index_axial=0, index_slice=0): 
        index_azimuthal = 0
        if self.randoms is not None and not isscalar(self.randoms): 
            randoms = self.randoms.to_nd_array()[index_axial,index_azimuthal,:,index_slice+1]
        else: 
            randoms = 0.0 
        if self.prompts is not None and not isscalar(self.prompts): 
            prompts = self.prompts.to_nd_array()[index_axial,index_azimuthal,:,index_slice+1]
        else: 
            prompts = 0.0 
        if self.sensitivity is not None and not isscalar(self.sensitivity): 
            sensitivity = self.sensitivity.to_nd_array()[index_axial,index_azimuthal,:,index_slice+1]
        else: 
            sensitivity = 1.0 
        if self.scatter is not None and not isscalar(self.scatter): 
            scatter = self.scatter.to_nd_array()
            prompts_duration = tile( self.prompts_duration / 1000.0, (self.prompts.data.shape[0], self.prompts.data.shape[1], \
                                                                  self.prompts.data.shape[2], 1) )
            scatter_duration = tile( self.scatter_duration / 1000.0, (self.prompts.data.shape[0], self.prompts.data.shape[1], \
                                                                  self.prompts.data.shape[2], 1) )
            scatter = scatter * prompts_duration / scatter_duration
            scatter = scatter[index_axial,index_azimuthal,:,index_slice+1]
        else: 
            scatter = 0.0 
        if has_pylab: 
            pylab.plot(prompts - randoms)
            pylab.hold(1)
            pylab.plot(sensitivity*scatter,'g') 
        else: 
            print "quick_inspect uses Pylab to display imaging data. Please install Pylab. "         
        
        
        

        
class PET_Dynamic_Scan(PET_Static_Scan): 
    """PET Dynamic Scan. This is useful for motion correction and for kinetic imaging, or both. """
    def __init__(self): 
        self._dynamic      = []                                 # Sequence of static scans, one per time bin 
        self.time_bins     = []                                 # Time binning 
        self.static        = None
        PET_Static_Scan.__init__(self)

    def import_listmode(self, hdr_filename, time_range_ms=[0,None], data_filename=None, motion_files_path=None, display_progress=False ): 
        """Load prompts data from a listmode file. """
        #Optionally load motion information: 
        if motion_files_path: 
            vNAV = load_vnav_mprage(motion_files_path) 
            self.__motion_events = vNAV 
            if time_range_ms[1] is not None: 
                raise "Either time_bins or motion_files_path should be defined, not both. "
            time_range_ms = self.__motion_events.extract_motion_events().sum()   #FIXME: this is not right, implement binning of sinogram according to motion events (this requires modifying the C code that does the binning and passing the right parameters: the list_mode trigger number corresponding to the beginning and end of each sinogram)

        print_debug("- Loading dynamic PET data from listmode file "+str(hdr_filename) )
        hdr = Interfile.load(hdr_filename) 
        #Extract information from the listmode header

        # 1) Guess the path of the listmode data file, if not specified or mis-specified; 
        #  1 - see if the specified listmode data file exists 
        if data_filename is not None: 
            data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)          # cross platform compatibility 
            if not os.path.exists(data_filename): 
                raise FileNotFound("listmode data",data_filename)  
        #  2 - if the listmode data file is not specified, try with the name (and full path) contained in the listmode header
        data_filename      = hdr['name of data file']['value']
        data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep)              # cross platform compatibility
        if not os.path.exists(data_filename): 
        #  3 - if it doesn't exist, look in the same path as the header file for the listmode data file with name specified in the listmode header file 
            data_filename = os.path.split(hdr_filename)[0]+os.path.sep+os.path.split(data_filename)[-1]  
            if not os.path.exists(data_filename): 
        #  4 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".l.hdr -> .l" 
                if hdr_filename.endswith(".l.hdr"): 
                    data_filename = hdr_filename.replace(".l.hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
        #  5 - if it doesn't exist, look in the same path as the header file for the listmode data file with same name as the listmode header file, replacing the extension: ".hdr -> .l" 
                elif hdr_filename.endswith(".hdr"): 
                    data_filename = hdr_filename.replace(".hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
           
        # 2) Determine the duration of the acquisition 
        n_packets              = hdr['total listmode word counts']['value'] 
        scan_duration          = hdr['image duration']['value']*1000            # milliseconds
        
        # 3) determine scanner parameters
        n_radial_bins          = hdr['number of projections']['value'] 
        n_angles               = hdr['number of views']['value'] 
        n_rings                = hdr['number of rings']['value'] 
        max_ring_diff          = hdr['maximum ring difference']['value']
        n_sinograms            = n_rings+2*n_rings*max_ring_diff-max_ring_diff**2-max_ring_diff

        # Determine the time binning 
        if time_range_ms[1] is None: 
            time_range_ms = int32(linspace(0, scan_duration, DEFAULT_N_TIME_BINS+1))
        elif isscalar(time_range_ms):    #time_bins in this case indicates the number of time bins
            time_range_ms = int32(linspace(0, scan_duration, time_range_ms+1)) 

        # Display information 
        print_debug(" - Number of packets:    %d       "%n_packets )
        print_debug(" - Scan duration:        %d [sec] "%(scan_duration/1000.0) )
        print_debug(" - Listmode data file:   %s       "%data_filename )
        print_debug(" - Listmode header file: %s       "%hdr_filename)
        print_debug(" - Number of time bins:  %d       "%(len(time_range_ms)-1) )
        print_debug(" - Time start:           %f [sec] "%(time_range_ms[ 0]/1000.0) )
        print_debug(" - Time end:             %f [sec] "%(time_range_ms[-1]/1000.0) )
        print_debug(" - time_range_ms:        %s       "%str(time_range_ms))
        print_debug(" - n_radial_bins:        %d       "%n_radial_bins)
        print_debug(" - n_angles:             %d       "%n_angles)
        print_debug(" - n_angles:             %d       "%n_sinograms)

        if display_progress: 
            progress_bar = ProgressBar()
            progress_callback = progress_bar.set_percentage
        else: 
            def progress_callback(value): 
                pass 
        
        # Load the listmode data 
        M = self.scanner.michelogram
        R = self.scanner.listmode.load_listmode(data_filename, n_packets, time_range_ms, self.binning, n_radial_bins, n_angles, n_sinograms, M.span, M.segments_sizes, M.michelogram_sinogram, M.michelogram_plane, progress_callback) 
        if display_progress: 
            progress_bar.set_percentage(100) 

#        self.dynamic_inflation = R['dynamic_inflation']

        N_time_bins  = R['N_time_bins'] 
        time_start   = R['time_start'] 
        time_end     = R['time_end'] 
        #self.time_bins = time_bins[time_start:N_time_bins+1]  #the actual time bins are less than the requested time bins, truncate time_bins 
        self.time_bins = time_range_ms

        # Make list of PET_Static_Scan objects, one per bin
        self._dynamic = [] 
        for t in range(N_time_bins): 
            PET_t = PET_Static_Scan() 
            PET_t.use_compression(self._use_compression)
            PET_t.use_gpu(self._use_gpu)
            PET_t.set_scanner(self.scanner.__class__) 
            PET_t.set_binning(self.binning) 
            PET_t._load_static_measurement(t) 
            # make list of static scans
            self._dynamic.append(PET_t) 
            # also make one attribut for each static scan
            setattr(self,"frame%d"%t,self._dynamic[t]) 
            # set activity shape and size and attenuation shape and size 
            PET_t.set_activity_size(self.activity_size)
            PET_t.set_activity_shape(self.activity_shape)
            PET_t.set_attenuation_size(self.activity_size)
            PET_t.set_attenuation_shape(self.activity_shape)

        # Make a global PET_Static_Scan object 
        self.static = PET_Static_Scan()
        self.static.use_compression(self._use_compression)
        self.static.use_gpu(self._use_gpu)
        self.static.set_scanner(self.scanner.__class__) 
        self.static.set_binning(self.binning) 
        self.static._load_static_measurement() 
        self.static.set_activity_size(self.activity_size)
        self.static.set_activity_shape(self.activity_shape)
        self.static.set_attenuation_size(self.activity_size)
        self.static.set_attenuation_shape(self.activity_shape)
        
        # Free structures listmode data
        self.scanner.listmode.free_memory()
        
        # Construct ilang model 
        self._construct_ilang_model()
    
    def import_prompts(self): 
        print "Not implemented: should load prompts for multiple time frames and also set self.static.prompts to the integral."
        
    def import_randoms(self): 
        print "Not implemented: should load randoms for multiple time frames and also set self.static.randoms to the integral."

    def export_prompts(self): 
        print "Not implemented." 
        
    def export_randoms(self): 
        print "Not implemented." 
        
    def get_prompts(self): 
        prompts = []
        for frame in self: 
            prompts.append(frame.prompts)
        return prompts
               
    def set_prompts(self, prompts_list): 
        N_time_bins = len(prompts_list)
        if len(self) == N_time_bins: 
            print "PET_Dynamic_Scan.set_prompts(): Number of sinograms matches current setup; no re-initialization."
            for t in range(N_time_bins): 
                self[t].set_prompts(prompts_list[t])                
        else: 
            print "PET_Dynamic_Scan.set_prompts(): Number of sinograms does not match current setup; re-initialization."
            self._dynamic = [] 
            total_prompts = prompts_list[0]*0
            for t in range(N_time_bins): 
                PET_t = PET_Static_Scan() 
                PET_t.use_compression(self._use_compression)
                PET_t.use_gpu(self._use_gpu)
                PET_t.set_scanner(self.scanner.__class__) 
                PET_t.set_binning(self.binning) 
                PET_t.set_prompts(prompts_list[t])
                self._dynamic.append(PET_t) 
                setattr(self,"frame%d"%t,self._dynamic[t]) 
                # FIXME: remove excess frames if new N_time_bins is less than previous 
                PET_t.set_activity_size(self.activity_size)
                PET_t.set_activity_shape(self.activity_shape)
                PET_t.set_attenuation_size(self.activity_size)
                PET_t.set_attenuation_shape(self.activity_shape)
                total_prompts += prompts_list[t]

            # Make a global PET_Static_Scan object 
            self.static = PET_Static_Scan()
            self.static.use_compression(self._use_compression)
            self.static.use_gpu(self._use_gpu)
            self.static.set_scanner(self.scanner.__class__) 
            self.static.set_binning(self.binning) 
            self.static.set_prompts(total_prompts)
            self.static.set_activity_size(self.activity_size)
            self.static.set_activity_shape(self.activity_shape)
            self.static.set_attenuation_size(self.activity_size)
            self.static.set_attenuation_shape(self.activity_shape)
        
            # Construct ilang model 
            self._construct_ilang_model()

    def get_randoms(self): 
        randoms = []
        for frame in self: 
            randoms.append(frame.randoms)
        return randoms

    def set_randoms(self,randoms): 
        N_time_bins = len(randoms_list)
        if len(self) == N_time_bins: 
            print "PET_Dynamic_Scan.set_randoms(): Number of sinograms matches current setup; no re-initialization."
            for t in range(N_time_bins): 
                self[t].set_randoms(randoms_list[t])                
        else: 
            print "PET_Dynamic_Scan.set_randoms(): Number of sinograms does not match current setup; re-initialization."
            self._dynamic = [] 
            total_randoms = randoms_list[0]*0
            for t in range(N_time_bins): 
                PET_t = PET_Static_Scan() 
                PET_t.use_compression(self._use_compression)
                PET_t.use_gpu(self._use_gpu)
                PET_t.set_scanner(self.scanner.__class__) 
                PET_t.set_binning(self.binning) 
                PET_t.set_randoms(randoms_list[t])
                self._dynamic.append(PET_t) 
                setattr(self,"frame%d"%t,self._dynamic[t]) 
                # FIXME: remove excess frames if new N_time_bins is less than previous 
                PET_t.set_activity_size(self.activity_size)
                PET_t.set_activity_shape(self.activity_shape)
                PET_t.set_attenuation_size(self.activity_size)
                PET_t.set_attenuation_shape(self.activity_shape)
                total_randoms += randoms_list[t]

            # Make a global PET_Static_Scan object 
            self.static = PET_Static_Scan()
            self.static.use_compression(self._use_compression)
            self.static.use_gpu(self._use_gpu)
            self.static.set_scanner(self.scanner.__class__) 
            self.static.set_binning(self.binning) 
            self.static.set_randoms(total_randoms)
            self.static.set_activity_size(self.activity_size)
            self.static.set_activity_shape(self.activity_shape)
            self.static.set_attenuation_size(self.activity_size)
            self.static.set_attenuation_shape(self.activity_shape)
        
            # Construct ilang model 
            self._construct_ilang_model()

    def get_scatter(self): 
        return self.static.scatter 

    def set_scatter(self, scatter, duration_ms=None): 
        self.scatter = scatter  
        self.static.set_scatter(scatter, duration_ms) 
        for frame in range(len(self)): 
            self[frame].set_scatter(scatter, duration_ms)

    def get_sensitivity(self): 
        return self.static.sensitivity 

    def set_sensitivity(self, sensitivity): 
        #FIXME: verify type: PET_projection or nd_array (the latter only in full sampling mode)
        self.sensitivity = sensitivity
        self.static.set_sensitivity(sensitivity) 
        for frame in range(len(self)): 
            self[frame].set_sensitivity(sensitivity)     

    def get_attenuation(self):
        return self.static.attenuation
    
    def set_attenuation(self, attenuation): 
        #FIXME: verify type
        self.static.set_attenuation(attenuation) 
        for frame in range(len(self)): 
            self[frame].set_attenuation(attenuation)  

    def get_attenuation_projection(self):
        return self.static.attenuation_projection 
    
    def set_attenuation_projection(self, attenuation_projection): 
        #FIXME: verify type
        self.attenuation_projection = attenuation_projection
        self.static.set_attenuation_projection(attenuation_projection) 
        for frame in range(len(self)): 
            self[frame].set_attenuation_projection(attenuation_projection)  

    def use_compression(self, use_it): 
        self._use_compression = use_it
        
        
    def set_activity_shape(self, activity_shape): 
        if not len(activity_shape) == 3: 
            print "Invalid activity shape"  #FIXME: raise invalid input error
        else: 
            self.activity_shape = activity_shape
            if hasattr(self,"static"): 
                if self.static is not None: 
                    self.static.set_activity_shape(activity_shape)
            for frame in range(len(self)):
                self[frame].set_activity_shape(activity_shape)

    def set_activity_size(self, activity_size): 
        if not len(activity_size) == 3: 
            print "Invalid activity size"  #FIXME: raise invalid input error
        else: 
            self.activity_size = activity_size
            if hasattr(self,"static"): 
                if self.static is not None: 
                    self.static.set_activity_size(activity_size)
            for frame in range(len(self)):
                self[frame].set_activity_size(activity_size)

    def set_attenuation_shape(self, attenuation_shape): 
        if not len(attenuation_shape) == 3: 
            print "Invalid attenuation shape"  #FIXME: raise invalid input error
        else: 
            self.attenuation_shape = attenuation_shape
            if hasattr(self,"static"): 
                if self.static is not None: 
                    self.static.set_attenuation_shape(attenuation_shape)
            for frame in range(len(self)):
                self[frame].set_attenuation_shape(attenuation_shape)

    def set_attenuation_size(self, attenuation_size): 
        if not len(attenuation_size) == 3: 
            print "Invalid attenuation size"  #FIXME: raise invalid input error
        else: 
            self.attenuation_size = attenuation_size
            if hasattr(self,"static"): 
                if self.static is not None: 
                    self.static.set_attenuation_size(attenuation_size)
            for frame in range(len(self)):
                self[frame].set_attenuation_size(attenuation_size)

    def brain_crop(self, bin_range=(100,240)):
        if self._use_compression is True: 
            print "Cropping currently only works with uncompressed data. "
            return 
        self.static.brain_crop(bin_range)
        for frame in range(len(self)):
            self[frame].brain_crop(bin_range)
        
    def osem_reconstruction(self, iterations=10, activity=None, subset_mode="random", subset_size=64, transformations=None):
        for frame in range(len(self)):
            if activity is not None: 
                activity_init = activity[frame]
            else:
                if hasattr(self[frame],'activity'): 
                    activity_init = self[frame].activity
                else: 
                    activity_init = None
            if transformations is not None: 
                transformation = transformations[frame]
            else: 
                transformation = None
            print "Reconstructing frame %d/%d"%(frame,len(self))
            activity_recon = self[frame].osem_reconstruction(iterations=iterations, \
                             activity=activity_init, subset_mode=subset_mode, subset_size=subset_size, transformation=transformation)
            self[frame].activity = activity_recon
    
    def osem_reconstruction_4D(self, iterations=10, activity=None, subset_mode="random", subset_size=64, transformations=None):
        if activity is None: 
            activity = self._make_Image3D_activity(ones(self.activity_shape,dtype=float32,order="F")) 
            subsets_generator = SubsetGenerator(self.binning.N_azimuthal, self.binning.N_axial)

        attenuation_projection = None 
        
        if self.sensitivity is None: 
            sensitivity = self.prompts.copy()
            sensitivity.data = 0.0*sensitivity.data+1
            self.set_sensitivity(sensitivity)
        
        for i in range(iterations):
            print i
            subsets_matrix = subsets_generator.new_subset(subset_mode, subset_size)
            #subsets_matrix = ones([11,252])
            activity = self.osem_step_4D(activity, subsets_matrix, attenuation_projection, transformations)
        return activity
    
    def osem_step_4D(self, activity, subsets_matrix, attenuation_projection, transformations): 
        print "Not implemented - please implement"
        return 
    
        epsilon = 1e-08
        if attenuation_projection is not None: 
            attenuation_projection = attenuation_projection.get_subset(subsets_matrix)
        else: 
            attenuation_projection = 1.0
        
        if self.sensitivity is not None: 
            sens_x_att  = self.sensitivity.get_subset(subsets_matrix) * attenuation_projection
        else: 
            sens_x_att   = attenuation_projection
    
        if self.randoms is not None: 
            mrandoms     = self.randoms.get_subset(subsets_matrix) / (sens_x_att + epsilon) 
        
        if self.scatter is not None: 
            mscatter     = (self.scatter.get_subset(subsets_matrix) + epsilon) / (attenuation_projection + epsilon)

        norm = self.backproject_activity(sens_x_att, transformation=transformation)
    
        projection = self.project_activity(activity, subsets_matrix = subsets_matrix, transformation=transformation)
    
        if self.randoms is not None: 
            if self.scatter is not None:
                update1 = self.backproject_activity(self.prompts.get_subset(subsets_matrix)/(projection + mrandoms + mscatter + epsilon), transformation=transformation)
            else: 
                update1 = self.backproject_activity(self.prompts.get_subset(subsets_matrix)/(projection + mrandoms + epsilon), transformation=transformation)
    
        else:
            if self.scatter is not None: 
                update1 = self.backproject_activity(self.prompts.get_subset(subsets_matrix)/(projection + mscatter + epsilon), transformation=transformation)
            else:
                update1 = self.backproject_activity(self.prompts.get_subset(subsets_matrix)/(projection + epsilon), transformation=transformation)
    
        #activity = activity - gradient1
        activity = (activity /(norm+epsilon)) * update1

        return activity    
    
    

    
    def get_activity_as_array(self):
        activities = []
        for frame in range(len(self)): 
            activity = self[frame].activity.data
            activities.append(activity)
        return asarray(activities)
    
    def get_2D_slices(self, slices=[62,63,64,65,66], azimuthal_index=5): 
        N_time_bins = len(self)

        pet = PET_Multi2D_Scan()
        pet.use_compression(False)
        pet.set_scanner(self.scanner.__class__)  
        pet.set_number_of_slices( N_time_bins )
        pet.set_activity_shape( (self.activity_shape[0], self.activity_shape[1]) )
        pet.set_activity_size( (self.activity_size[0], self.activity_size[1]) )
        
        # extract slice of prompts, scatter, randoms, sensitivity
        prompts = zeros([pet.binning.N_axial, pet.binning.N_u, N_time_bins + 1], dtype=float32)
        randoms = zeros([pet.binning.N_axial, pet.binning.N_u, N_time_bins + 1], dtype=float32)
        scatter = zeros([pet.binning.N_axial, pet.binning.N_u, N_time_bins + 1], dtype=float32)
        sensitivity = ones([pet.binning.N_axial, pet.binning.N_u, N_time_bins + 1],dtype=float32)
        attenuation_projection = ones([pet.binning.N_axial, pet.binning.N_u, N_time_bins + 1],dtype=float32)
        for t in range(N_time_bins): 
            if hasattr(self[t], "prompts"): 
                if self[t].prompts is not None: 
                    prompts[:,:,t+1] = self[t].prompts.to_nd_array()[:,azimuthal_index,:,slices].sum(0).squeeze()
        for t in range(N_time_bins): 
            if hasattr(self[t], "randoms"): 
                if self[t].randoms is not None: 
                    randoms[:,:,t+1] = self[t].randoms.to_nd_array()[:,azimuthal_index,:,slices].sum(0).squeeze()
        # FIXME: implement scatter scale in PET_Static_Scan and use it here; 
        # this will reduce memory for scatter by storing scatter projection only once. 
        for t in range(N_time_bins): 
            if hasattr(self[t], "scatter"): 
                if self[t].scatter is not None: 
                    scatter[:,:,t+1] = self[t].scatter.to_nd_array()[:,azimuthal_index,:,slices].sum(0).squeeze()
        for t in range(N_time_bins): 
            if hasattr(self[t], "sensitivity"): 
                if self[t].sensitivity is not None: 
                    sensitivity[:,:,t+1] = self[t].sensitivity.to_nd_array()[:,azimuthal_index,:,slices].mean(0).squeeze() 

        # project attenuation and extract slice of the projection
        # FIXME: implement .set_attenuation_projection() in PET_Static_SCAN and use in recon if loaded 
        if hasattr(self,"attenuation_projection"):
            if self.attenuation_projection is not None: 
                att = self.attenuation_projection.to_nd_array()[:,azimuthal_index,:,slices].mean(0).squeeze()
                for t in range(N_time_bins): 
                    attenuation_projection[:,:,t+1] = att 
        
        prompts = PET_Projection( pet.binning, prompts) 
        randoms = PET_Projection( pet.binning, randoms) 
        sensitivity = PET_Projection( pet.binning, sensitivity) 
        scatter = PET_Projection( pet.binning, scatter) 
        attenuation_projection = PET_Projection( pet.binning, attenuation_projection)
        
        pet.set_prompts(prompts)
        pet.set_scatter(scatter)
        pet.set_randoms(randoms)
        pet.set_sensitivity(sensitivity)
        pet.set_attenuation_projection(attenuation_projection)
        
        prompts_duration = ones([N_time_bins+1,])
        scatter_duration = ones([N_time_bins+1,])
        for t in range(N_time_bins): 
            prompts_duration[1+t] = self[t].prompts.get_duration()
            scatter_duration[1+t] = self[t].scatter.get_duration()
        pet.set_prompts_duration(prompts_duration)
        pet.set_scatter_duration(scatter_duration)
        
        return pet
        
    def _construct_ilang_model(self):
        # define the ilang probabilistic model 
        self.ilang_model = PET_Dynamic_Poisson(self)  
        # construct the Directed Acyclical Graph
        self.ilang_graph         = ProbabilisticGraphicalModel(['lambda','alpha','counts']) 
        #self.ilang_graph.set_nodes_given(['counts','alpha'],True) 
        #self.ilang_graph.add_dependence(self.ilang_model,{'lambda':'lambda','alpha':'alpha','z':'counts'}) 

    def __repr__(self): 
        """Display information about Dynamic_PET_Scan"""
        s = "Dynamic PET acquisition:  \n" 
        s = s+" - N_time_bins:                  %d \n"%len(self.time_bins)  
        s = s+" - Time_start:                   %s \n"%millisec_to_min_sec(self.time_bins[0])
        s = s+" - Time_end:                     %s \n"%millisec_to_min_sec(self.time_bins[-1])
        s = s+" - N_counts:                     %d \n"%self.static.prompts.get_integral() 
        s = s+" - N_locations:                  %d \n"%self.static.prompts.sparsity.get_N_locations()
        s = s+" - compression_ratio:            %d \n"%self.static.prompts.get_compression_ratio()
#        s = s+" - dynamic_inflation:            %d \n"%self.dynamic_inflation
        s = s+" - listmode_loss:                %d \n"%self.static.prompts.get_listmode_loss()
        s = s+" - Mean time bin duration:       %d [sec] \n"%0 #FIXME 
        if self.scanner  is not None: 
            s = s+" * Scanner: \n" 
            s = s+"     - Name:                     %s \n"%self.scanner.model 
            s = s+"     - Manufacturer:             %s \n"%self.scanner.manufacturer 
            s = s+"     - Version:                  %s \n"%self.scanner.version 
        if self.binning  is not None:
            s = s+" * Binning: \n" 
            s = s+"     - N_axial bins:             %d \n"%self.binning.N_axial 
            s = s+"     - N_azimuthal bins:         %d \n"%self.binning.N_azimuthal 
            #s = s+"     - Angles axial step:        %f \n"%self.binning.angles_axial 
            #s = s+"     - Angles azimuthal:         %f \n"%self.binning.angles_azimuthal 
            s = s+"     - Size_u:                   %f \n"%self.binning.size_u
            s = s+"     - Size_v:                   %f \n"%self.binning.size_v        
            s = s+"     - N_u:                      %s \n"%self.binning.N_u
            s = s+"     - N_v:                      %s \n"%self.binning.N_v
        return s

    def _repr_html_(self): 
        if not has_ipy_table: 
            return "Please install ipy_table."
        table_data = [['N_time_bins',len(self.time_bins)],
        ['Time_start',millisec_to_min_sec(self.time_bins[0])],
        ['Time_end',millisec_to_min_sec(self.time_bins[-1])],
        ['Duration',millisec_to_min_sec(self.time_bins[-1]-self.time_bins[0])],
        ['N_counts',pretty_print_large_number(self.static.prompts.get_integral() )],
        ['N_locations',pretty_print_large_number(self.static.prompts.sparsity.get_N_locations())],   
        ['compression_ratio',print_percentage(self.static.prompts.get_compression_ratio() )],
#        ['dynamic_inflation',self.dynamic_inflation],
        ['listmode_loss',self.static.prompts.get_listmode_loss()], ]
        if self.scanner: 
            table_data += [['Name',self.scanner.model],['Manufacturer',self.scanner.manufacturer],['Version',self.scanner.version], ]
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()

    def __iter__(self): 
        """This method makes the object iterable. """
        return iter(self._dynamic)
    
    def __getitem__(self,i): 
        """This method makes the object addressable like a list. """
        return self._dynamic[i] 

    def __len__(self):
        return len(self._dynamic)


        
        


        
class PET_Cyclic_Scan(PET_Dynamic_Scan): 
    """PET Cyclic Scan. This is useful for respiratory gated imaging and for cardiac gated imaging (or both)."""
    
    def __init__(self): 
        PET_Dynamic_Scan.__init__(self)
        
    def import_listmode(self, hdr_filename, time_range_matrix_ms, data_filename=None, display_progress=False ): 
        """Load cyclic measurement data from a listmode file. """
        print_debug("- Loading static PET data from listmode file "+str(hdr_filename) )
        hdr = Interfile.load(hdr_filename) 
        #Extract information from the listmode header

        # 1) Guess the path of the listmode data file, if not specified or mis-specified; 
        #  1 - see if the specified listmode data file exists 
        if data_filename is not None: 
            data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep) # cross platform compatibility 
            if not os.path.exists(data_filename): 
                raise FileNotFound("listmode data",data_filename)  
        #  2 - if the listmode data file is not specified, try with the name (and full path) 
        #      contained in the listmode header
        data_filename      = hdr['name of data file']['value']
        data_filename = data_filename.replace("/",os.path.sep).replace("\\",os.path.sep) # cross platform compatibility
        if not os.path.exists(data_filename): 
        #  3 - if it doesn't exist, look in the same path as the header file for the listmode data 
        #      file with name specified in the listmode header file 
            data_filename = os.path.split(hdr_filename)[0]+os.path.sep+os.path.split(data_filename)[-1]  
            if not os.path.exists(data_filename): 
        #  4 - if it doesn't exist, look in the same path as the header file for the listmode data 
        #      file with same name as the listmode header file, replacing the extension: ".l.hdr -> .l" 
                if hdr_filename.endswith(".l.hdr"): 
                    data_filename = hdr_filename.replace(".l.hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
        #  5 - if it doesn't exist, look in the same path as the header file for the listmode data 
        #      file with same name as the listmode header file, replacing the extension: ".hdr -> .l" 
                elif hdr_filename.endswith(".hdr"): 
                    data_filename = hdr_filename.replace(".hdr",".l") 
                    if not os.path.exists(data_filename):     
                        raise FileNotFound("listmode data",data_filename)  
           
        # 2) Determine duration of the acquisition 
        n_packets              = hdr['total listmode word counts']['value'] 
        scan_duration          = hdr['image duration']['value']*1000            # milliseconds
        
        # 3) determine scanner parameters
        n_radial_bins          = hdr['number of projections']['value'] 
        n_angles               = hdr['number of views']['value'] 
        n_rings                = hdr['number of rings']['value'] 
        max_ring_diff          = hdr['maximum ring difference']['value']
        n_sinograms            = n_rings+2*n_rings*max_ring_diff-max_ring_diff**2-max_ring_diff
        n_frames      = time_range_matrix_ms.shape[0] 
        n_cycles      = time_range_matrix_ms.shape[1] 

        # 4) Display information 
        print_debug(" - Number of packets:    %d       "%n_packets )
        print_debug(" - Scan duration:        %d [sec] "%(scan_duration/1000.0) )
        print_debug(" - Listmode data file:   %s       "%data_filename )
        print_debug(" - Listmode header file: %s       "%hdr_filename)
        print_debug(" - n_frames :            %d       "%n_frames)
        print_debug(" - n_cycles :            %d       "%n_cycles)
        print_debug(" - n_radial_bins:        %d       "%n_radial_bins)
        print_debug(" - n_angles:             %d       "%n_angles)
        print_debug(" - n_angles:             %d       "%n_sinograms)

        if display_progress: 
            progress_bar = ProgressBar()
            progress_callback = progress_bar.set_percentage
        else: 
            def progress_callback(value): 
                pass 
        
        # Load the listmode data 
        M = self.scanner.michelogram
        R = self.scanner.listmode.load_listmode_cyclic(data_filename, time_range_matrix_ms, self.binning, n_radial_bins, 
                                              n_angles, n_sinograms, M.span, M.segments_sizes, M.michelogram_sinogram,
                                              M.michelogram_plane, n_packets, progress_callback) 
        progress_callback(100) 

        self._dynamic = [] 
        for t in range(n_frames): 
            PET_t = PET_Static_Scan() 
            PET_t.use_compression(self._use_compression)
            PET_t.use_gpu(self._use_gpu)
            PET_t.set_scanner(self.scanner.__class__) 
            PET_t.set_binning(self.binning) 
            PET_t._load_static_measurement(t) 
            # make list of static scans
            self._dynamic.append(PET_t) 
            # also make one attribut for each static scan
            setattr(self,"frame%d"%t,self._dynamic[t]) 
            # set activity shape and size and attenuation shape and size 
            PET_t.set_activity_size(self.activity_size)
            PET_t.set_activity_shape(self.activity_shape)
            PET_t.set_attenuation_size(self.activity_size)
            PET_t.set_attenuation_shape(self.activity_shape)

        # Make a global PET_Static_Scan object 
        self.static = PET_Static_Scan()
        self.static.use_compression(self._use_compression)
        self.static.use_gpu(self._use_gpu)
        self.static.set_scanner(self.scanner.__class__) 
        self.static.set_binning(self.binning) 
        self.static._load_static_measurement() 
        self.static.set_activity_size(self.activity_size)
        self.static.set_activity_shape(self.activity_shape)
        self.static.set_attenuation_size(self.activity_size)
        self.static.set_attenuation_shape(self.activity_shape)

        # Free structures listmode data
        self.scanner.listmode.free_memory()
        
        # Construct ilang model 
        self._construct_ilang_model()

        

