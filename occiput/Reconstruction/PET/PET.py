
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Jan. 2015, Boston
# June 2015, Helsinki 
# Nov. 2015, Boston 


# If you are looking for PET reconstruction, this is where to start. 
# The objects defined here provide a abstractions for Static and Dynamic PET reconstruction, 
# abstracting the scanner geometries and vendor models and providing an interface to the 
# software tools for projection, backprojection and reconstruction. 


__all__ = ['PET_Static_Scan','PET_Dynamic_Scan','PET_Cyclic_Scan','Binning','PET_Projection_Sparsity','PET_Projection'] 


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
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 
import os
import h5py 
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

        



# FIXME: eliminate the class ROI; use transformation matrices in Image3D instead, for activity and attenuation volumes. 
# Use Core.Transform_6DOF or Core.Transform_Affine if required, to parameterize the projector and back_projector. 

class ROI(): 
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
                raise UnknownParameter('Parameter %s specified for the construction of ROI is not compatible. '%str(parameters)) 
        else: 
            raise UnknownParameter('Parameter %s specified for the construction of ROI is not compatible. '%str(parameters)) 
    
    def load_from_dictionary(self,dictionary):
        self.x       = dictionary['x']                           # Translation along x
        self.y       = dictionary['y']                           # Translation along y
        self.z       = dictionary['z']                           # Translation along z
        self.theta_x = dictionary['theta_x']                     # Rotation around x
        self.theta_y = dictionary['theta_y']                     # Rotation around y
        self.theta_z = dictionary['theta_z']                     # Rotation around z

    def __repr__(self):
        s = "PET volume location (ROI): \n"
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
        self.sensitivity   = None                        # sensitivity is a permanent parameter. 
        self.prompts       = None                        # measurement: prompts. Initialized as empty data structure.
        self.randoms       = None
        self.scatter       = None
        self._normalization = None                        # normalization volume - for all projections - memoize
        self._need_normalization_update = True            # If True, the normalization volume needs to be recomputed 
        self.use_compression(True)

        self.set_transform_scanner_to_world(Transform_Identity(map_from='scanner',map_to='world'))
        
        self._construct_ilang_model() 

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
                    print "Not able to compress once uncompressed. Please implement PET_Projection.uncompress_self() to "
                    print "enable this functionality. "
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
        self._use_compression = use_it 

    def get_scatter(self): 
        return self.scatter 

    def set_scatter(self, scatter): 
        self.scatter = scatter  

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
            if self.prompts is not None:  # FIXME: sensitivity loaded from interfile with some manufacturers has non-zero value 
                                          # where there are no detectors - set to zero where data is zero 
                                          #(good approx only for long acquisitions). See if there is a better 
                                          # way to handle this. 
                sensitivity.data[self.prompts.data==0]=0
            else: 
                print "Please load prompts before loading the sensitivity. " # FIXME: see comment two lines up 
        elif filetype is "mat": 
            print "Sensitivity from Matlab not yet implemented. All is ready, please spend 15 minutes and implement. "
            return 
        else: 
            print "File type unknown. "
            return 
        if self._use_compression is False: 
            sensitivity = sensitivity.uncompress_self() 
        self.set_sensitivity(sensitivity)

    def import_scatter(self, filename, datafile=''): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
            projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_scatter: file type unknown. "
            return 
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_scatter(projection) 

    def import_randoms(self, filename, datafile=''): 
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
             projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_randoms: file type unknown. "
            return
        if self._use_compression is False: 
            projection = projection.uncompress_self() 
        self.set_randoms(projection)

    def import_prompts(self, filename, datafile=''):
        filetype = guess_file_type_by_name(filename) 
        if filetype is "interfile_projection_header": 
            projection = import_interfile_projection(filename, self.binning, self.scanner.michelogram, datafile)
        elif filetype is "h5": 
             projection = import_PET_Projection(filename)
        else: 
            print "PET.import_prompts: file type unknown. "
            return
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
        self.set_attenuation(volume)

    def import_attenuation_projection(self, filename, datafile=''): 
        print "Not implemented. All is ready, please spend 10 minutes to implement it if you would like this feature. "
        # FIXME: support interfile, h5, matlab, nifti 


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
        time_bins = int32(linspace(time_range_ms[0],time_range_1,2))

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

        # Construct ilang model 
        self._construct_ilang_model()
    
    def quick_inspect(self): 
        if self.randoms is not None and not isscalar(self.randoms): 
            randoms = self.randoms.to_nd_array()[0,5,:,60]
        else: 
            randoms = 0.0 
        if self.prompts is not None and not isscalar(self.prompts): 
            prompts = self.prompts.to_nd_array()[0,5,:,60]
        else: 
            prompts = 0.0 
        if self.sensitivity is not None and not isscalar(self.sensitivity): 
            sensitivity = self.sensitivity.to_nd_array()[0,5,:,60]
        else: 
            sensitivity = 1.0 
        if self.scatter is not None and not isscalar(self.scatter): 
            scatter = self.scatter.to_nd_array()[0,5,:,60]
        else: 
            scatter = 0.0 
        if has_pylab: 
            pylab.plot(prompts - randoms)
            pylab.hold(1)
            pylab.plot(sensitivity*scatter,'g') 
        else: 
            print "quick_inspect uses Pylab to display imaging data. Please install Pylab. " 

    def _get_sparsity(self): 
        if self.prompts == None: 
            # in this case returns sparsity pattern for uncompressed projection 
            sparsity = PET_Projection_Sparsity(self.binning.N_axial, self.binning.N_azimuthal, self.binning.N_u, self.binning.N_v)
        else: 
            sparsity = self.prompts.sparsity
        return sparsity 

    def set_scale_activity(self,scale): 
        self.scale_activity = scale 


    def project_attenuation(self, attenuation, unit='inv_cm', roi=None, sparsity=None, subsets_matrix=None, exponentiate=True): 
        if isinstance(attenuation,ndarray): 
            attenuation_data = f_continuous(float32(attenuation))
        else: 
            attenuation_data = f_continuous(float32(attenuation.data))

        if not list(attenuation_data.shape) == list(self.attenuation_shape): 
            raise UnexpectedParameter("Attenuation must have the same shape as self.attenuation_shape")

        # By default, the center of the imaging volume is at the center of the scanner 
        if roi is None:  
            tx = 0.5*(self.attenuation_size[0] - self.attenuation_size[0] / self.attenuation_shape[0])
            ty = 0.5*(self.attenuation_size[1] - self.attenuation_size[1] / self.attenuation_shape[1])
            tz = 0.5*(self.attenuation_size[2] - self.attenuation_size[2] / self.attenuation_shape[2])
            roi = ROI((tx,ty,tz,0,0,0))

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
        if sparsity == None: 
            sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            subsets_matrix=self._subsets_generator.all_active()

        # Call the raytracer 
        projection_data = PET_project_compressed(attenuation_data,None,sparsity.offsets,sparsity.locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, 
            float32(self.binning.angles_axial), float32(self.binning.angles_azimuthal), 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            0.0,0.0,0.0, 
            roi.x, roi.y, roi.z, roi.theta_x, roi.theta_y, roi.theta_z, 
            0.0,0.0,0.0,0.0,0.0,0.0, 
            self.attenuation_projection_parameters.gpu_acceleration, self.attenuation_projection_parameters.N_samples,
            self.attenuation_projection_parameters.sample_step, self.attenuation_projection_parameters.background_attenuation, 
            0.0, self.attenuation_projection_parameters.truncate_negative_values, 
            self.attenuation_projection_parameters.direction, self.attenuation_projection_parameters.block_size) 
        
        # Fix scale and exponentiate 
        if exponentiate: 
            projection_data = exp(-float64(projection_data) * step_size) 
        else: 
            projection_data = projection_data * step_size

        # Create object PET_Projection: it contains the raw projection data and the description of the projection geometry 
        # and sparsity pattern. 
        time_bins = int32([0,0])  # Projection of the attenuation does not have timing information 
        projection = PET_Projection( self.binning, projection_data, sparsity.offsets, sparsity.locations, time_bins)  
        return projection 


    def backproject_attenuation(self, projection, unit="inv_cm", roi=None, sparsity=None, subsets_matrix=None): 
        if isinstance(projection,ndarray): 
            projection_data = float32(projection)
        else: 
            projection_data = float32(projection.data)

        # By default, the center of the imaging volume is at the center of the scanner 
        if roi is None:  
            tx = 0.5*(self.attenuation_size[0] - self.attenuation_size[0] / self.attenuation_shape[0])
            ty = 0.5*(self.attenuation_size[1] - self.attenuation_size[1] / self.attenuation_shape[1])
            tz = 0.5*(self.attenuation_size[2] - self.attenuation_size[2] / self.attenuation_shape[2])
            roi = ROI((tx,ty,tz,0,0,0)) 

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

        # Optionally specify a sparsity pattern. Otherwise use the sparsity pattern in 'projection'. 
        # If projection is raw data, by default use the sparsity pattern associated to the prompts. 
        # Note: if prompts have not been loaded, self._get_sparsity() assumes no compression. 
        if sparsity is None: 
            if hasattr(projection,"sparsity"): 
                sparsity = projection.sparsity 
            else: 
                sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            subsets_matrix=self._subsets_generator.all_active() 

        # Call ray-tracer 
        backprojection_data = PET_backproject_compressed(projection_data,None,sparsity.offsets, 
            sparsity.locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, float32(self.binning.angles_axial),
            float32(self.binning.angles_azimuthal), 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.attenuation_shape[0], self.attenuation_shape[1], self.attenuation_shape[2], 
            self.attenuation_size[0], self.attenuation_size[1], self.attenuation_size[2], 
            0.0, 0.0, 0.0, 
            roi.x, roi.y, roi.z, roi.theta_x, roi.theta_y, roi.theta_z, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            self.attenuation_backprojection_parameters.gpu_acceleration, self.attenuation_backprojection_parameters.N_samples,
            self.attenuation_backprojection_parameters.sample_step, 
            self.attenuation_backprojection_parameters.background_attenuation, 0.0, 
            self.attenuation_backprojection_parameters.direction, self.attenuation_backprojection_parameters.block_size) 

        backprojection_data = backprojection_data * step_size 
        
        # Set the correct scale - unit measure and return Image3D - FIXME: set scale for requested unit measure 
        return self._make_Image3D_attenuation(backprojection_data)  


    def project_activity(self, activity, unit="Bq/mm3", roi=None, sparsity=None, subsets_matrix=None): 
        if isinstance(activity,ndarray): 
            activity_data = f_continuous(float32(activity))
        else: 
            activity_data = f_continuous(float32(activity.data))

        if not list(activity_data.shape) == list(self.activity_shape): 
            raise UnexpectedParameter("Activity must have the same shape as self.activity_shape")

        # By default, the center of the imaging volume is at the center of the scanner 
        if roi  is None:  
            tx = 0.5*(self.activity_size[0] - self.activity_size[0] / self.activity_shape[0])
            ty = 0.5*(self.activity_size[1] - self.activity_size[1] / self.activity_shape[1])
            tz = 0.5*(self.activity_size[2] - self.activity_size[2] / self.activity_shape[2])
            roi = ROI((tx,ty,tz,0,0,0))

        # Optionally project with a sparsity pattern not equal to sparsity associated to the loaded prompts data 
        # Note: if prompts have not been loaded, self._get_sparsity() assumes no compression. 
        if sparsity is None: 
            sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            subsets_matrix=self._subsets_generator.all_active()
        
        scale = 1.0  #FIXME: change this according to the input unit measure - check how this is done in project_attenuation
        step_size_mm = self.activity_projection_parameters.sample_step 
        step_size = step_size_mm / scale

        # Call the raytracer 
        projection_data = PET_project_compressed(activity_data,None,sparsity.offsets,sparsity.locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, 
            float32(self.binning.angles_axial), float32(self.binning.angles_azimuthal), 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            0.0,0.0,0.0, 
            roi.x, roi.y, roi.z, roi.theta_x, roi.theta_y, roi.theta_z, 
            0.0,0.0,0.0,0.0,0.0,0.0, 
            self.activity_projection_parameters.gpu_acceleration, self.activity_projection_parameters.N_samples,
            self.activity_projection_parameters.sample_step, self.activity_projection_parameters.background_activity, 
            0.0, self.activity_projection_parameters.truncate_negative_values, 
            self.activity_projection_parameters.direction, self.activity_projection_parameters.block_size) 

        # Create object PET_Projection: it contains the raw projection data and the description of the projection geometry 
        # and sparsity pattern. 
        time_bins = int32([0,1000.0])     # 1 second - projection returns a rate - by design 
        projection_data = projection_data * step_size
        projection = PET_Projection( self.binning, projection_data, sparsity.offsets, sparsity.locations, time_bins)  

        # Optionally scale by sensitivity, attenuation, time, global sensitivity 
        #if attenuation is not None: 
        #    projection.data = projection.data * attenuation.data #attenuation.compress_as(projection).data
        #if self.sensitivity is not None: 
        #    projection.data = projection.data * self.sensitivity.data.reshape(projection_data.shape)
        #self.sensitivity.compress_as(projection).data
        return projection 


    def backproject_activity(self, projection, roi=None, sparsity=None, subsets_matrix=None): 
        if isinstance(projection,ndarray): 
            projection_data = float32(projection)
        else: 
            projection_data = float32(projection.data)

        # Optionally multiply by attenuation and sensitivity  
        #if attenuation is not None:  
        #    projection_data = projection_data * attenuation.data # attenuation.compress_as(self.prompts).data
        #if self.sensitivity is not None: 
        #    projection_data = projection_data * self.sensitivity.data.reshape(projection_data.shape) 
        #self.sensitivity.compress_as(self.prompts).data

        # By default, the center of the imaging volume is at the center of the scanner 
        if roi is None:  
            tx = 0.5*(self.activity_size[0] - self.activity_size[0] / self.activity_shape[0])
            ty = 0.5*(self.activity_size[1] - self.activity_size[1] / self.activity_shape[1])
            tz = 0.5*(self.activity_size[2] - self.activity_size[2] / self.activity_shape[2])
            roi = ROI((tx,ty,tz,0,0,0)) 

        # Optionally specify a sparsity pattern. Otherwise use the sparsity pattern in 'projection'. 
        # If projection is raw data, by default use the sparsity pattern associated to the prompts. 
        # Note: if prompts have not been loaded, self._get_sparsity() assumes no compression. 
        if sparsity is None: 
            if hasattr(projection,"sparsity"): 
                sparsity = projection.sparsity 
            else: 
                sparsity = self._get_sparsity() 

        # Optionally project only to a subset of projection planes  
        if subsets_matrix is None: 
            subsets_matrix=self._subsets_generator.all_active()

        scale = 1.0  #FIXME: change this according to the input unit measure - check how this is done in project_attenuation
        step_size_mm = self.activity_projection_parameters.sample_step 
        step_size = step_size_mm / scale

        # Call ray-tracer 
        backprojection_data = PET_backproject_compressed(projection_data,None,sparsity.offsets, 
            sparsity.locations, subsets_matrix, 
            self.binning.N_axial, self.binning.N_azimuthal, float32(self.binning.angles_axial),
            float32(self.binning.angles_azimuthal), 
            self.binning.N_u, self.binning.N_v, self.binning.size_u, self.binning.size_v, 
            self.activity_shape[0], self.activity_shape[1], self.activity_shape[2], 
            self.activity_size[0], self.activity_size[1], self.activity_size[2], 
            0.0, 0.0, 0.0, 
            roi.x, roi.y, roi.z, roi.theta_x, roi.theta_y, roi.theta_z, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            self.activity_backprojection_parameters.gpu_acceleration, self.activity_backprojection_parameters.N_samples,
            self.activity_backprojection_parameters.sample_step, 
            self.activity_backprojection_parameters.background_activity, 0.0, 
            self.activity_backprojection_parameters.direction, self.activity_backprojection_parameters.block_size) 

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
        self.set_prompts( PET_Projection( self.binning, Rp['counts'], Rp['offsets'], Rp['locations'], time_bins) )
        self.set_randoms( PET_Projection( self.binning, Rd['counts'], Rd['offsets'], Rd['locations'], time_bins) )
        self._construct_ilang_model() 
        

    def get_activity(self):
        return self.activity 

    def set_activity(self,activity): 
        self.ilang_graph.set_node_value('lambda',activity) 
        print "PET_Static_Scan.set_activity(): This is for integration with iLang - please implement"
        # FIXME: how about the roi ? 

    def get_attenuation(self):
        return self.attenuation 

    def set_attenuation(self,attenuation): 
        self.ilang_graph.set_node_value('alpha',attenuation) 
        # FIXME: how about the roi ? 
        # FIXME: setting activity and attenuation as members here is not in the spirit of iLang - memoization 
        self.attenuation = attenuation 

    def get_normalization(self, attenuation_times_sensitivity=None, roi=None, sparsity=None, duration_ms=None, subsets_matrix=None, epsilon=None): 
        # FIXME: memoization 
        if attenuation_times_sensitivity is None: 
            attenuation_times_sensitivity = self.sensitivity  #FIXME: include attenuation here - memoization mumap proj
        if isscalar(attenuation_times_sensitivity) or attenuation_times_sensitivity is None: 
            attenuation_times_sensitivity = ones(self.prompts.data.shape)
        if duration_ms is None: 
            duration_ms = self.prompts.get_duration() 
        duration_sec = duration_ms / 1000.0 
        alpha = self.scale_activity 
        normalization = self.backproject_activity(attenuation_times_sensitivity*duration_sec*alpha,roi,sparsity,subsets_matrix) 
        return normalization 

    def get_gradient_activity(self, activity, attenuation=None, unit_activity="Bq/mm3", roi_activity=None, 
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
        if sensitivity is None: 
            sensitivity = 1.0 
        att_sens = sensitivity * attenuation

        # Compute the firt term of the gradient: backprojection of the sensitivity of the scanner
        # If it is requested that the gradient is computed using all the projection measurements, use the 
        #memoized normalization. self.get_normalization() takes care of memoization.  
        gradient_term1 = self.get_normalization(att_sens, roi_activity, sparsity, duration_ms, subsets_matrix, epsilon=epsilon)

        # Compute the second term of the gradient: backprojection of the ratio between the measurement and the projection of 
        # current activity estimate... Ordinary Poisson to include scatter and randoms.
        projection = self.project_activity(activity, unit=unit_activity, roi=roi_activity, sparsity=sparsity,
                                           subsets_matrix=subsets_matrix) 
        gradient_term2 = self.backproject_activity( prompts / (projection+randoms/(att_sens*alpha*duration_sec+epsilon)+scatter/(attenuation*alpha*duration_sec+epsilon)+epsilon),
                              roi=roi_activity, sparsity=sparsity, subsets_matrix=subsets_matrix)

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

        attenuation_projection = self.project_attenuation(attenuation, unit='inv_cm', roi=None, sparsity=sparsity, subsets_matrix=subsets_matrix, exponentiate=True) 
        
        #FIXME: roi = None 
        pr_activity = self.project_activity(activity, roi=None, sparsity=sparsity,
                                           subsets_matrix=subsets_matrix)*sensitivity*attenuation_projection*duration_sec*alpha
        gradient = self.backproject_attenuation(pr_activity - prompts/ (randoms/(pr_activity+epsilon) +  scatter/(pr_activity/(sensitivity+epsilon)+epsilon) + 1 ), unit="inv_cm", roi=None, sparsity=sparisty, 
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
        # FIXME: use transformation - also notice that roi_activity is always set to None here
    
        for iteration in range(iterations): 
            [gradient1, gradient2, subsets_matrix] = self.get_gradient_activity(activity,self.attenuation_projection,
                                                                         roi_activity=None, 
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
        proj = PET_project_compressed(volume, None, offsets,locations, subsets_matrix, 
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

    def __del__(self):
        self.scanner.listmode.free_memory()
        del(self.scanner.listmode) 

        
        
class PET_Cyclic_Scan(PET_Static_Scan): 
    """PET Dynamic Cyclic Scan. """
    def load_listmode_file(self, hdr_filename, time_range_matrix_ms, data_filename=None, display_progress=False ): 
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
        self.set_span(11) 
        M = self.michelogram
        R = self.scanner.listmode.load_listmode_cyclic(data_filename, time_range_matrix_ms, self.binning, n_radial_bins, 
                                              n_angles, n_sinograms, M.span, M.segments_sizes, M.michelogram_sinogram,
                                              M.michelogram_plane, n_packets, progress_callback) 
        progress_callback(100) 

        self._dynamic = [] 
        for t in range(n_frames): 
            PET_t = PET_Static_Scan() 
            PET_t.set_scanner(self.scanner) 
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
            PET_t.set_span(self.michelogram.span)

        # Make a global PET_Static_Scan object 
#        self.static = PET_Static_Scan()
#        self.static.set_scanner(self.scanner) 
#        self.static.set_binning(self.binning) 
#        self.static._load_static_measurement() 

        # Load static measurement data 
        self._load_static_measurement() 

        # Construct ilang model 
        self._construct_ilang_model()
        
        #return self 



        
        
        
        
class PET_Dynamic_Scan(): 
    """PET Dynamic Scan. """
    def __init__(self): 
        self.scanner     = None                                 # PET scanner interface - by default set Generic PET scanner 
        self.set_scanner("Generic")                   

        self._dynamic      = []                                 # Sequence of static scans, one per time bin 
        self.time_bins     = []                                 # Time binning 
        self.static        = None

        self.projection_parameters     = ProjectionParameters()       
        self.backprojection_parameters = BackprojectionParameters()   
        self.set_gpu_acceleration(True)                         # change to self.disable_gpu_acceleration() to disable by default      

        self.activity      = None                               # memoization of activity  
        self.attenuation   = None  
        self.sensitivity   = None 
        #self.normalization = None                              # Normalization volume 
        #self._need_normalization_update = True                  # If True, the normalization volume needs to be recomputed 

        self.set_activity_shape(DEFAULT_ACTIVITY_SHAPE)  
        self.set_activity_size(DEFAULT_ACTIVITY_SIZE)   
        self.set_attenuation_shape(DEFAULT_ACTIVITY_SHAPE)  
        self.set_attenuation_size(DEFAULT_ACTIVITY_SIZE)   

        self._construct_ilang_model() 
        #self._display_node = DisplayNode() 

    def set_activity_shape(self, activity_shape): 
        if not len(activity_shape) == 3: 
            print "Invalid activity shape"  #FIXME: raise invalid input error
        else: 
            self.activity_shape = activity_shape 
            for frame in self._dynamic: 
                frame.set_activity_shape(activity_shape) 

    def set_activity_size(self, activity_size): 
        if not len(activity_size) == 3: 
            print "Invalid activity size"  #FIXME: raise invalid input error
        else: 
            self.activity_size = activity_size 
            for frame in self._dynamic: 
                frame.set_activity_size(activity_size)
                
    def set_attenuation_shape(self, attenuation_shape): 
        if not len(attenuation_shape) == 3: 
            print "Invalid attenuation shape"  #FIXME: raise invalid input error
        else: 
            self.attenuation_shape = attenuation_shape 
            for frame in self._dynamic: 
                frame.set_attenuation_shape(activity_shape)
                
    def set_attenuation_size(self, attenuation_size): 
        if not len(attenuation_size) == 3: 
            print "Invalid attenuation size"  #FIXME: raise invalid input error
        else: 
            self.attenuation_size = attenuation_size 
            for frame in self._dynamic: 
                frame.set_attenuation_size(activity_size)

    def _construct_ilang_model(self):
        # define the ilang probabilistic model 
        self.ilang_model = PET_Dynamic_Poisson(self)  
        # construct the Directed Acyclical Graph
        self.ilang_graph         = ProbabilisticGraphicalModel(['lambda','alpha','counts']) 
#        self.ilang_graph.set_nodes_given(['counts','alpha'],True) 
#        self.ilang_graph.add_dependence(self.ilang_model,{'lambda':'lambda','alpha':'alpha','z':'counts'}) 

    def set_binning(self, binning): 
        if isinstance(binning,Binning): 
            self.binning = binning
        else:    
            self.binning = Binning(binning)
        if self._dynamic is not None: 
            for static in self._dynamic: 
                static.set_binning(binning)
        return self.binning 

    def set_scanner(self,scanner): 
        self.scanner = scanner 

    def save_prompts_to_file(self,filename='./sinogram.h5'):
        self.get_prompts().save_to_file(filename) 

    def plot_motion_events(self): 
        self.__motion_events.plot_mean_displacement()
 
    def plot_motion_parameters(self): 
        self.__motion_events.plot_motion()
        
    def set_span(self, span, n_rings=64, max_ring_difference=60):     #FIXME: temporary, this will not be here but in the data loading machinery
        self.michelogram = Michelogram(n_rings=n_rings, span=span, max_ring_difference=max_ring_difference)
        parameters = {      "n_axial":                252,
                            "n_azimuthal":            self.michelogram.n_segments,
                            "angles_axial":           float32( linspace(0,pi-pi/252.0,252) ),
                            "angles_azimuthal":       float32( linspace(deg_to_rad(-26.816),deg_to_rad(26.816),span) ),
                            "size_u":                 437.98,
                            "size_v":                 255.00,
                            "n_u":                    344,
                            "n_v":                    self.michelogram.segments_sizes.max(),   }
        binning = Binning(parameters)
        self.set_binning(binning)

    def set_sensitivity(self, sensitivity): 
        #FIXME: verify type 
        self.sensitivity = sensitivity 

    def save_sensitivity_to_file(self, filename): 
        if self.sensitivity is None: 
            print "Sensitivity has not been loaded"
        else: 
            self.sensitivity.save_to_file(filename)

    def import_sensitivity_sinogram(self, headerfile, datafile):
        print "- implement me -" 

    def load_listmode_file(self, hdr_filename, time_range_ms=[0,None], data_filename=None, motion_files_path=None, display_progress=False ): 
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
        self.set_span(11) 
        M = self.michelogram
        R = self.scanner.load_listmode(data_filename, n_packets, time_range_ms, self.binning, n_radial_bins, n_angles, n_sinograms, M.span, M.segments_sizes, M.michelogram_sinogram, M.michelogram_plane, progress_callback) 
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
            PET_t.set_scanner(self.scanner) 
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
            PET_t.set_span(self.michelogram.span)

        # Make a global PET_Static_Scan object 
        self.static = PET_Static_Scan()
        self.static.set_scanner(self.scanner) 
        self.static.set_binning(self.binning) 
        self.static._load_static_measurement() 

        # Construct ilang model 
        self._construct_ilang_model()
        
        #return self 

    def set_gpu_acceleration(self, on=True): 
        if self._dynamic is not None: 
             for frame in self._dynamic: 
                 frame.set_gpu_acceleration(on) 
        if self.static is not None: 
            self.static.set_gpu_acceleration(on) 
    
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

    def __del__(self):
        """Delete interface when the object is deleted: interface needs explicit 
        deletion in order to manage C library memory deallocation  """
        self.scanner.free_memory()

        
        


