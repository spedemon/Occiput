
# NiftyRec - Ray-tracing tools
# Stefano Pedemonte
# Center for Medical Image Computing (CMIC), University College Lonson (UCL)
# 2009-2012, London
# Aalto University, School of Science
# Summer 2013, Helsinki
# Martinos Center for Biomedical Imaging, Harvard University/MGH
# Jan. 2014, Boston

from simplewrap import *
import numpy
import os, platform

__all__ = ['test_library_niftyrec_c','gpu_set','gpu_reset','gpu_list','gpu_exists',
'PET_project','PET_backproject','PET_project_compressed','PET_backproject_compressed',
'PET_compress_projection','PET_uncompress_projection','PET_initialize_compression_structure',
'SPECT_project_parallelholes','SPECT_backproject_parallelholes',
'CT_project_conebeam','CT_backproject_conebeam','CT_project_parallelbeam','CT_backproject_parallelbeam',
'ET_spherical_phantom','ET_cylindrical_phantom','ET_spheres_ring_phantom', 
'TR_grid_from_box_and_affine', 'TR_resample_grid', 'TR_resample_box', 'TR_gradient_grid', 'TR_gradient_box', 'TR_transform_grid',
 'INTERPOLATION_LINEAR','INTERPOLATION_POINT'] 

library_name = "_et_array_interface"
niftyrec_lib_paths = [localpath(), filepath(__file__), './', '/usr/local/niftyrec/lib/', 'C:/Prorgam Files/NiftyRec/lib/'] 

INTERPOLATION_LINEAR = 0
INTERPOLATION_POINT  = 1


####################################### Error handling: ########################################

class ErrorInCFunction(Exception): 
    def __init__(self,msg,status,function_name): 
        self.msg = str(msg) 
        self.status = status
        self.function_name = function_name
        if self.status == status_io_error(): 
            self.status_msg = "IO Error"
        elif self.status == status_initialisation_error(): 
            self.status_msg = "Error with the initialisation of the C library"
        elif self.status == status_parameter_error(): 
            self.status_msg = "One or more of the specified parameters are not right"
        elif self.status == status_unhandled_error(): 
            self.status_msg = "Unhandled error, likely a bug. "
        else: 
            self.status_msg = "Unspecified Error"
    def __str__(self): 
        return "'%s' returned by the C Function '%s' (error code %d). %s"%(self.status_msg,self.function_name,self.status,self.msg)

class ErrorGPU(Exception): 
    def __init__(self,msg): 
        self.msg = msg 
    def __str__(self): 
        return "%s"%str(self.msg)


def status_success(): 
    """Returns the value returned by the function calls to the library in case of success. """
    r = call_c_function( niftyrec_c.status_success, [{'name':'return_value',  'type':'uint', 'value':None}] ) 
    return r.return_value

def status_io_error(): 
    """Returns the integer value returned by the function calls to the library in case of IO error. """
    r = call_c_function( niftyrec_c.status_io_error, [{'name':'return_value',  'type':'uint', 'value':None}] ) 
    return r.return_value

def status_initialisation_error(): 
    """Returns the value returned by the function calls to the library in case of initialisation error. """
    r = call_c_function( niftyrec_c.status_initialisation_error, [{'name':'return_value',  'type':'uint', 'value':None}] ) 
    return r.return_value

def status_parameter_error(): 
    """Returns the value returned by the function calls to the library in case of parameter error. """
    r = call_c_function( niftyrec_c.status_parameter_error, [{'name':'return_value',  'type':'uint', 'value':None}] ) 
    return r.return_value

def status_unhandled_error(): 
    """Returns the value returned by the function calls to the library in case of unhandled error. """
    r = call_c_function( niftyrec_c.status_unhandled_error, [{'name':'return_value',  'type':'uint', 'value':None}] ) 
    return r.return_value


class LibraryNotFound(Exception): 
    def __init__(self,msg): 
        self.msg = msg 
    def __str__(self): 
        return "Library cannot be found: %s"%str(self.msg) 


####################################### Load library: ########################################
def test_library_niftyrec_c(): 
    """Test whether the C library niftyrec_c responds. """
    number = 101 # just a number
    descriptor = [  {'name':'input',  'type':'int', 'value':number},
                    {'name':'output', 'type':'int', 'value':None },  ]
    r = call_c_function( niftyrec_c.echo, descriptor ) 
    return r.output == number
    

# search for the library in the list of locations 'niftyrec_lib_paths' 
# and in the locations in the environment variables "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH" and "PATH"

if platform.system()=='Linux':
    sep = ":"
elif platform.system()=='Darwin':
    sep = ":"
elif platform.system()=='Windows':
    sep = ";"
if os.environ.has_key('LD_LIBRARY_PATH'): 
    niftyrec_lib_paths = niftyrec_lib_paths + os.environ['LD_LIBRARY_PATH'].split(sep)
if os.environ.has_key('DYLD_LIBRARY_PATH'): 
    niftyrec_lib_paths = niftyrec_lib_paths + os.environ['DYLD_LIBRARY_PATH'].split(sep)
if os.environ.has_key('PATH'): 
    niftyrec_lib_paths = niftyrec_lib_paths + os.environ['PATH'].split(sep)

(found,fullpath,path) = find_c_library(library_name,niftyrec_lib_paths) 
if found == NOT_FOUND: 
    print "The library %s cannot be FOUND, please make sure that the path to the NiftyRec libraries has been exported. "%fullpath
    print "1) Before launching Python, type the following in the terminal (the same terminal): "
    path = "'path to the niftyrec libraries'"
    if platform.system()=='Linux':
        print "export LD_LIBRARY_PATH=%s"%path
    elif platform.system()=='Darwin':
        print "export DYLD_LIBRARY_PATH=%s"%path
    elif platform.system()=='Windows':
        print "Add %s to the system PATH using Control Panel -> Advanced Settings -> System -> .."%path
    raise LibraryNotFound("NiftyRec"+" ("+str(library_name)+") ")
elif found == FOUND_NOT_LOADABLE: 
    print "The library %s cannot be LOADED, please make sure that the path to the NiftyRec libraries has been exported. "%fullpath
    print "1) Before launching Python, type the following in the terminal (the same terminal): "
    if platform.system()=='Linux':
        print "export LD_LIBRARY_PATH=%s"%path
    elif platform.system()=='Darwin':
        print "export DYLD_LIBRARY_PATH=%s"%path
    elif platform.system()=='Windows':
        print "Add %s to the system PATH using Control Panel -> Advanced Settings -> System -> .."%path
    raise LibraryNotFound("NiftyRec"+" ("+str(library_name)+") ")
else: 
    niftyrec_c = load_c_library(fullpath)

#################################### Create interface to the C functions: ####################################

def gpu_list(): 
    """List GPUs and get information. """
    MAX_GPUs  = 1000
    INFO_SIZE = 5
    descriptor = [{'name':'N',     'type':'array',  'value':None,   'dtype':int32,  'size':(1,) }, 
                  {'name':'info',  'type':'array',  'value':None,   'dtype':int32,  'size':(MAX_GPUs,INFO_SIZE) }, ]
    r = call_c_function( niftyrec_c.et_array_list_gpus, descriptor) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'gpu_list' was unsuccessful.",r.status,'niftyrec_c.et_array_list_gpus') 
    N = r.dictionary['N'][0]
    info = r.dictionary['info'] 
    gpus = []
    for i in range(N): 
        gpus.append({'id':info[i,0], 'gflops':info[i,1], 'multiprocessors':info[i,2], 'clock':info[i,3], 'globalmemory':info[i,4] })
    return gpus 

def gpu_exists(gpu_id): 
    for gpu in gpu_list(): 
        if gpu['id']==gpu_id: 
            return True
    return False 

def gpu_set(gpu_id=0): 
    """Set GPU (when multiple GPUs are installed in the system). """
    if not gpu_exists(gpu_id): 
        gpus = gpu_list()
        raise ErrorGPU("The execution of 'gpu_set' was unsuccessful - the requested GPU ID does not exist. The following GPUs have been found: %s"%str(gpus))
    descriptor = [{'name':'id',   'type':'int',    'value':gpu_id }, ] 
    r = call_c_function( niftyrec_c.et_array_set_gpu_pointer, descriptor)
    #print "STATUS:",r.status
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'gpu_set' was unsuccessful.",r.status,'niftyrec_c.et_array_set_gpu')

def gpu_reset(): 
    """Reset the currently selected GPU. """ 
    r = call_c_function( niftyrec_c.et_array_reset_gpu, [])
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'gpu_reset' was unsuccessful.",r.status,'niftyrec_c.et_array_reset_gpu')
    #print "STATUS:",r.status

def PET_project(activity,attenuation,binning,use_gpu=0): #FIXME: in this and all other functions, replace 'binning' object with (only the required) raw variables 
    """PET projection; output projection data is compressed. """
    descriptor = [{'name':'activity',            'type':'array',  'value':activity}, 
                  {'name':'activity_size_x',     'type':'int',    'value':activity.shape[0]}, 
                  {'name':'activity_size_y',     'type':'int',    'value':activity.shape[1]}, 
                  {'name':'activity_size_z',     'type':'int',    'value':activity.shape[2]}, 

                  {'name':'attenuation',         'type':'array',  'value':attenuation}, 
                  {'name':'attenuation_size_x',  'type':'int',    'value':attenuation.shape[0]}, 
                  {'name':'attenuation_size_y',  'type':'int',    'value':attenuation.shape[1]}, 
                  {'name':'attenuation_size_z',  'type':'int',    'value':attenuation.shape[2]}, 

                  {'name':'use_gpu',             'type':'int',    'value':use_gpu}, ]
    r = call_c_function( niftyrec_c.PET_project, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_project' was unsuccessful.",r.status,'niftyrec_c.PET_project')
    return r.dictionary['projection']


def PET_backproject(projection_data,attenuation,binning,use_gpu=0): 
    """PET back-projection; input projection data is compressed. """
    descriptor = [    ] 
    r = call_c_function( niftyrec_c.PET_backproject, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_backproject' was unsuccessful.",r.status,'niftyrec_c.PET_backproject')
    return r.dictionary 

    
def PET_project_compressed(activity, attenuation, offsets, locations, active, 
N_axial, N_azimuthal, angles_axial, angles_azimuthal, N_u, N_v, size_u, size_v, 
activity_size_x, activity_size_y, activity_size_z, attenuation_size_x, attenuation_size_y, attenuation_size_z, 
T_activity_x, T_activity_y, T_activity_z, R_activity_x, R_activity_y, R_activity_z, 
T_attenuation_x, T_attenuation_y, T_attenuation_z, R_attenuation_x, R_attenuation_y, R_attenuation_z, 
use_gpu, N_samples, sample_step, background, background_attenuation, truncate_negative_values,direction,block_size): 
    """PET projection; output projection data is compressed. """
    N_locations = locations.shape[1] 
    #accept attenuation=None: 
    if attenuation  is None: 
        attenuation = numpy.zeros((0,0,0))
    descriptor = [{'name':'projection',             'type':'array',   'value':None,   'dtype':float32,  'size':(N_locations) }, 
                  {'name':'activity',               'type':'array',   'value':activity}, 
                  {'name':'N_activity_x',           'type':'uint',    'value':activity.shape[0]}, 
                  {'name':'N_activity_y',           'type':'uint',    'value':activity.shape[1]}, 
                  {'name':'N_activity_z',           'type':'uint',    'value':activity.shape[2]}, 
                  {'name':'activity_size_x',        'type':'float',   'value':activity_size_x}, 
                  {'name':'activity_size_y',        'type':'float',   'value':activity_size_y}, 
                  {'name':'activity_size_z',        'type':'float',   'value':activity_size_z}, 
                  
                  {'name':'T_activity_x',           'type':'float',   'value':T_activity_x}, 
                  {'name':'T_activity_y',           'type':'float',   'value':T_activity_y}, 
                  {'name':'T_activity_z',           'type':'float',   'value':T_activity_z}, 
                  {'name':'R_activity_x',           'type':'float',   'value':R_activity_x}, 
                  {'name':'R_activity_y',           'type':'float',   'value':R_activity_y}, 
                  {'name':'R_activity_z',           'type':'float',   'value':R_activity_z}, 

                  {'name':'attenuation',            'type':'array',   'value':attenuation}, 
                  {'name':'N_attenuation_x',        'type':'uint',    'value':attenuation.shape[0]}, 
                  {'name':'N_attenuation_y',        'type':'uint',    'value':attenuation.shape[1]}, 
                  {'name':'N_attenuation_z',        'type':'uint',    'value':attenuation.shape[2]}, 
                  {'name':'attenuation_size_x',     'type':'float',   'value':attenuation_size_x}, 
                  {'name':'attenuation_size_y',     'type':'float',   'value':attenuation_size_y}, 
                  {'name':'attenuation_size_z',     'type':'float',   'value':attenuation_size_z}, 
 
                  {'name':'T_attenuation_x',        'type':'float',   'value':T_attenuation_x}, 
                  {'name':'T_attenuation_y',        'type':'float',   'value':T_attenuation_y}, 
                  {'name':'T_attenuation_z',        'type':'float',   'value':T_attenuation_z}, 
                  {'name':'R_attenuation_x',        'type':'float',   'value':R_attenuation_x}, 
                  {'name':'R_attenuation_y',        'type':'float',   'value':R_attenuation_y}, 
                  {'name':'R_attenuation_z',        'type':'float',   'value':R_attenuation_z}, 

                  {'name':'N_axial',                'type':'uint',    'value':N_axial}, 
                  {'name':'N_azimuthal',            'type':'uint',    'value':N_azimuthal}, 
                  {'name':'angles_axial',           'type':'array',   'value':angles_axial}, 
                  {'name':'angles_azimuthal',       'type':'array',   'value':angles_azimuthal},
                  {'name':'N_u',                    'type':'uint',    'value':N_u}, 
                  {'name':'N_v',                    'type':'uint',    'value':N_v}, 
                  {'name':'size_u',                 'type':'float',   'value':size_u}, 
                  {'name':'size_v',                 'type':'float',   'value':size_v}, 

                  {'name':'N_locations',            'type':'uint',    'value':N_locations},
                 
                  {'name':'offsets',                'type':'array',   'value':offsets}, 
                  {'name':'locations',              'type':'array',   'value':locations}, 
                  {'name':'active',                 'type':'array',   'value':active}, 

                  {'name':'N_samples',                'type':'uint',    'value':N_samples}, 
                  {'name':'sample_step',              'type':'float',   'value':sample_step}, 
                  {'name':'background',               'type':'float',   'value':background},
                  {'name':'background_attenuation',   'type':'float',   'value':background_attenuation},  
                  {'name':'truncate_negative_values', 'type':'uint',    'value':truncate_negative_values},  

                  {'name':'use_gpu',                  'type':'uint',    'value':use_gpu}, 
                  {'name':'direction',                'type':'uint',    'value':direction}, 
                  {'name':'block_size',               'type':'uint',    'value':block_size}, 
                  ]
    r = call_c_function( niftyrec_c.PET_project_compressed, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_project_compressed' was unsuccessful.",r.status,'niftyrec_c.PET_project_compressed')
    return r.dictionary["projection"]



def PET_project_compressed_test(activity, attenuation, N_axial, N_azimuthal, offsets, locations, active): 
    N_locations = locations.shape[1]
    #accept attenuation=None: 
    if attenuation  is None: 
        attenuation = numpy.zeros((0,0,0))
    descriptor = [{'name':'projection',             'type':'array',   'value':None,   'dtype':float32,  'size':(N_locations),  }, 
                  {'name':'activity',               'type':'array',   'value':activity}, 
                  {'name':'N_activity_x',           'type':'uint',    'value':activity.shape[0]}, 
                  {'name':'N_activity_y',           'type':'uint',    'value':activity.shape[1]}, 
                  {'name':'N_activity_z',           'type':'uint',    'value':activity.shape[2]}, 

                  {'name':'attenuation',            'type':'array',   'value':attenuation}, 
                  {'name':'N_attenuation_x',        'type':'uint',    'value':attenuation.shape[0]}, 
                  {'name':'N_attenuation_y',        'type':'uint',    'value':attenuation.shape[1]}, 
                  {'name':'N_attenuation_z',        'type':'uint',    'value':attenuation.shape[2]}, 

                  {'name':'N_axial',                'type':'uint',    'value':N_axial}, 
                  {'name':'N_azimuthal',            'type':'uint',    'value':N_azimuthal}, 
                  {'name':'N_locations',            'type':'uint',    'value':N_locations},
                 
                  {'name':'offsets',                'type':'array',   'value':offsets}, 
                  {'name':'locations',              'type':'array',   'value':locations}, 
                  {'name':'active',                 'type':'array',   'value':active}, ]
    r = call_c_function( niftyrec_c.PET_project_compressed_test, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_project_compressed_test' was unsuccessful.",r.status,'niftyrec_c.PET_project_compressed_test')
    return r.dictionary["projection"]


                           

def PET_backproject_compressed(projection_data, attenuation, offsets, locations, active, 
N_axial, N_azimuthal, angles_axial, angles_azimuthal, N_u, N_v, size_u, size_v, 
N_activity_x, N_activity_y, N_activity_z, 
activity_size_x, activity_size_y, activity_size_z, 
attenuation_size_x, attenuation_size_y, attenuation_size_z, 
T_activity_x, T_activity_y, T_activity_z, R_activity_x, R_activity_y, R_activity_z, 
T_attenuation_x, T_attenuation_y, T_attenuation_z, R_attenuation_x, R_attenuation_y, R_attenuation_z, 
use_gpu, N_samples, sample_step, background, background_attenuation, direction, block_size): 
    """PET back-projection; input projection data is compressed. """
    N_locations = locations.shape[1] 
    #accept attenuation=None: 
    if attenuation  is None: 
        attenuation = numpy.zeros((0,0,0))
    descriptor = [{'name':'back_projection',        'type':'array',   'value':None,   'dtype':float32,  'size':(N_activity_x,N_activity_y,N_activity_z),  'order':"F"  }, 
                  {'name':'N_activity_x',           'type':'uint',    'value':N_activity_x}, 
                  {'name':'N_activity_y',           'type':'uint',    'value':N_activity_y}, 
                  {'name':'N_activity_z',           'type':'uint',    'value':N_activity_z}, 
                  {'name':'activity_size_x',        'type':'float',   'value':activity_size_x}, 
                  {'name':'activity_size_y',        'type':'float',   'value':activity_size_y}, 
                  {'name':'activity_size_z',        'type':'float',   'value':activity_size_z}, 
                  
                  {'name':'T_activity_x',           'type':'float',   'value':T_activity_x}, 
                  {'name':'T_activity_y',           'type':'float',   'value':T_activity_y}, 
                  {'name':'T_activity_z',           'type':'float',   'value':T_activity_z}, 
                  {'name':'R_activity_x',           'type':'float',   'value':R_activity_x}, 
                  {'name':'R_activity_y',           'type':'float',   'value':R_activity_y}, 
                  {'name':'R_activity_z',           'type':'float',   'value':R_activity_z}, 

                  {'name':'attenuation',            'type':'array',   'value':attenuation}, 
                  {'name':'N_attenuation_x',        'type':'uint',    'value':attenuation.shape[0]}, 
                  {'name':'N_attenuation_y',        'type':'uint',    'value':attenuation.shape[1]}, 
                  {'name':'N_attenuation_z',        'type':'uint',    'value':attenuation.shape[2]}, 
                  {'name':'attenuation_size_x',     'type':'float',   'value':attenuation_size_x}, 
                  {'name':'attenuation_size_y',     'type':'float',   'value':attenuation_size_y}, 
                  {'name':'attenuation_size_z',     'type':'float',   'value':attenuation_size_z}, 

                  {'name':'T_attenuation_x',        'type':'float',   'value':T_attenuation_x}, 
                  {'name':'T_attenuation_y',        'type':'float',   'value':T_attenuation_y}, 
                  {'name':'T_attenuation_z',        'type':'float',   'value':T_attenuation_z}, 
                  {'name':'R_attenuation_x',        'type':'float',   'value':R_attenuation_x}, 
                  {'name':'R_attenuation_y',        'type':'float',   'value':R_attenuation_y}, 
                  {'name':'R_attenuation_z',        'type':'float',   'value':R_attenuation_z}, 

                  {'name':'N_axial',                'type':'uint',    'value':N_axial}, 
                  {'name':'N_azimuthal',            'type':'uint',    'value':N_azimuthal}, 
                  {'name':'angles_axial',           'type':'array',   'value':angles_axial}, 
                  {'name':'angles_azimuthal',       'type':'array',   'value':angles_azimuthal},
                  {'name':'N_u',                    'type':'uint',    'value':N_u}, 
                  {'name':'N_v',                    'type':'uint',    'value':N_v}, 
                  {'name':'size_u',                 'type':'float',   'value':size_u}, 
                  {'name':'size_v',                 'type':'float',   'value':size_v}, 

                  {'name':'N_locations',            'type':'uint',    'value':N_locations},
                 
                  {'name':'offsets',                'type':'array',   'value':offsets}, 
                  {'name':'locations',              'type':'array',   'value':locations}, 
                  {'name':'active',                 'type':'array',   'value':active}, 
                  {'name':'projection_data',        'type':'array',   'value':projection_data}, 
           
                  {'name':'use_gpu',                'type':'uint',    'value':use_gpu}, 
                  {'name':'N_samples',              'type':'uint',    'value':N_samples}, 
                  {'name':'sample_step',            'type':'float',   'value':sample_step}, 
                  {'name':'background_activity',    'type':'float',   'value':background}, 
                  {'name':'background_attenuation', 'type':'float',   'value':background_attenuation}, 
                  {'name':'direction',              'type':'uint',    'value':direction},
                  {'name':'block_size',             'type':'uint',    'value':block_size},  ]
    r = call_c_function( niftyrec_c.PET_backproject_compressed, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_backproject_compressed' was unsuccessful.",r.status,'niftyrec_c.PET_backproject_compressed')
    return r.dictionary['back_projection']



def PET_initialize_compression_structure(N_axial, N_azimuthal, N_u, N_v): 
    """Obtain 'offsets' and 'locations' arrays for fully sampled PET compressed projection data. """
    descriptor = [{'name':'N_axial',              'type':'uint',         'value':N_axial         }, 
                  {'name':'N_azimuthal',          'type':'uint',         'value':N_azimuthal     }, 
                  {'name':'N_u',                  'type':'uint',         'value':N_u             },  
                  {'name':'N_v',                  'type':'uint',         'value':N_v             }, 
                  {'name':'offsets',              'type':'array',        'value':None,   'dtype':int32,     'size':(N_azimuthal,N_axial),           'order': 'F'  }, 
                  {'name':'locations',            'type':'array',        'value':None,   'dtype':uint16,    'size':(3,N_u*N_v*N_axial*N_azimuthal), 'order':'F'},
                ] 

    r = call_c_function( niftyrec_c.PET_initialize_compression_structure, descriptor )
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_initialize_compression_structure' was unsuccessful.",r.status,'niftyrec_c.PET_initialize_compression_structure')
    return [r.dictionary['offsets'],r.dictionary['locations']]


def PET_compress_projection(offsets, data, locations, N_u, N_v): 
    """Find the zero entries in fully sampled PET projection data and compress it."""
    N_locations = locations.shape[1]
    N_axial     = offsets.shape[1]
    N_azimuthal = offsets.shape[0] 
    descriptor = [  {'name':'N_locations',             'type':'uint',      'value':N_locations   }, 
                    {'name':'N_axial',                 'type':'uint',      'value':N_axial       }, 
                    {'name':'N_azimuthal',             'type':'uint',      'value':N_azimuthal   }, 
                    {'name':'N_u',                     'type':'uint',      'value':N_u           }, 
                    {'name':'N_v',                     'type':'uint',      'value':N_v           },  
                    {'name':'offsets',                 'type':'array',     'value':offsets       }, 
                    {'name':'data',                    'type':'array',     'value':data          }, 
                    {'name':'locations',               'type':'array',     'value':locations     },   
                    {'name':'projection',              'type':'array',     'value':None,   'dtype':float32,  'size':(N_locations,),   'order':'F' },  
                 ] 
    r = call_c_function( niftyrec_c.PET_compress_projection, descriptor )
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_compress_projection' was unsuccessful.",r.status,'niftyrec_c.PET_compress_projection')
    return r.dictionary['projection']


def PET_uncompress_projection(offsets, data, locations, N_u, N_v): 
    """Uncompress compressed PET projection data. """
    N_locations = locations.shape[1]
    N_axial     = offsets.shape[1]
    N_azimuthal = offsets.shape[0] 
    descriptor = [  {'name':'N_locations',             'type':'uint',      'value':N_locations   }, 
                    {'name':'N_axial',                 'type':'uint',      'value':N_axial       }, 
                    {'name':'N_azimuthal',             'type':'uint',      'value':N_azimuthal   }, 
                    {'name':'N_u',                     'type':'uint',      'value':N_u           }, 
                    {'name':'N_v',                     'type':'uint',      'value':N_v           },  
                    {'name':'offsets',                 'type':'array',     'value':offsets       }, 
                    {'name':'data',                  'type':'array',     'value':data            }, 
                    {'name':'locations',               'type':'array',     'value':locations     },   
                    {'name':'projection',              'type':'array',     'value':None,   'dtype':float32,  'size':(N_v * N_u * N_azimuthal * N_axial, ),   'order':'F' },  
                 ] 
    r = call_c_function( niftyrec_c.PET_uncompress_projection, descriptor )
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'PET_uncompress_projection' was unsuccessful.",r.status,'niftyrec_c.PET_uncompress_projection')
    return r.dictionary['projection']




def ET_spherical_phantom(voxels,size,center,radius,inner_value,outer_value): 
    """Create a spherical phantom. """
    descriptor = [{'name':'image',                 'type':'array', 'value':None,   'dtype':float32,  'size':(voxels[0],voxels[1],voxels[2]),  'order':"F"  }, 
                  {'name':'Nx',                    'type':'uint',  'value':voxels[0]}, 
                  {'name':'Ny',                    'type':'uint',  'value':voxels[1]}, 
                  {'name':'Nz',                    'type':'uint',  'value':voxels[2]}, 
                  {'name':'sizex',                 'type':'float', 'value':size[0]},
                  {'name':'sizey',                 'type':'float', 'value':size[1]},
                  {'name':'sizez',                 'type':'float', 'value':size[2]},
                  {'name':'centerx',               'type':'float', 'value':center[0]},
                  {'name':'centery',               'type':'float', 'value':center[1]},
                  {'name':'centerz',               'type':'float', 'value':center[2]},
                  {'name':'radius',                'type':'float', 'value':radius},
                  {'name':'inner_value',           'type':'float', 'value':inner_value},
                  {'name':'outer_value',           'type':'float', 'value':outer_value},  ] 
    r = call_c_function( niftyrec_c.ET_spherical_phantom, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'ET_spherical_phantom' was unsuccessful.",r.status,'niftyrec_c.ET_spherical_phantom')
    return r.dictionary['image']



def ET_cylindrical_phantom(voxels,size,center,radius,length,axis,inner_value,outer_value): 
    """Create a cylindrical phantom. """
    descriptor = [{'name':'image',                 'type':'array', 'value':None,   'dtype':float32,  'size':(voxels[0],voxels[1],voxels[2]),  'order':"F"  }, 
                  {'name':'Nx',                    'type':'uint',  'value':voxels[0]}, 
                  {'name':'Ny',                    'type':'uint',  'value':voxels[1]}, 
                  {'name':'Nz',                    'type':'uint',  'value':voxels[2]}, 
                  {'name':'sizex',                 'type':'float', 'value':size[0]},
                  {'name':'sizey',                 'type':'float', 'value':size[1]},
                  {'name':'sizez',                 'type':'float', 'value':size[2]},
                  {'name':'centerx',               'type':'float', 'value':center[0]},
                  {'name':'centery',               'type':'float', 'value':center[1]},
                  {'name':'centerz',               'type':'float', 'value':center[2]},
                  {'name':'radius',                'type':'float', 'value':radius},
                  {'name':'length',                'type':'float', 'value':length},
                  {'name':'axis',                  'type':'uint',  'value':axis},
                  {'name':'inner_value',           'type':'float', 'value':inner_value},
                  {'name':'outer_value',           'type':'float', 'value':outer_value},  ] 
    r = call_c_function( niftyrec_c.ET_cylindrical_phantom, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'ET_cylindrical_phantom' was unsuccessful.",r.status,'niftyrec_c.ET_cylindrical_phantom')
    return r.dictionary['image']



def ET_spheres_ring_phantom(voxels,size,center,ring_radius,min_sphere_radius,max_sphere_radius,N_spheres=6,inner_value=1.0,outer_value=0.0,taper=0,axis=0): 
    """Create a phantom with a ring of spheres of variable radius. """
    descriptor = [{'name':'image',                 'type':'array', 'value':None,   'dtype':float32,  'size':(voxels[0],voxels[1],voxels[2]),  'order':"F"  }, 
                  {'name':'Nx',                    'type':'uint',  'value':voxels[0]}, 
                  {'name':'Ny',                    'type':'uint',  'value':voxels[1]}, 
                  {'name':'Nz',                    'type':'uint',  'value':voxels[2]}, 
                  {'name':'sizex',                 'type':'float', 'value':size[0]},
                  {'name':'sizey',                 'type':'float', 'value':size[1]},
                  {'name':'sizez',                 'type':'float', 'value':size[2]},
                  {'name':'centerx',               'type':'float', 'value':center[0]},
                  {'name':'centery',               'type':'float', 'value':center[1]},
                  {'name':'centerz',               'type':'float', 'value':center[2]},
                  {'name':'ring_radius',           'type':'float', 'value':ring_radius}, 
                  {'name':'min_sphere_radius',     'type':'float', 'value':min_sphere_radius}, 
                  {'name':'max_sphere_radius',     'type':'float', 'value':max_sphere_radius}, 
                  {'name':'N_spheres',             'type':'uint',  'value':N_spheres}, 
                  {'name':'inner_value',           'type':'float', 'value':inner_value},
                  {'name':'outer_value',           'type':'float', 'value':outer_value},  
                  {'name':'taper',                 'type':'float', 'value':taper},  
                  {'name':'ring_axis',             'type':'uint',  'value':axis}, ] 
    r = call_c_function( niftyrec_c.ET_spheres_ring_phantom, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'ET_spheres_ring_phantom' was unsuccessful.",r.status,'niftyrec_c.ET_spheres_ring_phantom')
    return r.dictionary['image']

    

def SPECT_project_parallelholes(activity,cameras,attenuation=None,psf=None,background=0.0, background_attenuation=0.0, use_gpu=1, truncate_negative_values=0): 
    """SPECT projection; parallel-holes geometry. """
    #accept attenuation=None and psf=None: 
    if attenuation  is None: 
        attenuation = numpy.zeros((0,0,0)) 
    if psf is None: 
        psf = numpy.zeros((0,0,0)) 
    N_projections = cameras.shape[0]
    descriptor = [{'name':'activity',               'type':'array',   'value':activity }, 
                  {'name':'activity_size',          'type':'array',   'value':int32(activity.shape) }, 
                  {'name':'projection',             'type':'array',   'value':None,   'dtype':float32,  'size':(activity.shape[0],activity.shape[1],N_projections),  'order':"F"  }, 
                  {'name':'projection_size',        'type':'array',   'value':int32([N_projections, activity.shape[0], activity.shape[1]]) }, 
                  {'name':'cameras',                'type':'array',   'value':cameras,                  'order':"F" }, 
                  {'name':'cameras_size',           'type':'array',   'value':int32(cameras.shape) }, 
                  {'name':'psf',                    'type':'array',   'value':psf,                      'order':"F" }, 
                  {'name':'psf_size',               'type':'array',   'value':int32(psf.shape) }, 
                  {'name':'attenuation',            'type':'array',   'value':attenuation }, 
                  {'name':'attenuation_size',       'type':'array',   'value':int32(attenuation.shape) }, 
                  {'name':'background',             'type':'float',   'value':background }, 
                  {'name':'background_attenuation', 'type':'float',   'value':background_attenuation }, 
                  {'name':'use_gpu',                'type':'int',     'value':use_gpu }, 
                  {'name':'truncate_negative_values','type':'int',    'value':truncate_negative_values },  ]

    r = call_c_function( niftyrec_c.SPECT_project_parallelholes, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'SPECT_project_parallelholes' was unsuccessful.",r.status,'niftyrec_c.SPECT_project_parallelholes')
    return r.dictionary['projection']


def SPECT_backproject_parallelholes(projection, cameras, attenuation=None,psf=None,background=0.0, background_attenuation=0.0, use_gpu=1, truncate_negative_values=0): 
    """SPECT backprojection; parallel-holes geometry. """
    #accept attenuation=None and psf=None: 
    if attenuation  is None: 
        attenuation = numpy.zeros((0,0,0)) 
    if psf is None: 
        psf = numpy.zeros((0,0,0)) 
    N_projections = cameras.shape[0]
    descriptor = [{'name':'projection',             'type':'array',   'value':projection }, 
                  {'name':'projection_size',        'type':'array',   'value':int32(projection.shape) }, 
                  {'name':'backprojection',         'type':'array',   'value':None,   'dtype':float32,  'size':(projection.shape[0],projection.shape[1],projection.shape[0]),  'order':"F"  }, 
                  {'name':'backprojection_size',    'type':'array',   'value':int32( [projection.shape[0],projection.shape[1],projection.shape[0]] ) }, 
                  {'name':'cameras',                'type':'array',   'value':cameras,                  'order':"F" }, 
                  {'name':'cameras_size',           'type':'array',   'value':int32(cameras.shape) }, 
                  {'name':'psf',                    'type':'array',   'value':psf,                      'order':"F" }, 
                  {'name':'psf_size',               'type':'array',   'value':int32(psf.shape) }, 
                  {'name':'attenuation',            'type':'array',   'value':attenuation }, 
                  {'name':'attenuation_size',       'type':'array',   'value':int32(attenuation.shape) }, 
                  {'name':'background',             'type':'float',   'value':background }, 
                  {'name':'background_attenuation', 'type':'float',   'value':background_attenuation }, 
                  {'name':'use_gpu',                'type':'int',     'value':use_gpu }, 
                  {'name':'truncate_negative_values','type':'int',    'value':truncate_negative_values },  ]

    r = call_c_function( niftyrec_c.SPECT_backproject_parallelholes, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'SPECT_backproject_parallelholes' was unsuccessful.",r.status,'niftyrec_c.SPECT_backproject_parallelholes')
    return r.dictionary['backprojection']
    


def CT_project_conebeam(attenuation,camera_trajectory,source_trajectory,use_gpu=0): 
    """Transmission imaging projection; cone-beam geometry. """
    descriptor = [    ] 
    r = call_c_function( niftyrec_c.CT_project_conebeam, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'CT_project_conebeam' was unsuccessful.",r.status,'niftyrec_c.CT_project_conebeam')
    return r.dictionary     


def CT_backproject_conebeam(projection_data,camera_trajectory,source_trajectory,use_gpu=0): 
    """Transmission imaging back-projection; cone-beam geometry. """
    descriptor = [    ] 
    r = call_c_function( niftyrec_c.CT_backproject_conebeam, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'CT_backproject_conebeam' was unsuccessful.",r.status,'niftyrec_c.CT_backproject_conebeam')
    return r.dictionary             


def CT_project_parallelbeam(attenuation,camera_trajectory,source_trajectory,use_gpu=0): 
    """Transmission imaging projection; parallel-beam geometry. """
    descriptor = [    ] 
    r = call_c_function( niftyrec_c.CT_project_parallelbeam, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'CT_project_parallelbeam' was unsuccessful.",r.status,'niftyrec_c.CT_project_parallelbeam')
    return r.dictionary         


def CT_backproject_parallelbeam(attenuation,camera_trajectory,source_trajectory,use_gpu=0): 
    """Transmission imaging back-projection; parallel-beam geometry. """
    descriptor = [    ] 
    r = call_c_function( niftyrec_c.CT_backproject_parallelbeam, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'CT_backproject_parallelbeam' was unsuccessful.",r.status,'niftyrec_c.CT_backproject_parallelbeam')
    return r.dictionary         








def TR_grid_from_box_and_affine(box_min, box_max, box_n, affine_box2grid=None): 
    """Create 3D grid from box and affine transformation. """
    if affine_box2grid is None: 
        affine_box2grid=numpy.eye(4,dtype=float32) 
    descriptor = [{'name':'grid',                'type':'array',  'value':None,                       'dtype':float32,   'size':(box_n[0],box_n[1],box_n[2],3),   'order':"F"    }, 
                  {'name':'affine_box2grid',     'type':'array',  'value':affine_box2grid, },
                  {'name':'box_min_x',           'type':'float',  'value':box_min[0], }, 
                  {'name':'box_min_y',           'type':'float',  'value':box_min[1], }, 
                  {'name':'box_min_z',           'type':'float',  'value':box_min[2], }, 
                  {'name':'box_max_x',           'type':'float',  'value':box_max[0], }, 
                  {'name':'box_max_y',           'type':'float',  'value':box_max[1], }, 
                  {'name':'box_max_z',           'type':'float',  'value':box_max[2], }, 
                  {'name':'box_n_x',             'type':'uint',   'value':box_n[0],   }, 
                  {'name':'box_n_y',             'type':'uint',   'value':box_n[1],   }, 
                  {'name':'box_n_z',             'type':'uint',   'value':box_n[2],   }, ]
    r = call_c_function( niftyrec_c.TR_grid_from_box_and_affine, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'TR_grid_from_box_and_affine' was unsuccessful.",r.status,'niftyrec_c.TR_grid_from_box_and_affine')
    return r.dictionary['grid']  


def TR_resample_grid(image_array, grid_array, affine_index2grid=None, background=0.0, use_gpu=1, interpolation_mode=INTERPOLATION_LINEAR): 
    """Resample the image at locations specified by grid_array (array of 3D locations) and given the affine transformation that maps 
    image array indexes to grid coordinates.  """
    if affine_index2grid is None: 
        affine_index2grid=numpy.eye(4,dtype=float32) 
    resampled_shape = numpy.asarray([grid_array.shape[0],grid_array.shape[1],grid_array.shape[2]])
    descriptor = [{'name':'resampled_array',     'type':'array',  'value':None,                       'dtype':float32,   'size':(resampled_shape[0],resampled_shape[1],resampled_shape[2]),   'order':"F"    }, 
                  {'name':'image_array',         'type':'array',  'value':image_array,                'dtype':float32}, 
                  {'name':'affine',              'type':'array',  'value':affine_index2grid,          'dtype':float32}, 
                  {'name':'grid_array',          'type':'array',  'value':grid_array,                 'dtype':float32},     
                  {'name':'Nx',                  'type':'uint',   'value':image_array.shape[0]}, 
                  {'name':'Ny',                  'type':'uint',   'value':image_array.shape[1]}, 
                  {'name':'Nz',                  'type':'uint',   'value':image_array.shape[2]}, 
                  {'name':'Nx_grid',             'type':'uint',   'value':resampled_shape[0]}, 
                  {'name':'Ny_grid',             'type':'uint',   'value':resampled_shape[1]}, 
                  {'name':'Nz_grid',             'type':'uint',   'value':resampled_shape[2]}, 
                  {'name':'background',          'type':'float',  'value':background}, 
                  {'name':'use_gpu',             'type':'uint',   'value':numpy.uint32(use_gpu)},
                  {'name':'interpolation_mode',  'type':'uint',   'value':numpy.uint32(interpolation_mode)},
                   ]
    r = call_c_function( niftyrec_c.TR_resample_grid, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'TR_resample_grid' was unsuccessful.",r.status,'niftyrec_c.TR_resample_grid')
    return r.dictionary['resampled_array']
    

def TR_resample_box(image_array, box_min, box_max, box_n, affine_index2grid=None, background=0, use_gpu=1, interpolation_mode=INTERPOLATION_LINEAR): 
    pass 
    
def TR_gradient_grid(image_array, grid_array, affine_index2grid=None, background=0, use_gpu=1, interpolation_mode=INTERPOLATION_LINEAR):
    pass 
    
def TR_gradient_box(image_array, box_min, box_max, box_n, affine_index2grid=None, background=0, use_gpu=1, interpolation_mode=INTERPOLATION_LINEAR):
    pass 


def TR_transform_grid(grid_array, affine_from_grid, use_gpu=1): 
    """Transform 3D grid according to affine transformation. """
    if affine_from_grid  is None: 
        affine_from_grid = numpy.eye(4,dtype=float32) 
    descriptor = [{'name':'transformed_array',   'type':'array',  'value':None,                       'dtype':float32,   'size':(grid_array.shape[0],grid_array.shape[1],grid_array.shape[2],3),   'order':"F"    }, 
                  {'name':'grid_array',          'type':'array',  'value':grid_array,                 'dtype':float32},     
                  {'name':'Nx',                  'type':'uint',   'value':grid_array.shape[0]}, 
                  {'name':'Ny',                  'type':'uint',   'value':grid_array.shape[1]}, 
                  {'name':'Nz',                  'type':'uint',   'value':grid_array.shape[2]}, 
                  {'name':'affine',              'type':'array',  'value':affine_from_grid,           'dtype':float32}, 
                  {'name':'use_gpu',             'type':'uint',   'value':numpy.uint32(use_gpu)}, ]
    r = call_c_function( niftyrec_c.TR_transform_grid, descriptor ) 
    if not r.status == status_success(): 
        raise ErrorInCFunction("The execution of 'TR_transform_grid' was unsuccessful.",r.status,'niftyrec_c.TR_transform_grid')
    return r.dictionary['transformed_array'] 

