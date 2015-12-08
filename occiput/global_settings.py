
# occiput - Inference Language 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging, Boston
# Dec 2013, Boston, MA, USA 




# Print with 3 levels of verbosity 

__verbose = 1
def set_verbose_high(): 
    """Print everything - DEBUG mode"""
    global __verbose; __verbose = 2    
def set_verbose_medium(): 
    """Print runtime information"""
    global __verbose; __verbose = 1
def set_verbose_low(): 
    """Print only important messages"""
    global __verbose; __verbose = 0
def set_verbose_no_printing(): 
    """Do not print messages at all"""
    global __verbose; __verbose = -1        
def get_verbose_level():
    return __verbose
def print_debug(msg):
    """Use this for DEBUG Information"""
    if __verbose >= 2: 
        print msg
def print_runtime(msg):
    """Use this for messages useful at runtime"""
    if __verbose >= 1: 
        print msg
def print_important(msg):
    """Use this for important messages"""
    if __verbose >= 0: 
        print msg


# Other print options 

import contextlib as __contextlib
import numpy as __numpy
@__contextlib.contextmanager
def printoptions(*args, **kwargs):
    original = __numpy.get_printoptions()
    __numpy.set_printoptions(*args, **kwargs)
    yield 
    __numpy.set_printoptions(**original)



# GPU enables / disable

__use_gpu = True 
def enable_gpu():
    global __use_gpu; __use_gpu = True
def disable_gpu():
    global __use_gpu; __use_gpu = False 
def is_gpu_enabled():
    global __use_gpu
    return __use_gpu



# Default background of images

__background = 0.0
def set_default_background(bg):
    global __background; __background = bg
def get_default_background():
    global __background; 
    return __background


