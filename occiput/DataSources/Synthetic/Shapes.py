
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston


#__all__ = ['uniform_sphere','uniform_cylinder','uniform_spheres_ring']



from occiput.Core import Image3D as __Image3D
try: 
    from NiftyPy.NiftyRec import ET_spherical_phantom as __ET_spherical_phantom
    from NiftyPy.NiftyRec import ET_cylindrical_phantom as __ET_cylindrical_phantom
    from NiftyPy.NiftyRec import ET_spheres_ring_phantom as __ET_spheres_ring_phantom
except: 
    has_NiftyPy = False
    print "Please install NiftyPy"
else: 
    has_NiftyPy = True

class InstallationError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

        

def uniform_sphere(voxels=[256,256,256],size=[1,1,1],center=[0.5,0.5,0.5],radius=0.2,inner_value=1.0,outer_value=0.0): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    if not has_NiftyPy: 
        raise InstallationError("Please install NiftyPy to enable function uniform_sphere(..)")
    return __Image3D(__ET_spherical_phantom(voxels,size,center,radius,inner_value,outer_value)) 

def uniform_cylinder(voxels=[256,256,256],size=[1,1,1],center=[0.5,0.5,0.5],radius=0.3,length=0.7,axis=1,inner_value=1.0,outer_value=0.0): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    if not has_NiftyPy: 
        raise InstallationError("Please install NiftyPy to enable function uniform_cylinder(..)")
    return __Image3D(__ET_cylindrical_phantom(voxels,size,center,radius,length,axis,inner_value,outer_value)) 
    
def uniform_spheres_ring(voxels=[256,256,256],size=[1,1,1],center=[0.5,0.5,0.5],ring_radius=0.2,min_sphere_radius=0.01,max_sphere_radius=0.1,N_spheres=6,inner_value=1.0,outer_value=0.0,taper=0,axis=0): 
    """Create volume (3D numpy array) with uniform value within a spherical region. """
    if not has_NiftyPy: 
        raise InstallationError("Please install NiftyPy to enable function uniform_spheres_ring(..)")
    return __Image3D(__ET_spheres_ring_phantom(voxels,size,center,ring_radius,min_sphere_radius,max_sphere_radius,N_spheres,inner_value,outer_value,taper,axis))

