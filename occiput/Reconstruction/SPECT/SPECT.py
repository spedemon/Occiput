
# occiput 
# Stefano Pedemonte 
# April 2014 
# April 2015 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA 


from occiput.Reconstruction.SPECT import Scintillators 
from occiput.Reconstruction.SPECT import Collimators 
from occiput.DataSources.FileSources import import_nifti 

# Import NiftyPy ray-tracers
from occiput.Core.NiftyPy_wrap import SPECT_project_parallelholes, SPECT_backproject_parallelholes, has_NiftyPy
from occiput.Core import Image3D
from occiput.Core.Errors import FileNotFound, UnknownParameter, UnexpectedParameter

from occiput.Visualization import ProgressBar, svgwrite, has_svgwrite, ipy_table, has_ipy_table
from occiput.Visualization.Colors import *

# Import DisplayNode for IPython integration
from DisplayNode import DisplayNode 

# Import various other libraries 
from numpy import *
from numpy.random import randint 
from numpy import asfortranarray, asfortranarray
from scipy import optimize 
import scipy
import scipy.signal


# Default parameters 
DEFAULT_ITERATIONS  = 20
DEFAULT_SUBSET_SIZE = 32
EPS                 = 1e-9






class SPECT_Projection(): 
    """SPECT projection object. """
    def __init__(self, data): 
        self.data = data        
    
    def get_data(self):
        """Returns the raw projection data (note that is can be accessed also as self.data ). """
        return data

    def save_to_file(self, filename): 
        h5f = h5py.File(filename, 'w')
        h5f.create_dataset('data', data=self.data)
        #h5f.create_dataset('size_x', data=size_x)
        #h5f.create_dataset('size_y', data=size_y)
        h5f.close()

    def get_integral(self): 
        return self.data.sum() 

    def to_image(self,data,index=0,scale=None,absolute_scale=False): 
        from PIL import Image
        a = float32(data[:,:,index].reshape((data.shape[0],data.shape[1])))   
        if scale is None:
            a = 255.0*(a)/(a.max()+1e-12)
        else: 
            if absolute_scale: 
                a = scale*(a) 
            else: 
                a = scale*255.0*(a)/(a.max()+1e-12)
        return Image.fromarray(a).convert("RGB")
        
    def display_in_browser(self,axial=True,azimuthal=False,index=0,scale=None): 
        self.display(axial=axial,azimuthal=azimuthal,index=index,scale=scale,open_browser=True)

    def display(self,scale=None,open_browser=False): 
        data = self.data
        d = DisplayNode()
        images = []
        progress_bar = ProgressBar(height='6px', width='100%%', background_color=LIGHT_GRAY, foreground_color=GRAY) 
        if scale is not None:
            scale = scale*255.0/(data.max()+1e-12)
        else: 
            scale = 255.0/(data.max()+1e-12) 
        N_projections = self.data.shape[2]
        N_x = self.data.shape[0]
        N_y = self.data.shape[1]
        print "SPECT Projection   [N_projections: %d   N_x: %d   N_y: %d]"%(N_projections,N_x,N_y)
        for i in range( N_projections ): 
                images.append( self.to_image(data,i,scale=scale,absolute_scale=True) ) 
                progress_bar.set_percentage(i*100.0/N_projections)                         
        progress_bar.set_percentage(100.0)
        return d.display('tipix', images, open_browser)
        
    def _repr_html_(self): 
        return self.display()._repr_html_()
        
        




class SubsetGenerator():  
    def __init__(self,N_positions):
        self._N_positions = N_positions

    def new_subset(self,mode,subset_size): 
        if mode=='random': 
            return self._random_no_replacement(subset_size) 
        elif mode=='ordered': 
            raise UnexpectedParameter("'%s' subset selection mode not yet supported."%str(mode))
        else: 
            raise UnexpectedParameter("'mode' parameter %s not recognised."%str(mode))
        
    def all_active(self): 
        return ones((self._N_positions),dtype=uint32) 

    def _random_no_replacement(self,subset_size): 
        if subset_size>=self._N_positions: 
            return self.all_active() 
        M = zeros((self._N_positions),dtype=int32) 
        n = 0
        while n<subset_size:
            active = randint(self._N_positions)
            if M[active] == 0: 
                M[active] = 1
                n+=1
        return M




def deg_to_rad(deg):
    return deg*pi/180.0

def rad_to_deg(rad):
    return rad*180.0/pi
    






class SPECT_Static_Scan(object):
    def __init__(self): 
        self._name         = "Generic SPECT Scanner"     
        self._scanner_type = "SPECT" 
        self._manufacturer = "No manufacturer"
        self._version = "0.0"
        # scanner parameters are named with 'self._p_xxx' 
        self._p_gantry_angular_positions      = 180   #[adim,integer]
        self._p_gantry_angular_position_first = 0.0   #[degrees]
        self._p_gantry_angular_position_last  = 358.0 #[degrees]     
        self._subset_generator = SubsetGenerator(self._p_gantry_angular_positions)                     
        self._p_scan_time_sec = 600.0                 #[seconds]
        self._p_radius_mm     = 300.0                 #[mm]
        self._p_n_pix_x       = 128                   #[adim]
        self._p_n_pix_y       = 128                   #[adim]
        self._p_pix_size_x_mm = 2.5                   #[mm]
        self._p_pix_size_y_mm = 2.5                   #[mm]
        self.set_background_activity(0.0) 
        self.set_background_attenuation(0.0) 
        self.set_use_gpu(True) 
        self.set_truncate_negative(False)
        self.set_scintillator( Scintillators.Ideal() )
        self.set_collimator( Collimators.LEHR() ) 
        self._measurement = None 
        self._need_update_norm = True
        self._attenuation = None 

    def get_name(self):
        return self._name

    def get_type(self):
        return self._scanner_type

    def get_manufacturer(self):
        return self._manufacturer

    def get_version(self): 
        return self._version

    def _get_parameters(self):
        parameters = {}
        dic = self.__dict__
        for k in dic.keys(): 
            if k.startswith('_p_'):
                parameters[k[3:]]=dic[k]        
        return parameters 

    def get_gantry_angular_positions(self): 
        return (self._p_gantry_angular_position_first, self._p_gantry_angular_position_last, self._p_gantry_angular_positions)

    def set_gantry_angular_positions(self, first_position_deg, last_position_deg, N_positions): 
        if not ( isscalar(first_position_deg) and isscalar(last_position_deg) and isscalar(N_positions) ): 
            raise UnexpectedParameter('Expected scalar values.')
        if not isinstance(N_positions,type(1)): 
            raise UnexpectedParameter('Expected an integer value.') 
        self._p_gantry_angular_position_first = first_position_deg
        self._p_gantry_angular_position_last =  last_position_deg
        self._p_gantry_angular_positions = N_positions
        self._subset_generator = SubsetGenerator(self._p_gantry_angular_positions)

    def get_scan_time(self): 
        return self._p_scan_time_sec

    def set_scan_time(self,scan_time_sec): 
        if not isscalar(scan_time_sec): 
            raise UnexpectedParameter('Expected a scalar value.')
        self._p_scan_time_sec = scan_time_sec

    def get_radius(self): 
        return self._p_radius_mm 

    def set_radius(self,radius_mm): 
        if not isscalar(radius_mm): 
            raise UnexpectedParameter('Expected a scalar value.') 
        self._p_radius_mm = radius_mm

    def get_n_pixels(self): 
        return (self._p_n_pix_x, self._p_n_pix_y)

    def set_n_pixels(self,n_pixels_x,n_pixels_y): 
        if (not isscalar(n_pixels_x)) or (not isscalar(n_pixels_y)): 
            raise UnexpectedParameter('Expected integer scalar values.') #FIXME: make sure it is integer 
        self._p_n_pix_x = n_pixels_x
        self._p_n_pix_y = n_pixels_y
        self._need_update_norm = True

    def get_pixel_size(self): 
        return (self._p_pix_size_x_mm, self._p_pix_size_y_mm)

    def set_pixel_size(self,pixel_size_x,pixel_size_y): 
        if (not isscalar(pixel_size_x)) or (not isscalar(pixel_size_y)): 
            raise UnexpectedParameter('Expected scalar values.') 
        self._p_pix_size_x_mm = pixel_size_x
        self._p_pix_size_y_mm = pixel_size_y

    def get_scintillator(self):
        return self._scintillator 

    def set_scintillator(self,scintillator): 
        if not isinstance(scintillator,Scintillators.BaseScintillatorSPECT): 
            raise UnexpectedParameter('Expected an instance of BaseScintillatorSPECT')
        self._scintillator = scintillator 
        self.__make_psf() 
        self._need_update_norm = True

    def get_collimator(self): 
        return self._collimator 

    def set_collimator(self,collimator): 
        if not isinstance(collimator,Collimators.BaseCollimatorSPECT): 
            raise UnexpectedParameter('Expected an instance of BaseCollimatorSPECT')
        self._collimator = collimator 
        self.__make_psf() 
        self._need_update_norm = True
        
    def set_background_activity(self,value): 
        self._background_activity    = value 
        
    def get_background_activity(self,value): 
        return self._background_activity
    
    def set_background_attenuation(self,value): 
        self._background_attenuation = value 
        
    def get_background_attenuation(self,value): 
        return self._background_attenuation
        
    def set_use_gpu(self, value): 
        self._use_gpu = value 
    
    def set_truncate_negative(self,value): 
        self._truncate_negative = value 

    def get_camera_positions(self): 
        return float32(linspace(deg_to_rad(self._p_gantry_angular_position_first),deg_to_rad(self._p_gantry_angular_position_last),self._p_gantry_angular_positions).reshape((self._p_gantry_angular_positions,1))) 

    def project(self, activity, attenuation=None, cameras=None, psf=None, subsets_array=None): 
        if isinstance(activity,ndarray): 
            activity = float32(activity)
        else: 
            activity = float32(activity.data)
        if attenuation is None:
            attenuation = self._attenuation 
            if attenuation is not None: 
                if isinstance(attenuation,ndarray):
                    attenuation = float32(attenuation)
                else: 
                    attenuation = float32(attenuation.data)
        if cameras  is None: 
            cameras = self.get_camera_positions()
        # subsets: 
        if subsets_array is not None: 
            cameras=cameras[where(subsets_array)]
        if psf  is None: 
            psf=self._psf
        proj = SPECT_project_parallelholes(activity, cameras, attenuation, psf, self._background_activity, self._background_attenuation, self._use_gpu, self._truncate_negative)
        return SPECT_Projection(proj) 


    def backproject(self, projection, attenuation=None,  cameras=None, psf=None, subsets_array=None):
        if isinstance(projection,ndarray): 
            projection = float32(projection)
        else: 
            projection = float32(projection.data)
        if attenuation is None:
            attenuation = self._attenuation 
            if attenuation is not None: 
                if isinstance(attenuation,ndarray):
                    attenuation = float32(attenuation)
                else: 
                    attenuation = float32(attenuation.data)
        if cameras  is None: 
            cameras = self.get_camera_positions()
        # subsets: 
        if subsets_array is not None: 
            cameras=cameras[where(subsets_array)]
        if psf  is None: 
            psf=self._psf   
        backproj = SPECT_backproject_parallelholes(projection, cameras, attenuation, psf, self._background_activity, self._background_attenuation, self._use_gpu, self._truncate_negative)
        return Image3D(backproj)

    def scan(self,activity_Bq,scan_time_sec=None): 
        if scan_time_sec  is None: 
            scan_time_sec = self.get_scan_time() 
        sinogram = 0
        return sinogram
        
    def __make_probabilistic_graphical_model(self): 
        pass 
        from occiput.Visualization import Graph 
        self.graph = Graph({'nodes':[{'name': 'activity', 'type': 0},{'name': 'counts', 'type': 0}],'links':[{'source': 'activity', 'target': 'counts', 'type': 't1'},]})

    def __make_psf(self): 
        self._psf = None  #FIXME 

    def set_psf(self, fwhm0_mm=0.5, depth_dependence=0.0, n_pixels=5): 
        radius = self.get_radius()
        N = self._p_n_pix_x
        [pixx,pixy] = self.get_pixel_size() 
        psf = zeros([n_pixels,n_pixels,N]) 

        def gaussian(fwhm,size): 
            x = arange(0, size, 1, float)
            y = x[:,newaxis]
            x0 = y0 = (size-1) / 2.0
            return exp(-4*log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

        for i in range(N): 
            distance = radius-N/2*pixx+i*pixx
            if distance<0: 
                distance=0
            fwhm_mm = fwhm0_mm + depth_dependence*distance
            fwhm = fwhm_mm/pixy
            psf[:,:,i] = gaussian(fwhm,n_pixels) 
        self._psf = float32(psf) 

    def get_normalization(self): 
        if self._need_update_norm: 
            self._compute_normalisation() 
        return self._norm 
        
    def _compute_normalisation(self): 
        subsets_array = self._subset_generator.all_active()
        self._norm = self.backproject(ones(( self._p_n_pix_x,self._p_n_pix_y,self._p_gantry_angular_positions ),dtype=float32, order="F") ).data 
        self._need_update_norm = False 

    def set_attenuation(self, attenuation): 
        self._attenuation = attenuation 

    def load_attenuation_from_file(self,attenuation_file):
        attenuation = import_nifti(attenuation_file).data
        self.set_attenuation(attenuation)
        
    def estimate_activity(self, attenuation=None, psf=None, iterations=DEFAULT_ITERATIONS, subset_size=DEFAULT_SUBSET_SIZE, subset_mode='random', method="EM", m=32, factr=0.01, pgtol = 1e-16, maxfun = 10000, smoothing=0.0, activity=None): 
        progress_bar = ProgressBar() 
        progress_bar.set_percentage(0.1)
        if activity==None: 
            activity = ones((self._p_n_pix_x,self._p_n_pix_y,self._p_n_pix_x),dtype=float32, order="F")
        if attenuation is None:
            attenuation = self._attenuation 
            if attenuation is not None: 
                if isinstance(attenuation,ndarray):
                    attenuation = float32(attenuation)
                else: 
                    attenuation = float32(attenuation.data)
        if psf  is None: 
            psf=self._psf
        if method=="EM": 
            print "Reconstruction method: EM"
            for i in range(iterations): 
                # Subsets: 
                if subset_size is None:
                    subsets_array = None
                    subset_size = self._p_gantry_angular_positions
                elif subset_size >= self._p_gantry_angular_positions: 
                    subsets_array = None
                    subset_size = self._p_gantry_angular_positions
                else: 
                    subsets_array = self._subset_generator.new_subset(subset_mode, subset_size)
                if subsets_array is not None: 
                    proj = self.project(activity,attenuation=attenuation,psf=psf,subsets_array=subsets_array).data
                    measurement = self._measurement[:,:,where(subsets_array)].reshape((self._p_n_pix_x,self._p_n_pix_y,subset_size))
                    measurement = asfortranarray(ascontiguousarray(measurement))
                    P = (measurement+EPS)/(proj+EPS)
                    norm = self.backproject(ones(( self._p_n_pix_x,self._p_n_pix_y, subset_size),dtype=float32, order="F"), attenuation=attenuation, psf=psf, subsets_array=subsets_array).data 
                    update = (self.backproject( P, attenuation=attenuation, psf=psf, subsets_array=subsets_array).data+EPS) / (norm +EPS) 
                else: 
                    proj = self.project(activity, attenuation=attenuation, psf=psf).data
                    P = (self._measurement+EPS)/(proj+EPS)   
                    norm = self.get_normalization()  
                    update = (self.backproject( P, attenuation=attenuation, psf=psf).data+EPS) / (norm +EPS) 
                activity = activity * update #* self.get_mask().data

                progress_bar.set_percentage((i+1)*100.0/iterations) 
                #print "Iteration: %d    max act: %f    min act: %f    max proj: %f    min proj: %f    max norm: %f    min norm: %f"%(i, activity.max(), activity.min(), proj.max(), proj.min(), norm.data.max(), norm.data.min() )
            progress_bar.set_percentage(100.0)
        elif method=="LBFGS": 
            print "Reconstruction method: LBFGS-B"
            bounds = [(None,None)] * activity.size
            for i in range(0,activity.size):
                bounds[i] = (0, None)
            args = [activity.shape,smoothing]
            activity0 = float64(activity.reshape(activity.size))
           # print "SIZE ACTIVITY0: ",activity0.shape
            activity_rec,f,d = optimize.fmin_l_bfgs_b(self.get_likelihood, activity0, fprime=self.get_gradient_activity, m=m, factr=factr, pgtol=pgtol, args=args,maxfun=maxfun,iprint=0,bounds=bounds)     
            activity = float32(activity_rec.reshape(activity.shape))
            progress_bar.set_percentage(100.0)
        else: 
            raise UnexpectedParameter("Reconstruction method %s unknown"%method)
        return Image3D(activity)
            
    def get_likelihood(self,activity, activity_size, smoothing=0.0):
                """Returns the likelihood value - given the activity. This at the moment implements the Poisson likelihood only. """
                eps = 1e-16
                sinogram = self._measurement
                psf = self._psf
                attenuation = None
                sinosize = sinogram.size
                activity = activity.reshape(activity_size)
                if any(activity<0): 
                    return 0
                #if any(activity<0): 
                #    activity[activity<0]=0
                print "MIN MAX activity: ",activity.min(), activity.max()
                proj = self.project(activity).data      #FIXME: optionally use subsets 
                print "MIN MAX proj: ",proj.min(), proj.max()
                log_proj = log(proj+eps)
                print "MIN MAX log_proj: ",log_proj.min(), log_proj.max()
                log_proj[isinf(log_proj)]=0 
                d = (-proj.reshape(sinosize) + sinogram.reshape(sinosize) * log_proj.reshape(sinosize)).sum()
                print "func:",d
                return -float64(d)

    def get_gradient_activity(self,activity, activity_size, smoothing=0.0):
                """Returns the derivative of the log-likelihood with respect to the activity. At the moment, this 
                implements only the Poisson likelihood. """
                eps = 1e-16
                sinogram = self._measurement
                psf = self._psf
                attenuation = None
                sinoshape = sinogram.shape
                activity = activity.reshape(activity_size)
                if any(activity<0): 
                    activity[activity<0]=0
                proj = self.project(activity).data + eps       #FIXME: optionally use subsets 
                norm = self.get_normalization()  
                back = self.backproject(sinogram.reshape(sinoshape)/proj.reshape(sinoshape)).data+eps
                grad = (back-norm)
                #print "SIZE0: ",activity_size, proj.shape, norm.shape, grad.shape
                kernel = asarray([[0,1,0],[1,-4,1],[0,1,0]])
                a = activity.reshape([activity_size[0],activity_size[1],activity_size[2]])
                prior = asfortranarray(smoothing * scipy.signal.convolve2d(a,kernel,'same'))
                G = grad.reshape(activity.size) + prior.reshape(activity.size)
                G = G.reshape(activity_size[0]*activity_size[1]*activity_size[2])
                return -float64(G)

    def get_likelihood_log(self, activity_log, activity_size, smoothing=0.0): 
                """Returns the likelihood value - given the log of the activity. 
                This at the moment implements the Poisson likelihood only. This is useful to optimize the likelihood 
                or the posterior with respect to the log of the activity - e.g. to include the non negativity constraint. """
                return self.get_likelihood(exp(activity_log), activity_size, smoothing)

    def get_gradient_log_activity(self, activity_log, activity_size, smoothing=0.0):
                """Returns the derivative of the log-likelihood with respect to the log of the activity. At the moment, this 
                implements only the Poisson likelihood. This is useful to optimize the likelihood or the posterior with respect to the 
                log of the activity - e.g. to include the non negativity constraint. """
                activity = exp(activity_log)
                g = self.get_gradient_activity(activity, activity_size, smoothing)
                return float64(g*activity)

    def volume_render(self,volume,scale=1.0): 
        # FIXME: use the VolumeRenderer object in occiput.Visualization (improve it), the following is a quick fix: 
        if isinstance(volume,ndarray): 
            volume = float32(volume)
        else: 
            volume = float32(volume.data)
        proj = self.project(volume).data
        proj[where(proj>proj.max()/scale )]=proj.max()/scale
        return SPECT_Projection(proj)

    def load_measurement_file(self,filename): 
        pass 

    def set_measurement(self, measurement): 
        if not ( self._p_n_pix_x == measurement.shape[0] and self._p_n_pix_y == measurement.shape[1] and self._p_gantry_angular_positions == measurement.shape[2] ):
            raise UnexpectedParameter('Measurement size is not compatible with n_pix_x, n_pix_y, gantry_angular_positions. ')
        self._measurement = measurement 

    def load_measurement_from_file(self,measurement_file):
        measurement = import_nifti(measurement_file).data
        self.set_measurement(measurement)
        
    def get_measurement(self): 
        return Volume(self._measurement)

    def display_measurement(self): 
        return SPECT_Projection(self._measurement)

    def _make_svg(self): 
        if not has_svgwrite: 
            self._svg_string = None 
            return self._svg_string 

        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('SPECT.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)

        # DETECTOR 
        # collimator 
        rect = dwg.add(dwg.rect(insert=(12, 30), size=(8, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.5).stroke('black',width=0.3,opacity=0.001)

        # scintillator
        rect = dwg.add(dwg.rect(insert=(9, 30), size=(3, 40), rx=0.5, ry=0.5))
        rect.fill('green',opacity=0.1).stroke('none',width=0.3,opacity=0.001)

        # photomultipliers 
        for i in range(8): 
            rect = dwg.add(dwg.rect(insert=(1, 31.2+i*4.8), size=(8, 4), rx=0.3, ry=0.3))
            rect.fill('grey',opacity=0.25).stroke('none',width=0.3,opacity=0.001)
        
        # IMAGING VOLUME
        rect = dwg.add(dwg.rect(insert=(30, 30), size=(40, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.02).stroke('grey',width=0.3,opacity=0.02)
        
        # GEOMETRIC NOTATIONS 
        # circle, gantry rotation 
        circle = dwg.add(dwg.circle(center=(50, 50), r=30))
        circle.fill('none').stroke('grey', width=0.1).dasharray([0.5, 0.5]) 
        # center
        circle = dwg.add(dwg.circle(center=(50, 50), r=0.5))
        circle.fill('grey',opacity=0.1).stroke('grey', width=0.1)    
        line = dwg.add(dwg.line(start=(50-1,50), end=(50+1,50)))
        line.stroke('grey', width=0.1) 
        line = dwg.add(dwg.line(start=(50,50-1), end=(50,50+1)))
        line.stroke('grey', width=0.1) 
        
        #line = dwg.add(dwg.polyline([(10, 10), (10, 100), (100, 100), (100, 10), (10, 10)],stroke='black', fill='none'))
        self._svg_string = dwg.tostring() 
        return self._svg_string

    def _repr_svg_(self): 
        self._make_svg()
        return self._svg_string    





class Gantry(): 
    def __init__(self): 
        self.svg_string = self.make_svg()

    def make_svg(self):
        if not has_svgwrite: 
            self._svg_string = None 
            return self._svg_string 

        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('test.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)

        # DETECTOR 
        # collimator 
        rect = dwg.add(dwg.rect(insert=(12, 30), size=(8, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.5).stroke('black',width=0.3,opacity=0.001)

        # scintillator
        rect = dwg.add(dwg.rect(insert=(9, 30), size=(3, 40), rx=0.5, ry=0.5))
        rect.fill('green',opacity=0.1).stroke('none',width=0.3,opacity=0.001)

        # photomultipliers 
        for i in range(8): 
            rect = dwg.add(dwg.rect(insert=(1, 31.2+i*4.8), size=(8, 4), rx=0.3, ry=0.3))
            rect.fill('grey',opacity=0.25).stroke('none',width=0.3,opacity=0.001)
        
        # IMAGING VOLUME
        rect = dwg.add(dwg.rect(insert=(30, 30), size=(40, 40), rx=0.5, ry=0.5))
        rect.fill('grey',opacity=0.02).stroke('grey',width=0.3,opacity=0.02)
        
        # GEOMETRIC NOTATIONS 
        # circle, gantry rotation 
        circle = dwg.add(dwg.circle(center=(50, 50), r=30))
        circle.fill('none').stroke('grey', width=0.1).dasharray([0.5, 0.5]) 
        # center
        circle = dwg.add(dwg.circle(center=(50, 50), r=0.5))
        circle.fill('grey',opacity=0.1).stroke('grey', width=0.1)    
        line = dwg.add(dwg.line(start=(50-1,50), end=(50+1,50)))
        line.stroke('grey', width=0.1) 
        line = dwg.add(dwg.line(start=(50,50-1), end=(50,50+1)))
        line.stroke('grey', width=0.1) 
        
        #line = dwg.add(dwg.polyline([(10, 10), (10, 100), (100, 100), (100, 10), (10, 10)],stroke='black', fill='none'))
        return dwg.tostring() 

    def _repr_svg_(self): 
        return self.svg_string 




class GE_Infinia(SPECT_Static_Scan):
    def __init__(self): 
        SPECT.__init__(self)
        self._name = "GE Infinia SPECT Scanner with LEHR collimator"  
