
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Jan. 2015, Boston
# Feb. 2015, Helsinki
# Nov. 2015, Boston 


# This file defines a PET projection: the data structure that contains projection data. 
# In occiput.io, the projection data is of two types: sparse or non-sparse. 
# The sparse projection data is memory efficient and leads to efficient projection and 
# back-projection. Low memory is necessary in order to use hundreds or thousands of sinograms 
# for dynamic and kinetic imaging. There is a single object to represent projection data, 
# with a flag that sets the type: sparse or non-sparse. The ray-tracer knows how to handle the projection data object.  


import occiput 
from occiput.Visualization import ipy_table, has_ipy_table, svgwrite, has_svgwrite 
from occiput.Core.Print import rad_to_deg, deg_to_rad, array_to_string
from occiput.Core.NiftyPy_wrap import PET_compress_projection, PET_uncompress_projection, PET_initialize_compression_structure #, PET_compress_projection_array 

import h5py 
import copy
import os
import inspect 
from numpy import isscalar, linspace, int32, uint32, uint16, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan, tile
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 
from numpy.random import poisson



__all__ = ["DEFAULT_BINNING","Binning", "PET_Projection_Sparsity", "PET_Projection",
          "PET_compress_projection", "PET_uncompress_projection", "PET_initialize_compression_structure","display_PET_Projection_geometry",] #"PET_compress_projection_array"] 

    
def display_PET_Projection_geometry(): 
    filename = inspect.getfile(occiput).strip("__init__.pyc")+"Data"+os.sep+"occiput_ray_tracer_PET.png"
    from IPython.display import Image
    img = Image(filename=filename, width=1000) 
    return img 



EPS = 1e-10
DEFAULT_BINNING = {      "n_axial":                120,
                         "n_azimuthal":            5,
                         "angles_axial":           float32( linspace(0,pi-pi/120.0,120) ), 
                         "angles_azimuthal":       float32( [-0.5, -0.25, 0.0, 0.25, 0.5] ),
                         "size_u":                 256.0,
                         "size_v":                 256.00,
                         "n_u":                    128,
                         "n_v":                    128,   }



class Binning(): 
    """PET detectors binning. """
    def __init__(self,parameters=None): 
        self._initialised = 0
        self.name = "Unknown binning name"
        if parameters is None: 
            self.load_from_dictionary(DEFAULT_BINNING)  
        elif type(parameters) == dict:  
            self.load_from_dictionary(parameters)        
        else: 
            raise UnknownParameter('Parameter %s specified for the construction of Binning is not compatible. '%str(parameters)) 
            
    def load_from_dictionary(self,dictionary):
        self.N_axial                = dictionary['n_axial']                    # Number of axial bins 
        self.N_azimuthal            = dictionary['n_azimuthal']                # Number of azimuthal bins
        self.angles_axial           = dictionary['angles_axial']               # Axial angles radians
        self.angles_azimuthal       = dictionary['angles_azimuthal']           # Azimuthal angles radians
        self.size_u                 = dictionary['size_u']                     # Size of the detector plane, axis u, [adimensional]
        self.size_v                 = dictionary['size_v']                     # Size of the detector plane, axis v,  [adimensional]
        self.N_u                    = dictionary['n_u']                        # Number of pixels of the detector plan, axis u 
        self.N_v                    = dictionary['n_v']                        # Number of pixels of the detector plan, axis v 
        self.angles_axial           = dictionary['angles_axial']    
        
    def get_angles(self, subsets_matrix=None): 
        angles = zeros((2,self.N_azimuthal,self.N_axial))
        angles[0,:,:] = tile(self.angles_axial,(self.N_azimuthal,1))
        angles[1,:,:] = tile(self.angles_azimuthal,(self.N_axial,1)).T
        
        if subsets_matrix is not None: 
            subsets_matrix = int32(subsets_matrix)
            N = subsets_matrix.sum()
            angles2 = angles.reshape(2,self.N_azimuthal*self.N_axial).copy()
            angles = zeros([2, 1, N])
            angles[0,0,:] = ( angles2[0,subsets_matrix.flatten()==1] ).reshape([N,])
            angles[1,0,:] = ( angles2[1,subsets_matrix.flatten()==1] ).reshape([N,])
        return float32(angles)

    def __repr__(self): 
        s = "PET Binning: \n"        
        s = s+" - N_axial_bins:             %d \n"%self.N_axial 
        s = s+" - N_azimuthal_bins:         %d \n"%self.N_azimuthal 
        s = s+" - Size_u:                   %f \n"%self.size_u
        s = s+" - Size_v:                   %f \n"%self.size_v        
        s = s+" - N_u:                      %d \n"%self.N_u
        s = s+" - N_v:                      %d \n"%self.N_v
        s = s+" - Angles axial:             %s \n"%array_to_string(self.angles_axial)
        s = s+" - Angles azimuthal:         %s \n"%array_to_string(self.angles_azimuthal)
        return s

    def _repr_html_(self):       #FIXME: add angles_axial and angles_azimuthal 
        if not has_ipy_table: 
            return "Please install ipy_table."
        table_data = [['N_axial',self.N_axial],['N_azimuthal',self.N_azimuthal],['Size_u',self.size_u],['Size_v',self.size_v],['N_u',self.N_u],['N_v',self.N_v]] 
        N_per_row = 12
        for i in range(len(self.angles_azimuthal)/N_per_row+1): 
            i_start = i*N_per_row
            i_end = (i+1)*N_per_row 
            if i_end>=len(self.angles_azimuthal): 
                i_end=len(self.angles_azimuthal)
            if i==0: 
                table_data.append([ "Angles_azimuthal [deg]", array_to_string(rad_to_deg(self.angles_azimuthal[i_start:i_end]),"%3.4f ") ]) 
            else: 
                table_data.append([ "", array_to_string(rad_to_deg(self.angles_azimuthal[i_start:i_end]),"%3.4f ") ])
        i_start=0
        i_end  =len(self.angles_axial)-1
        if i_end-i_start >= 4: 
            table_data.append([ "Angles_axial [deg]", str(rad_to_deg(self.angles_axial[i_start]))+" "+str(rad_to_deg(self.angles_axial[i_start+1]))+"  ...  "+str(rad_to_deg(self.angles_axial[i_end-1]))+"   "+str(rad_to_deg(self.angles_axial[i_end])), ])
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic_left')
        #table = ipy_table.set_column_style(0, color='lightBlue')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        return table._repr_html_()
        
    def _make_svg(self): 
        if not has_svgwrite: 
            self._svg_string = None
            return self._svg_string

        w = '100%'
        h = '100%'
        
        dwg = svgwrite.Drawing('test.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)
        
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

    def display_geometry(self):
        return display_PET_Projection_geometry()
        


class PET_Projection_Sparsity(): 
    """Represents the sparsity pattern of projection data and implements methods to compress and uncompress data according to the sparsity pattern of self. 
    If offsets and locations are not specified as parameters of the constructor, the data structure represents dense projection data. """
    def __init__(self, N_axial, N_azimuthal, N_u, N_v , offsets=None, locations=None): 
        if offsets is None or locations is None: 
            # Uncompressed 
            [offsets,locations] = PET_initialize_compression_structure(N_axial,N_azimuthal,N_u,N_v) 
            offsets = offsets.transpose().reshape((N_azimuthal, N_axial))
        self.offsets     = offsets
        self.locations   = locations
        self.N_axial     = int32(N_axial)
        self.N_azimuthal = int32(N_azimuthal)
        self.N_u         = N_u
        self.N_v         = N_v

    def compress_data(self, data): 
        """Compresses an array of uncompressed projection data, according to the sparsity pattern of self (self.offsets and self.locations). """
        compressed_data = PET_compress_projection(self.offsets, data, self.locations, self.N_u, self.N_v)
        #print "Compression done"
        return compressed_data

    def uncompress_data(self, data): 
        """Uncompresses a data array that contains projection compressed according to the sparsity pattern of self (self.offsets and self.locations). """
        uncompressed_data = PET_uncompress_projection(self.offsets, data, self.locations, self.N_u, self.N_v)
        uncompressed_data.reshape([self.offsets.shape[1], self.offsets.shape[0], self.N_u, self.N_v]) 
        #print "Uncompression done"
        return uncompressed_data

    def same_sparsity_as(self, obj): 
        """Returns True if the argument is an object of type PET_Projection_Sparsity with the same identical sparsity pattern as self. """
        if not hasattr(obj,'offsets'): 
            return False 
        if not hasattr(obj,'locations'): 
            return False 
        return obj.offsets == self.offsets and obj.locations == self.locations

    def is_compressed(self): 
        return not self.is_uncompressed() 

    def is_uncompressed(self): 
        """Returns True if the sparsity pattern represents an uncompressed sampled projection. """
        return self.get_N_locations() == self.N_u * self.N_v * self.N_axial * self.N_azimuthal

    def get_N_locations(self): 
        """Returns the number of non-zero entries of the sparse projection data. """
        return self.locations.shape[1]

    # Overload math operators 
    def __and__(self, other): 
        """Overload the AND operator. """
        if not isinstance(other,self.__class__): 
            raise Exception("Binary operation must be with another object of the same type - PET_Projection_Sparsity. ")
        print "This is not implemented, it should be implemented at low level. "

    def __xor__(self, other): 
        if not isinstance(other,self.__class__): 
            raise Exception("Binary operation must be with another object of the same type - PET_Projection_Sparsity. ")
        print "This is not implemented, it should be implemented at low level. "
    
    def __or__(self, other): 
        if not isinstance(other,self.__class__): 
            raise Exception("Binary operation must be with another object of the same type - PET_Projection_Sparsity. ")
        print "This is not implemented, it should be implemented at low level. "

    def get_locations_per_plane(self):
        """Returns the number of active locations per projection plane."""
        shape = self.offsets.shape
        offsets = self.offsets.transpose().flatten()
        v = offsets.copy()
        v[0:-1] = offsets[1::]
        v[-1]   = self.locations.shape[1]
        sizes = v - offsets 
        return sizes.reshape(shape)

    def get_subset(self, subset_matrix):
        subset_matrix = int32(subset_matrix)
        indexes = subset_matrix == 1
        N = int32(subset_matrix.sum())
        locations = self.locations
        offsets   = self.offsets
        N_axial = self.N_axial
        N_azimuthal = self.N_azimuthal
        N_u = self.N_u
        N_v = self.N_v

        offsets_return = zeros((1,N), dtype=int32, order="F")
        sizes = self.get_locations_per_plane().reshape(indexes.shape)
        
        N_locations = sizes[indexes].sum()
        locations_return = zeros((3,N_locations), dtype=uint16, order="F")
    
        ii = 0
        previous_offset = 0
        for iax in range(N_axial):
            for iaz in range(N_azimuthal): 
                if subset_matrix[iaz,iax]: 
                    ii = ii + 1
                    size = sizes[iaz,iax]
                    if ii<N: 
                        offsets_return[0,ii] = offsets_return[0,ii-1]+size
                    locations_return[:,offsets_return[0,ii-1]:offsets_return[0,ii-1]+size] = locations[:,offsets[iaz,iax]:offsets[iaz,iax] + size]
        return PET_Projection_Sparsity(N, 1, N_u, N_v , offsets = offsets_return, locations = locations_return)

    def display_geometry(self):
        return display_PET_Projection_geometry()


class PET_Projection(): 
    """Sparse PET projection object. If offsets=None or locations=None in the constructor, the projection is uncompressed. """
    def __init__(self, binning, data=None, offsets=None, locations=None, time_bins=None, subsets_matrix=None): 
        """Sparse PET projection object. If offsets=None or locations=None in the constructor, the projection is uncompressed. """
        # The data structure that represents a projection is composed of: 
        #1) Sparsity structure 
        #2) Time binning (list of time flags)
        #3) Size of projection planes - in units of [mm]
        #4) Angular steps in axial and azimuthal directions - in units of [radians] 
        #5) The projection data (e.g. counts)
        self.binning = binning 
        self.angles = self.binning.get_angles(subsets_matrix=subsets_matrix)
        
        sparsity = PET_Projection_Sparsity(binning.N_axial, binning.N_azimuthal, binning.N_u, binning.N_v, offsets, locations)
        if subsets_matrix is not None: 
            self.sparsity = sparsity.get_subset(subsets_matrix)
        else: 
            self.sparsity = sparsity

        if time_bins is None: 
            time_bins = int32([0,1000])   # 1 second if not specified 
        self.time_bins = time_bins
        
        N_axial = int32(self.angles.shape[2])
        N_azimuthal = int32(self.angles.shape[1])
        
        if data is None: 
            data = zeros((N_axial, N_azimuthal, binning.N_u, binning.N_v))
        elif isscalar(data): 
            data = data + zeros((N_axial, N_azimuthal, binning.N_u, binning.N_v))
        if self.sparsity.is_uncompressed(): 
            data = data.reshape((N_axial, N_azimuthal, binning.N_u, binning.N_v))
        self.data = data     
          
    def get_binning(self): 
        """Returns the binning parameters (that is an instance of Binning). """
        return self.binning 

    def get_time_bins(self):
        """Returns the time bins (list of N+1 time flags in case of N time bins) associated to the projection. """
        return self.time_bins 
  
    def get_sparsity(self): 
        """Returns an instance of PET_Projection_Sparsity, which describes the sparsity pattern of the projection data. """
        return self.sparsity 
    
    def get_data(self):
        """Returns the (compressed) raw projection data (note that is can be accessed also as self.data ). """
        return self.data

    def get_angles(self): 
        return self.angles

    def compress_as_self(self, projection):
        """Returns a PET_Projection object obtained by sampling the PET_Projection object given as argument 
        according to the sparsity of self. The argument may be a compressed or uncompressed PET_Projection object. """
        # FIXME: take a PET_Proejction_sparsity as argument optionally
        # FIXME: verify that the two PET_Projection objects have the same binning 
        if isinstance(projection, PET_Projection): 
            if not projection.is_uncompressed(): 
                print "Re-compressing projection data, possible data loss. "
                uncompressed_projection = projection.uncompress_self() 
                uncompressed_data = uncompressed_projection.data
            else: 
                uncompressed_projection = projection
                uncompressed_data = uncompressed_projection.data
        else: 
            uncompressed_data = projection
            # FIXME: verify that size of the array is correct for offsets and locations 
        compressed_data = self.sparsity.compress_data(uncompressed_data) 
        return PET_Projection( self.get_binning(), compressed_data, self.sparsity.offsets, self.sparsity.locations, self.get_time_bins() )

    def compress_as(self, projection):
        """Returns a PET_Projection object obtained by sampling the projection data of self 
        according to the sparsity of the PET_Proejction object given as argument. 
        The argument may be a compressed or uncompressed PET_Projection object. """
        # FIXME: take a PET_Proejction_sparsity as argument optionally
        # FIXME: verify that the two PET_Projection objects have the same binning 
        uncompressed_data = self.sparsity.uncompress_data(self.data) 
        data = projection.sparsity.compress_data(uncompressed_data)
        return PET_Projection( projection.get_binning(), data, projection.sparsity.offsets, projection.sparsity.locations, projection.get_time_bins() )
        
    def compress_self(self): 
        if self.is_compressed(): 
            return self
        else: 
            print "Find the zeros and compress. "
            print "Not implemented. Please implement me. This needs low level implementation. Right now sparsity comes only from listmode data. "
            #offsets, locations, data = PET_compress_projection_array(self.data, self.N_axial, self.N_azimuthal, self.N_u, self.N_v)
    #return PET_Projection( self.get_binning(), data, offsets, locations, self.get_time_bins() )

    def uncompress_self(self): 
        """Returns an instance of PET_Projection obtained by zero-filling self. The projection object returned is not sparse. """
        if self.is_uncompressed(): 
            return self
        uncompressed_data = self.sparsity.uncompress_data(self.data) 
        return PET_Projection( self.get_binning(), uncompressed_data, offsets=None, locations=None, time_bins=self.get_time_bins() )

    def save_to_file(self, filename): 
        h5f = h5py.File(filename, 'w')
        h5f.create_dataset('offsets', data=self.sparsity.offsets)
        h5f.create_dataset('locations', data=self.sparsity.locations)
        h5f.create_dataset('data', data=self.data)
        h5f.create_dataset('time_bins', data=self.time_bins)
        binning = self.get_binning() 
        h5f.create_dataset('n_axial', data=binning.N_axial)
        h5f.create_dataset('n_azimuthal', data=binning.N_azimuthal)
        h5f.create_dataset('angles_axial', data=binning.angles_axial)
        h5f.create_dataset('angles_azimuthal', data=binning.angles_azimuthal)
        h5f.create_dataset('size_u', data=binning.size_u)
        h5f.create_dataset('size_v', data=binning.size_v)
        h5f.create_dataset('n_u', data=binning.N_u)
        h5f.create_dataset('n_v', data=binning.N_v)
        h5f.close()

    def same_sparsity_as(self, obj): 
        """Returns True if the argument is an object of type PET_Projection_Sparsity with the same identical sparsity pattern as self. """
        if not hasattr(obj,'sparsity'): 
            print "PET.py - same_sparsity(): Object of wrong type"  # FIXME:raise input error
        return self.sparsity.same_sparsity_as(obj.sparsity) 

    def is_compressed(self): 
        return not self.is_uncompressed() 

    def is_uncompressed(self): 
        """Returns True if the sparsity pattern represents an uncompressed projection. """
        return self.sparsity.is_uncompressed() 

    def get_N_locations(self): 
        """Returns the number of non-zero entries of the sparse projection data. """
        return self.sparsity.get_N_locations() 

    def max(self):
        return self.data.max() 
    
    def min(self):
        return self.data.min()

    def get_time_start(self):
        return self.time_bins[0]
    
    def get_time_end(self): 
        return self.time_bins[-1] 

    def get_duration(self): 
        return self.get_time_end() - self.get_time_start() 

    def get_integral(self): 
        return self.data.sum() 

    def get_compression_ratio(self):
        return (1.0 * self.get_N_locations()) / (self.sparsity.N_u*self.sparsity.N_v*self.sparsity.N_axial*self.sparsity.N_azimuthal)
        
    def get_listmode_loss(self): 
        return (1.0 * self.get_integral()) / self.get_N_locations()

    def get_projection_info(self):
        """Returns information related to the compression, such as the compression ratio. It returns a dictionary. """
        return self.sparsity.get_projection_info() 

    def to_nd_array(self): 
        """Uncompress (if compressed) and convert into an nd array. """
        if not self.is_uncompressed(): 
            uncompressed = self.uncompress_self() 
        else: 
            uncompressed = self 
        shape = (uncompressed.sparsity.N_axial, uncompressed.sparsity.N_azimuthal, uncompressed.sparsity.N_u, uncompressed.sparsity.N_v ) 
        return uncompressed.data.reshape(shape)

    def _is_3D_data(self): 
        if self.sparsity.N_axial == 1 or self.sparsity.N_azimuthal == 1: 
            return True 
        return False

    def _is_4D_data(self): 
        if self.sparsity.N_axial != 1 and self.sparsity.N_azimuthal != 1: 
            return True 
        return False

    def to_image(self,data,index_axial=0,index_azimuthal=0,scale=None,absolute_scale=False): 
        from PIL import Image
        a = float32(data[index_axial,index_azimuthal,:,:].reshape((data.shape[2],data.shape[3])))   
        if scale is None:
            a = 255.0*(a)/(a.max()+1e-12)
        else: 
            if absolute_scale: 
                a = scale*(a) 
            else: 
                a = scale*255.0*(a)/(a.max()+1e-12)
        im = Image.fromarray(a).convert("RGB")
        return im.rotate(90) 
        
    def display_in_browser(self,axial=True,azimuthal=False,index=0,scale=None): 
        self.display(axial=axial,azimuthal=azimuthal,index=index,scale=scale,open_browser=True)

    def display(self,axial=True,azimuthal=None,index=0,scale=None,open_browser=False): 
        """If 'axial' is set to True, then it displays projections along the axial direction, with azimuthal index given by 'index'. 
        If 'axial' is set to False, then it displays projections along the azimuthal direction, with axial index given by 'index'. """
        data = self.to_nd_array() 
        binning = self.get_binning()
        d = DisplayNode()
        images = []
        progress_bar = ProgressBar(height='6px', width='100%%', background_color=LIGHT_GRAY, foreground_color=GRAY) 
        if scale is not None:
            scale = scale*255.0/(data.max()+1e-12)
        else: 
            scale = 255.0/(data.max()+1e-12) 
        if axial is None and azimuthal is None:   # if axial and azimuthal are not specified, select them automatically based on data shape 
            axial = (self.sparsity.N_axial > 1) 
            azimuthal = (self.sparsity.N_azimuthal > 1)
        else: 
            if axial == None: 
                axial = False
            if azimuthal == None: 
                azimuthal = False 
        if axial and not azimuthal: 
            for i in range( self.sparsity.N_axial ): 
                images.append( self.to_image(data,i,index,scale=scale,absolute_scale=True) ) 
                progress_bar.set_percentage(i*100.0/self.sparsity.N_axial) 
        elif azimuthal and not axial: 
            for i in range( binning.N_azimuthal ): 
                images.append( self.to_image(data,index,i,scale=scale,absolute_scale=True) )
                progress_bar.set_percentage(i*100.0/binning.N_azimuthal)     
        else: 
            for i in range(self.sparsity.N_axial): 
                images_az = []
                for j in range(binning.N_azimuthal): 
                    images_az.append( self.to_image(data,i,j,scale=scale,absolute_scale=True) )
                images.append(images_az) 
                progress_bar.set_percentage(i*100.0/self.sparsity.N_axial)                         
        progress_bar.set_percentage(100.0) 
        return d.display('tipix', images, open_browser) 

    def _repr_html_(self): 
        return self.display()._repr_html_() 


    def _is_same_sparsity_(self, other): 
        if isinstance(other,self.__class__): 
            # FIXME: make sure that the two projections are compatible: same geometry and sparsity pattern
            return True 
        return False 

    def copy(self): 
        return copy.deepcopy(self)

    
    def crop(self, bins_range):
        """Crop projection"""
        if self.is_compressed(): 
            print "Cropping of PET_Projection currently only works with uncompressed data. Please implement for compressed data. "
            return 
        A = bins_range[0]
        B = bins_range[1]
        data = self.data[:,:,A:B,:]
        self.binning.size_u = (1.0*self.binning.size_u) / self.binning.N_u * (B-A)
        self.binning.N_u = B-A
        return PET_Projection(self.binning, data=data, time_bins=self.time_bins) 
        
    def apply_noise_Poisson(self):
        self.data = poisson(self.data)
        
    # Overload math operators 

    def __add__(self, other): 
        """Overload the '+' (sum) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other):             
            other = other.data 
        out = self.copy()
        out.data = out.data + other
        return out

    def __sub__(self, other): 
        """Overload the '-' (subtraction) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other): 
            other = other.data 
        out = self.copy()
        out.data = out.data - other
        return out

    def __mul__(self, other): 
        """Overload the '*' (multiplication) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other): 
            other = other.data 
        out = self.copy()
        out.data = out.data * other
        return out

    def __div__(self, other): 
        """Overload the '/' (division) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other): 
            other = other.data   
        out = self.copy()
        out.data = out.data / other
        return out

    def __radd_(self, other): 
        return self.__add__(other) 

    def __rmul_(self, other): 
        return self.__add__(other)

    def __rsub__(self, other): 
        """Overload the '-' (subtraction) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other): 
            other = other.data 
        out = self.copy()
        out.data = other-out.data 
        return out

    def __rdiv__(self, other): 
        """Overload the '/' (division) operator. """
        if other is None: #This makes reconstruction code cleaner, when e.g. scatter and randoms are not specified
            return self
        if self._is_same_sparsity_(other): 
            other = other.data 
        out = self.copy()
        out.data = other/out.data 
        return out

    def get_subset(self, subsets_matrix):
        indexes = subsets_matrix.flatten() == 1
        data=self.data.swapaxes(0,1).reshape((self.sparsity.N_axial*self.sparsity.N_azimuthal,self.binning.N_u,self.binning.N_v))[indexes,:,:]
        #sparsity = self.sparsity.get_subset(subsets_matrix)
        #offsets = sparsity.offsets
        #locations = sparsity.locations 
        return PET_Projection(self.binning, data, self.sparsity.offsets, self.sparsity.locations, self.time_bins, subsets_matrix)


    def display_geometry(self):
        return display_PET_Projection_geometry()

