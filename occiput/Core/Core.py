
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec. 2014, Boston, MA 
# April. 2014, Boston, MA


import numpy
from PIL import Image
import nibabel
import json 
import DisplayNode 
import copy

from occiput.Visualization import MultipleVolumes, MultipleVolumesNiftyPy, VolumeRenderer
from occiput.Visualization import ipy_table, has_ipy_table
from occiput.Core import transformations as tr 
from occiput import global_settings 
from occiput.global_settings import printoptions
from occiput.Core.NiftyPy_wrap import transform_grid, grid_from_box_and_affine, resample_image_on_grid, INTERPOLATION_LINEAR, INTERPOLATION_POINT
from occiput.Core.Conversion import nipy_to_occiput, nifti_to_occiput, occiput_to_nifti, occiput_from_array





class Transform_Affine(object): 
    def __init__(self,data=None,map_from="", map_to=""): 
        self.set_map_from(map_from)
        self.set_map_to(map_to)
        self.set_data(data)

    def can_inverse_left_multiply(self,obj): 
        return self.__space_of_obj(obj) == self.map_to

    def can_left_multiply(self,obj):       
        return self.__space_of_obj(obj) == self.map_from
        
    def left_multiply(self,obj): 
        # if obj is a nd-array (of size [4xN]): 
        if isinstance(obj,numpy.ndarray): 
            return numpy.dot(self.data,obj) 
        # verify if object and affine are in the same space
        if not self.can_left_multiply(obj): 
            if hasattr(obj,'map_from'): 
                X = obj.map_from
            else: 
                X = "X"
            print "Affine is incompatible with the given object: composing [%s,%s] with [%s,%s] is not allowed. "%( self.map_to,self.map_from,self.__space_of_obj(obj),X )
        data = numpy.dot(self.data,obj.data)
        # return an affine transformation is the input was an affine transformation 
        if isinstance(obj,Transform_Affine): 
            return Transform_Affine(data=data, map_from=obj.map_from, map_to=self.map_to)
        else:
            return data

    def __space_of_obj(self,obj): 
        if hasattr(obj,'space'): 
            return obj.space
        elif hasattr(obj,'map_to'): 
            return obj.map_to
        else: 
            return ""

    def export_dictionary(self): 
        return {'map_to':self.map_to, 'map_from':self.map_from, 'data':self.data.tolist()}

    def export_json(self): 
        return json.dumps(self.export_dictionary())

    def save_to_file(self,filename): 
        with open(filename, 'w') as fid: 
            json.dump(self.export_json(), fid)

    def load_from_file(self,filename): 
        with open(filename, 'r') as fid: 
            self.import_json(json.load(fid)) 

    def import_json(self,string): 
        self.import_dictionary(json.loads(string))
        
    def import_dictionary(self,dict):
        self.map_to = dict['map_to']
        self.map_from = dict['map_from']
        self.data  = numpy.asarray(dict['data'])

    def is_6DOF(self):
        #tra, rot, axis = self.to_translation_rotation() 
        #if tra is None: 
        #    return False
        #mat = numpy.dot(tr.translation_matrix(tra), tr.rotation_matrix(rot,axis)) 
        #return tr.is_same_transform(mat,self.data)
        return True # FIXME
        
    # FIXME: implement the following functions (see is_6DOF() )
    def is_9DOF(self): 
        tra, scale, rot, rot_axis = self.to_translation_scale_rotation()
        if tra is None: 
            return False   
        tra_mat = tr.translation_matrix(tra)
        scale_mat = numpy.diag([scale[0],scale[1],scale[2],1.0])
        rot_mat = tr.rotation_matrix(rot,rot_axis)
        mat = numpy.dot(tra_mat, numpy.dot(scale_mat, rot_mat ) )
        return tr.is_same_transform(mat,self.data)

    def is_rigid(self):
        #return self.is_9DOF()
        return True

    def is_rotation(self): 
        return False # FIXME

    def is_translation(self):
        return False # FIXME

    def is_scale(self): 
        return False # FIXME

    def determinant(self): 
        try: 
            det = numpy.linalg.det(self.data)
        except: 
            det = None    
        return det

    def to_translation_rotation(self): 
        tra = tr.translation_from_matrix(self.data)
        mat = self.data.copy()
        mat[0:3,3]=0
        try: 
            rot,axis,point = tr.rotation_from_matrix(mat)
        except ValueError: 
            #print "Not a rotation matrix. "
            return None, None, None
        return tra, rot, axis

    def from_translation_rotation(self, tra, rot, axis=[0,0,1]): 
        return numpy.dot(tr.translation_matrix(tra), tr.rotation_matrix(rot,axis) )

    def to_translation_scale_rotation(self): 
        # FIXME
        mat = self.data.copy()
        tra  = tr.translation_from_matrix(mat)
        tra_mat = tr.translation_matrix(tra)
        mat = numpy.dot(numpy.linalg.inv(tra_mat), mat)
        factor, origin, direction = tr.scale_from_matrix(mat)
        scale_mat = tr.scale_matrix(factor, origin, direction)
        scale = numpy.diag(scale_mat) 
        mat = numpy.dot(numpy.linalg.inv(scale_mat), mat) 
        try: 
            rot,rot_axis,point = tr.rotation_from_matrix(mat) 
        except ValueError: 
            #print "Not a rotation matrix. "
            return None, None, None, None
        #print "Rotation axis: ",rot_axis
        return tra, scale, rot, rot_axis
 
    def to_quaternion_translation(self): 
        tra = tr.translation_from_matrix(self.data)
        qua = tr.quaternion_from_matrix(self.data)
        return qua, tra 

    def from_quaternion_translation(self,quaternion,translation): 
        rot = quaternion_matrix(quaternion)
        return numpy.dot(translation, rot)

    def derivative_parameters(self, gradient_transformed_image ): 
        pass 
        # FIXME: implement (implemented in case of 6DOF in the derived class Transform_6DOF)
        
    def __get_inverse(self):
        inverse = Transform_Affine(data = self.__inverse, map_to = self.map_from, map_from = self.map_to)
        return inverse

    def __compute_inverse(self):
        self.__inverse = numpy.linalg.inv(self.data)

    def __get_data(self):
        return self.__data
    
    def set_data(self,data): 
        if not isinstance(data, numpy.ndarray): 
            if data  is None: 
                data = numpy.eye(4)
            else:
                # FIXME: raise exception
                return 
        if not data.size == 16: 
            # FIXME: raise exception 
            return 
        self.__data = data 
        self.__compute_inverse() 

    def __get_map_to(self):
        return self.__map_to
    
    def set_map_to(self, map_to): 
        self.__map_to = map_to
    
    def __get_map_from(self):
        return self.__map_from
    
    def set_map_from(self, map_from): 
        self.__map_from = map_from

    def __repr__(self): 
        with printoptions(precision=4, suppress=True):
            s = "Transformation "+"\n-from: "+str(self.map_from)+"\n-to:   "+str(self.map_to)+"\n-matrix: \n"+self.__data.__repr__()
        return s

    def _repr_html_(self):
        if not has_ipy_table: 
            return "Please install ipy_table."
        s = 'Transformation'
        table_data = [[s,'','','','','',''],
                      ['map_from','map_to','determinant','is_rigid','is_6DOF','is_rotation','is_translation'],
                      [self.map_from, self.map_to, str(self.determinant()), self.is_rigid(), self.is_6DOF(), self.is_rotation(), self.is_translation() ], ] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic')
        table = ipy_table.set_cell_style(0,0, column_span=7)
        table = ipy_table.set_cell_style(0,0, align='center') 
        table = ipy_table.set_row_style(1, color='#F7F7F7')
        table = ipy_table.set_row_style(2, color='#F7F7F7')
        table = ipy_table.set_row_style(0, color='#FFFF99')
        table = ipy_table.set_row_style(1, bold = True)
        table = ipy_table.set_row_style(1, align = 'center')
        table = ipy_table.set_row_style(1, height = 30)
        table = ipy_table.set_global_style(float_format="%3.3f") 
        table2 = ipy_table.make_table(self.data)
        table3 = ipy_table.make_table(['   '])
        table3 = ipy_table.set_row_style(0, no_border='all')
        s = '<div><style>table {float: left; margin-right:10px;}</style> %s %s %s </div>'%( table._repr_html_(), table3._repr_html_() ,table2._repr_html_() )
        return s

    data     = property(__get_data, set_data)
    map_to   = property(__get_map_to, set_map_to)
    map_from = property(__get_map_from, set_map_from)
    inverse  = property(__get_inverse )





class Transform_Identity(Transform_Affine): 
    def __init__(self,map_from="",map_to=""): 
        Transform_Affine.__init__(self, numpy.eye(4),map_from,map_to)     


#class Transform_Scale(Transform_Affine): 
#    def __init__(self, factor, origin=None, direction=None, map_from="",map_to="" ): 
#        mat = tr.scale_matrix(factor, origin, direction)
#        Transform_Affine.__init__(self, mat, map_from,map_to)
class Transform_Scale(Transform_Affine):  #FIXME: figure out what transformations.py means by scale matrix 
    def __init__(self, scale_xyz, map_from="",map_to="" ): 
        mat = numpy.diag([scale_xyz[0],scale_xyz[1],scale_xyz[2],1]) 
        Transform_Affine.__init__(self, mat, map_from, map_to)


class Transform_Rotation(Transform_Affine): 
    def __init__(self, angle, direction, point=None, map_from="",map_to=""): 
        mat = tr.rotation_matrix(angle, direction, point)
        Transform_Affine.__init__(self, mat,map_from,map_to)

    def derivative_parameters(self, gradient_transformed_image ): 
        pass 

class Transform_Translation(Transform_Affine): 
    def __init__(self, direction, map_from="",map_to=""): 
        mat = tr.translation_matrix(direction)     
        Transform_Affine.__init__(self, mat,map_from,map_to)

    def derivative_parameters(self, gradient_transformed_image ): 
        pass

class Transform_6DOF(Transform_Affine): 
    def __init__(self, translation, rotation_angle, rotation_direction, rotation_point=None, map_from="",map_to=""): 
        rot = tr.rotation_matrix(rotation_angle, rotation_direction, rotation_point) 
        tra = tr.translation_matrix(translation) 
        mat = numpy.dot(tra,rot)
        Transform_Affine.__init__(self, mat,map_from,map_to) 

    def derivative_parameters(self, gradient_transformed_image, grid_transformed_image ): 
        dRx = []
        dRy = []
        dRz = []
        return [0,0,0,0,0,0]
        






class GridND(object):
    def __init__(self, data=None, space="", is_uniform=True, is_affine=True, is_axis_aligned=True): 
        self.ndim = None
        self.__set_data(data)
        self.space = space 
        self.__clear_cache() 
        self.__is_uniform = is_uniform 
        self.__is_affine  = is_affine 
        self.__is_axis_aligned = is_axis_aligned 

    def min(self): 
        if self.__min  is None: 
            self.__min = eval('self.data%s'%('.min(0)'*self.ndim))
        return self.__min

    def max(self): 
        if self.__max  is None:
            self.__max = eval('self.data%s'%('.max(0)'*self.ndim))
        return self.__max 

    def span(self): 
        if self.__span  is None:
            self.__span = self.max() - self.min() 
        return self.__span 

    def center(self,use_corners_only=True): 
        if self.__center  is None or (self.__use_corners_only != use_corners_only): 
            if use_corners_only: 
                corners = self.corners()
                center = corners.mean(1)
            else: 
                center = None 
                #FIXME: implement 
            self.__center = center
            self.__use_corners_only = use_corners_only
        return self.__center

    def mean_distance_from_point(self, point, use_corners_only=True): 
        if use_corners_only: 
            corners = self.corners()
            dist = corners - numpy.tile( numpy.asarray(point).reshape(3,1),[1,corners.shape[1]] ) 
            dist =  numpy.sqrt((dist * dist)).sum(1) / corners.shape[1]
        else:
            dist = None
            #FIXME: implement 
        return dist 

    def mean_dist_from_center(self): 
        if self.__mean_dist_center  is None: 
            center = self.center() 
            self.__mean_dist_center = self.mean_distance_from_point(center) 
        return self.__mean_dist_center
        
    def get_shape(self):
        return self.data.shape

    def is_uniform(self): 
        """Returns True if the grid is uniform. """
        return self.__is_uniform
        #FIXME: change the flags when grid is transformed 
    
    def is_affine(self): 
        """Returns True if the grid is the affine transformation of a uniform grid. """
        return self.__is_affine 

    def is_axis_aligned(self): 
        """Returns True if the grid is uniform and aligned to the x,y,z axis. """
        return self.__is_axis_aligned 

    def transform(self, affine_from_grid): 
        return transform_grid(self, affine_from_grid )

    def corners(self, homogeneous_coords=False): 
        if self.__corners  is None: 
            n_corners = 2**self.ndim 
            corners = []
            for i in range(n_corners): 
                b = eval( '['+bin(i)[2:].zfill(self.ndim).replace('0','0,').replace('1','1,') +"]" )
                b = (numpy.asarray(self.data.shape[0:self.ndim])-1) * numpy.asarray(b)
                s = str(b.tolist())[1:-1]
                corner = eval('self.data[%s,:]'%s)
                corners.append(corner)
            corners = numpy.asarray(corners).transpose() 
            if homogeneous_coords: 
                corners2 = numpy.ones((4,n_corners))
                corners2[0:3,:] = corners
                corners = corners2
            self.__corners = corners
        return self.__corners

    def __get_data(self):
        return self.__data 
    
    def __set_data(self,data): 
        self.__data = data
        if data is not None: 
            self.ndim = self.data.ndim-1 
        else: 
            self.ndim = None
        self.__clear_cache() 
        
    def __get_shape(self): 
        return self.data.shape 
    
    def __clear_cache(self): 
        self.__min  = None 
        self.__max  = None 
        self.__span = None
        self.__center = None
        self.__use_corners_only = None
        self.__mean_dist_center = None
        self.__corners = None

    def __repr__(self):  
        with printoptions(precision=4):
            if self.is_axis_aligned():
                type = 'axis aligned'
            elif self.is_affine():
                type = 'affine'
            elif self.is_uniform():
                type = 'uniform'
            else: 
                type = ''
            if type!='': 
                type = " (%s)"%type
            s = "Grid%dD %s: "%(int(self.ndim),type)
            s = s + "\n - space -------------------> " + str(self.space)
            s = s + "\n - shape -------------------> " + str(list(self.shape)) 
            s = s + "\n - min ---------------------> " + str(self.min()) 
            s = s + "\n - max ---------------------> " + str(self.max()) 
            s = s + "\n - span --------------------> " + str(self.span()) 
            s = s + "\n - center ------------------> " + str(self.center()) 
            s = s + "\n - mean dist from center: --> " + str(self.mean_dist_from_center()) 
            s = s + "\n"
        return s
    
    def _repr_html_(self): 
        if not has_ipy_table: 
            return "Please install ipy_table."
        if self.is_axis_aligned():
            type = 'axis aligned'
        elif self.is_affine():
            type = 'affine'
        elif self.is_uniform():
            type = 'uniform'
        else: 
            type = ''
        if type!='': 
            type = " (%s)"%type
        s = "Grid%dD %s: "%(int(self.ndim),type)
        def pretty(list):
            return str(['{0:.2f}'.format(flt) for flt in list])
        table_data = [[s,'','','','','',''],
                      ['space','shape','min','max','span','center','spread'],
                      [str(self.space), pretty(self.shape) ,pretty(self.min()) ,pretty(self.max()) ,pretty(self.span()) ,pretty(self.center()) ,pretty(self.mean_dist_from_center())], ] 
        table = ipy_table.make_table(table_data)
        table = ipy_table.apply_theme('basic')
        table = ipy_table.set_cell_style(0,0, column_span=7)
        table = ipy_table.set_cell_style(0,0, align='center') 
        table = ipy_table.set_row_style(1, color='#F7F7F7')
        table = ipy_table.set_row_style(2, color='#F7F7F7')
        table = ipy_table.set_row_style(0,color='#CCFF99')
        table = ipy_table.set_row_style(1, bold = True)
        table = ipy_table.set_row_style(1, align = 'center')
        table = ipy_table.set_global_style(float_format="%3.3f")        
        s = table._repr_html_() 
        return s

    data = property(__get_data, __set_data )
    shape = property(__get_shape )

    
class Grid3D(GridND):
    def __init__(self, data=None, space=""): 
        GridND.__init__(self, data, space)
        self.ndim = 3






class ImageND(object): 
    def __init__(self, data=None, affine=None, space="", mask_flag=0):
        self.ndim   = None
        self.space = ""
        if isinstance(data,str): 
            self.load_from_file(data)
        else: 
            self.set_data(data)
        self.set_affine(affine)
        self.set_space(space)

        self.background = global_settings.get_default_background() 
        self.use_gpu    = global_settings.is_gpu_enabled() 
        
        self.set_mask_flag(mask_flag)

    def set_mask_flag(self, is_mask): 
        self.__is_mask = is_mask
    
    def is_mask(self): 
        """Returns True if the image is a mask (i.e. the values should be interpreted as discrete labels - take note 
        that, however, that this reports the state of a flag; the data can be non-integer. )"""
        return self.__is_mask 

    def set_lookup_table(self,lut): 
        """Set a lookup table for visualization. """
        self.__lut = lut 

    def get_lookup_table(self): 
        try: 
            lut = self.__lut
        except: 
            lut = None 
        return lut

    def set_data(self,data): 
        if isinstance(data,numpy.ndarray): 
            self.data = data
        else:   
            self.data = numpy.asarray([]) 
        self.ndim = self.data.ndim

    def min(self): 
        return self.data.min() #FIXME: pass on optional parameters
    
    def max(self): 
        return self.data.max() #FIXME: pass on optional parameters

    def set_affine(self,affine): 
        if isinstance(affine,Transform_Affine): 
            self.affine = affine
        elif isinstance(affine,numpy.ndarray): 
            self.affine = Transform_Affine(data=affine)
        elif affine is None: 
            self.affine = Transform_Identity()  
        else: 
            print "'affine' must be an instance of Affine"
        self.affine.map_to = self.space        
        
    def set_space(self,space): 
        self.affine.map_to = space
        self.space = space 
        self.affine.map_from = "index" 

    def get_shape(self): 
        return self.data.shape 

    def get_world_grid(self, n_points=None): 
        if n_points  is None: 
            n_points = self.get_shape() 
        grid = grid_from_box_and_affine(self.get_world_grid_min(), self.get_world_grid_max(), n_points) 
        return grid

    def get_world_grid_min(self): 
        d = self.get_data()
        s = numpy.asarray(d.shape)-1
        corners = numpy.asarray([ [0,0,0,1], [s[0],0,0,1], [0,s[1],0,1], [s[0],s[1],0,1], [0,0,s[2],1], [s[0],0,s[2],1], [0,s[1],s[2],1], [s[0],s[1],s[2],1] ]).transpose()
        corners = self.affine.left_multiply(corners)
        corners=corners[0:3,:]
        m = corners.min(1)
        return m
         
    def get_world_grid_max(self): 
        d = self.get_data()
        s = numpy.asarray(d.shape)-1
        corners = numpy.asarray([ [0,0,0,1], [s[0],0,0,1], [0,s[1],0,1], [s[0],s[1],0,1], [0,0,s[2],1], [s[0],0,s[2],1], [0,s[1],s[2],1], [s[0],s[1],s[2],1] ]).transpose()
        corners = self.affine.left_multiply(corners)
        corners=corners[0:3,:]
        M = corners.max(1)
        return M

    def get_pixel_grid(self): 
        print "get_pixel_grid"
        n_points = numpy.uint32(self.shape)
        grid = grid_from_box_and_affine(numpy.asarray([0,0,0]), n_points-1, n_points) 
        return grid

    def transform(self,affine): 
        if not isinstance(affine,Transform_Affine): 
            print "Transformation must be an instance of Affine."
            # FIXME raise error
            return None
        self.affine = affine.left_multiply(self.affine)
        #return self 

    def save_to_file(self, filename): 
        nii = occiput_to_nifti(self)
        nibabel.save(nii,filename)

    def __get_shape(self):
        return self.data.shape 
        
    def __get_size(self): 
        return self.data.size 

    def copy(self): 
        return copy.copy(self)

    shape = property(__get_shape)
    size  = property(__get_size)


class Image3D(ImageND):
    def __init__(self, data=None, affine=None, space=""): 
        ImageND.__init__(self, data, affine, space)
        self.ndim = 3 

    def compute_resample_on_grid(self,grid,affine_grid_to_world=None, verify_mapping=True,interpolation_mode=INTERPOLATION_LINEAR): 
        resampled_data = resample_image_on_grid(self, grid, affine_grid_to_world, verify_mapping, self.background, self.use_gpu, interpolation_mode)
        # create new Image3D object
        return Image3D(data=resampled_data) #FIXME: perhaps just return the raw resampled data
    
    def compute_resample_in_box(self,box_min,box_max,box_n): 
        pass 

    def compute_resample_in_space_of_image(self,image): 
        pass 

    def compute_gradient_on_grid(self, grid, affine_grid_to_world=None, verify_mapping=True): 
        resampled_data = resample_image_on_grid(self, grid, affine_grid_to_world, verify_mapping, self.background, self.use_gpu)
        # create new Image3D object
        #print resampled_data.max()
        gradient_data = numpy.gradient(resampled_data) #FIXME: use NiftyPy
        return gradient_data
        
    def compute_gradient_in_box(self,box): 
        pass 
    
    def compute_gradient_in_space_of_image(self,image): 
        pass 

    def compute_gradient_pixel_space(self):
        gradient_data = numpy.gradient(self.data)  #FIXME: use NiftyPy
        return gradient_data

    def compute_smoothed(self, smoothing): 
        pass 

    def get_data(self):
        return self.data 

    def export_image(self,index,axis=0,normalise=True,scale_factor=None,shrink=None,rotate=0): 
        pass 

    def has_data(self): 
        return self.data is not None

    def min(self): 
        return self.data.min()
    
    def max(self): 
        return self.data.max()

    def display_in_browser(self, axis=0): 
        D = self.display(axis) 
        D.display_in_browser() 

    def display(self, axis=0,shrink=None,rotate=None,subsample_slices=None,scales=None,open_browser=None): 
        # The following is a quick fix: use MultipleVolumesNiftyPy if the image has small size, 
        # MultipleVolumes otherwise. MultipleVolumesNiftyPy makes use of the GPU but crashes with large images. 
        # NOTE that MultipleVolumes produces a different visualisation from MultipleVolumesNiftyPy: 
        # it displays the raw imaging data, without accounting for the transformation to world space. 
        # Modify MultipleVolumesNiftyPy so that it processes the data in sequence if it is too large for the GPU. 
        # Then get rid of  MultipleVolumes. 
        #if self.size <= 256**3: 
        #    D = MultipleVolumesNiftyPy([self],axis,open_browser=open_browser) 
        #else: 
        D = MultipleVolumes([self],axis=axis,shrink=shrink,rotate=rotate,subsample_slices=subsample_slices,scales=scales,open_browser=open_browser)
        return D 

    def display_with(self, other_images, axis=0,shrink=None,rotate=None,subsample_slices=None,scales=None,open_browser=None): 
        if isinstance(other_images, type( () )): 
            other_images = list(other_images)
        elif isinstance(other_images, type( [] )): 
            other_images = other_images
        else: 
            other_images = [other_images,]
        D = MultipleVolumes([self,]+other_images,axis=axis,shrink=shrink,rotate=rotate,subsample_slices=subsample_slices,scales=scales,open_browser=open_browser)
        return D 

    def display_slice(self,axis=0,index=None,open_browser=False): 
        if index is None: 
            if axis==0: 
                index = numpy.uint32(self.shape[0]/2)
            if axis==1: 
                index = numpy.uint32(self.shape[1]/2)
            if axis==2: 
                index = numpy.uint32(self.shape[2]/2)
            if axis==0: 
                a = self.data[index,:,:].reshape((self.shape[1],self.shape[2]))  
            if axis==1: 
                a = self.data[:,index,:].reshape((self.shape[0],self.shape[2]))  
            if axis==2: 
                a = self.data[:,:,index].reshape((self.shape[0],self.shape[1]))          
        D = DisplayNode.DisplayNode()
        im = Image.fromarray(a).convert("RGB").rotate(90)
        D.display('image', im, open_browser) 
        return D 
    
    def volume_render(self,scale=1.0): 
        V = VolumeRenderer(self,scale=scale) 
        return V.display() 

    def _repr_html_(self): 
        display = self.display()
        if display: 
            return display._repr_html_()


    # Overload math operators
    def _is_same_type(self, other): 
        return isinstance(other,self.__class__)

    def _is_in_same_space(self, other): 
        # FIXME: implement
        return True 

    def _is_on_same_grid(self, other): 
        # FIXME: implement
        return True 

    def __add__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = out.data+other.data 
                    return out
                else: 
                    #FIXME: implement 
                    print "SUM of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "SUM of images not in the same space. It cannot be done. "
        else: 
            out=self.copy()
            out.data = out.data+other
            return out
        return None 

    def __sub__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = out.data-other.data 
                    return out
                else: 
                    #FIXME: implement 
                    print "SUB of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "SUB of images not in the same space. It cannot be done. "
        else: 
            out=self.copy()
            out.data = out.data-other
            return out
        return None
    
    def __mul__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = out.data*other.data 
                    return out
                else: 
                    #FIXME: implement 
                    print "MUL of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "MUL of images not in the same space. It cannot be done. "
        else: 
            out=self.copy()
            out.data = out.data*other
            return out
        return None
    
    def __div__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = out.data/other.data 
                    return out
                else: 
                    #FIXME: implement 
                    print "DIV of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "DIV of images not in the same space. It cannot be done. "
        else: 
            out=self.copy()
            out.data = out.data/other
            return out
        return None
            
    def __radd_(self, other): 
        return self.__add__(other) 

    def __rmul_(self, other): 
        return self.__mul__(other)

    def __rsub__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = other.data - out.data
                    return out
                else: 
                    #FIXME: implement 
                    print "SUB of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "SUB of images not in the same space. It cannot be done. "
        else: 
            out=self.copy()
            out.data = other-out.data
            return out
        return None

    def __rdiv__(self, other): 
        if self._is_same_type(other): 
            if self._is_in_same_space(other): 
                if self._is_on_same_grid(other): 
                    out = self.copy() 
                    out.data = other.data / out.data
                    return out
                else: 
                    #FIXME: implement 
                    print "DIV of images on different grids; the right hand side image must be resampled, please implement this."
            else: 
                # FIXME: raise error 
                print "DIV of images not in the same space. It cannot be done. "
        else: 
            out = self.copy()
            out.data = other / out.data
            return out
        return None
    




