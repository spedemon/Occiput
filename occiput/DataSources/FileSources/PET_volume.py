
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2014, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA
# April 2014
# Nov. 2015 


# Import an interfile volume as an Image3D and export. 


from occiput.Core import Image3D, Transform_Identity, Transform_6DOF, Transform_Affine, Transform_Scale
from interfile import Interfile
import h5py
import os
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose  


__all__ = ["import_interfile_volume","export_interfile_volume"]

def import_interfile_volume_data(headerfile='', datafile=''):  #FIXME: this should be in the Interfile package
        F = Interfile.load(headerfile)
        if F.has_key('matrix size[1]'): 
            Nx    = F['matrix size[1]']['value']
            Ny    = F['matrix size[2]']['value']
            Nz    = F['matrix size[3]']['value']        
        else:
            Nx    = F['matrix size [1]']['value']
            Ny    = F['matrix size [2]']['value']
            Nz    = F['matrix size [3]']['value']    
        if datafile == '': 
            datafile1 = headerfile.replace(headerfile.split(os.sep)[-1],F['name of data file']['value'])
            datafile2 = headerfile.replace('.v.hdr','.v') 
            datafile2 = datafile2.replace('.h33','.v')
            datafile3 = headerfile.replace('.h33','.v')
            try: 
                data = fromfile(datafile1,dtype=float32)
            except: 
                try: 
                    data = fromfile(datafile2,dtype=float32)
                except: 
                    try: 
                        data = fromfile(datafile3,dtype=float32)
                    except: 
                        print "Data file not found."
        else: 
            data = fromfile(datafile,dtype=float32)
        data = data.reshape([Nz,Ny,Nx])
        data = asfortranarray(data.swapaxes(0,2))
        data = transpose(data,[1,0,2])
        data = data[::-1,:,:]
        return data    

def import_interfile_volume(headerfile='', datafile=''): 
        # Load ndarray data 
        data = import_interfile_volume_data(headerfile, datafile) 
        # Load other information - e.g. pixels size 
        F = Interfile.load(headerfile)
        if F.has_key('scale factor (mm/pixel) [1]'): 
            pixsize_x = F['scale factor (mm/pixel) [1]']['value']
            pixsize_y = F['scale factor (mm/pixel) [2]']['value']
            pixsize_z = F['scale factor (mm/pixel) [3]']['value']
        # Create Image3D 
        T_pix_to_world = Transform_Scale(int32([pixsize_x,pixsize_y,pixsize_z]), map_from='pixels_PET', map_to='world') 
        volume = Image3D(data=data, affine=T_pix_to_world, space='world')  
        return volume


def export_interfile_volume(data_file_name, data): 
        data = data[::-1,:,:]
        data = transpose(data,[1,0,2])
        data = data.swapaxes(0,2)
        data = asarray(data, dtype=float32, order='C')
        data.tofile(data_file_name)




