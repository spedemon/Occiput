
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2014, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA
# April 2014
# Nov. 2015 

# Import an interfile projection as a PET_Projection object of occiput and export. 


__all__ = ["import_interfile_projection", "export_interfile_projection", "import_PET_Projection"] 


from occiput.Reconstruction.PET.PET_projection import PET_Projection, Binning, PET_Projection_Sparsity
from interfile import Interfile
import h5py
import os
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose  


def import_interfile_projection_data(headerfile='', datafile='',load_time=False):   #FIXME: this should be in the Interfile package
        F = Interfile.load(headerfile)
        
        if F.has_key('matrix size[1]'): 
            N_planes    = F['matrix size[1]']['value']
            N_axial     = F['matrix size[2]']['value']
            N_sinograms = F['matrix size[3]']['value']        
        else:
            N_planes    = F['matrix size [1]']['value']
            N_axial     = F['matrix size [2]']['value']
            N_sinograms = F['matrix size [3]']['value']    
        if datafile == '': 
            datafile1 = headerfile.replace(headerfile.split(os.sep)[-1],F['name of data file']['value'])
            datafile2 = headerfile.replace('.s.hdr','.s') 
            datafile2 = datafile2.replace('.h33','.a')
            datafile3 = headerfile.replace('.h33','.s')
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
        data = data.reshape([N_sinograms,N_axial,N_planes])
        if load_time: 
            try: 
                duration = int32([0,F['image duration']['value']])*1000
            except: 
                print "Unable to load image (sinogram) duration. This may determine an incorrect scale and use of randoms and scatter when reconstructing. Set .time_bins manually. "
                duration = int32([0,0])
        else: 
            duration = int32([0,0])
        return data, duration

def import_interfile_projection(headerfile, binning, michelogram, datafile='' ,invert=False, vmin=0.00, vmax=1e10,load_time=False): 
        data, duration = import_interfile_projection_data(headerfile, datafile, load_time=load_time) 
        N_planes    = data.shape[2]
        N_axial     = data.shape[1]
        N_sinograms = data.shape[0]

        N_azim      = michelogram.n_segments
        z_elements = int32([michelogram.segments_sizes[0]])
        for i in range((N_azim-1)/2): 
            z_elements = concatenate([z_elements,int32([michelogram.segments_sizes[i+1],michelogram.segments_sizes[i+1]])])
        
        projection = zeros([binning.N_axial, binning.N_azimuthal, binning.N_u,  binning.N_v], dtype=float32, order="C")
        indexes_azim = [5,6,4,7,3,8,2,9,1,10,0] 
        index0 = 0

        for i in range(N_azim): 
            index1 = index0 + z_elements[i]
            data_azim = data[index0:index1,:,:] 
            for j in range(N_axial): 
                plane = zeros([projection.shape[2],projection.shape[3]],order='F')
                offset = (projection.shape[3]-(index1-index0))/2
                plane[:,offset:offset+index1-index0] = data_azim[:,j,:].squeeze().transpose()
                projection[j,indexes_azim[i],:,:] = fliplr(plane)
            index0 += z_elements[i]
            
        # flip azimuthally - this makes it consistent with Occiput's routines that import from Petlink32 listmode 
        projection2 = zeros(projection.shape,dtype=float32,order="C")
        projection2[:,5,:,:] = projection[:,5,:,:] 
        projection2[:,4,:,:] = projection[:,6,:,:] 
        projection2[:,3,:,:] = projection[:,7,:,:] 
        projection2[:,2,:,:] = projection[:,8,:,:] 
        projection2[:,1,:,:] = projection[:,9,:,:] 
        projection2[:,0,:,:] = projection[:,10,:,:] 
        projection2[:,6,:,:] = projection[:,4,:,:] 
        projection2[:,7,:,:] = projection[:,3,:,:] 
        projection2[:,8,:,:] = projection[:,2,:,:] 
        projection2[:,9,:,:] = projection[:,1,:,:] 
        projection2[:,10,:,:] = projection[:,0,:,:]
        
        sparsity = PET_Projection_Sparsity(binning.N_axial, binning.N_azimuthal, binning.N_u,  binning.N_v)
        ## invert, except where values are 0 
        # if there are broken detectors, disable them (set sensitivity to 0) 
        if invert: 
            projection_fixed = projection2.copy() 
            projection_fixed[projection_fixed<vmin]=0.0 
            projection_fixed[projection_fixed>vmax]=vmax 

            projection_inv = zeros(projection2.shape)
            projection_inv[projection_fixed!=0] = 1.0/projection_fixed[projection_fixed!=0]

            P = PET_Projection( binning, projection_inv, sparsity.offsets, sparsity.locations, duration) 
        else:
            P = PET_Projection( binning, projection2, sparsity.offsets, sparsity.locations, duration) 
        return P 


def export_interfile_projection(sinogram_data_file, projection_data, binning, michelogram, invert=False): 
    # projection_data has size  N_axial, N_azimuthal, N_u, N_v      dtype=float32     order="C"
    # export as  N_sinograms, N_axial, N_planes 
    
    # FIXME: need to flip azimuthally as in import_interfile 

    N_planes = binning.N_u
    N_axial  = binning.N_axial
    N_sinograms = int32(michelogram.segments_sizes[0]) + 2*int32(michelogram.segments_sizes[1::]).sum()
    data = zeros([N_sinograms, N_axial, N_planes], dtype=float32, order='C') 
    
    N_azim      = michelogram.n_segments
    z_elements = int32([michelogram.segments_sizes[0]])
    for i in range((N_azim-1)/2): 
        z_elements = concatenate([z_elements,int32([michelogram.segments_sizes[i+1],michelogram.segments_sizes[i+1]])])
    indexes_azim = [5,6,4,7,3,8,2,9,1,10,0]
    index0=0 

    for i in range(N_azim):
        index1 = index0 + z_elements[i]
        data_azim = zeros( [index1-index0,data.shape[1],data.shape[2]], dtype=float32, order='C' )
        for j in range(N_axial): 
            plane = fliplr( projection_data[j,indexes_azim[i],:,:] )
            offset = (projection_data.shape[3]-(index1-index0))/2
            data_azim[:,j,:] = plane[:,offset:offset+index1-index0].transpose() 
        data[index0:index1,:,:] = data_azim
        index0 += z_elements[i]

    if invert: 
        data2 = zeros([N_sinograms, N_axial, N_planes], dtype=float32, order='C') 
        data2[data!=0] = 1.0 / data[data!=0]
    else: 
        data2 = data
    data2.tofile(sinogram_data_file)



def import_PET_Projection( filename ): 
    h5f = h5py.File(filename,'r') 
    offsets   = asarray( h5f['offsets'],order='F' ) 
    locations = asarray( h5f['locations'],order='F'  )
    try: 
        data      = asarray( h5f['data'],order='F'  )
    except: 
        data      = asarray( h5f['static_data'],order='F'  ) # compatibility with first version - deprecate at some point
    try: 
        time_bins = asarray( h5f['time_bins'],order='F'  )
    except: 
        time_bins = int32([0,0])                   # compatibility with first version - deprecate at some point
    binning = Binning(    {'n_axial':                int32( h5f['n_axial'] ), 
                           'n_azimuthal':            int32( h5f['n_azimuthal'] ), 
                           'angles_axial':           float32( h5f['angles_axial'] ), 
                           'angles_azimuthal':       float32( h5f['angles_azimuthal'] ), 
                           'size_u':                 float32( h5f['size_u'] ), 
                           'size_v':                 float32( h5f['size_v'] ),
                           'n_u':                    int32( h5f['n_u'] ),
                           'n_v':                    int32( h5f['n_v'] ) })
    h5f.close() 
    return PET_Projection(binning, data, offsets, locations, time_bins) 

