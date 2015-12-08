
# ## Basic SPECT reconstruction example

import time
import sys
from occiput.Reconstruction.SPECT import SPECT_Static_Scan 
from occiput.DataSources.Synthetic.Shapes import uniform_spheres_ring
from pylab import *
ion()

def demo_SPECT_reconstruction(datatype="phantom"): 

    ### 1) Create an instance of a SPECT static scan
    spect = SPECT_Static_Scan()
    spect.set_use_gpu(True)
    
    ### 2) Set the geometry of the scanner 
    spect.set_n_pixels(128,128)
    spect.set_gantry_angular_positions(0.0, 360.0, 59) 
    
    ### 3) Load the projection data and the attenuation map  
    if datatype == "phantom": 
        phantom_activity    = uniform_spheres_ring([128,128,128],inner_value=1000.0) 
        phantom_attenuation = uniform_spheres_ring([128,128,128],inner_value=0.2) 
        measurement = spect.project(phantom_activity.data, attenuation=phantom_attenuation.data)
        spect.set_measurement(measurement.data) 
        spect.set_attenuation(phantom_attenuation.data)
    else: 
        spect.load_measurement_from_file('../data/spect/projection.nii')
        spect.load_attenuation_from_file('../data/spect/attenuation.nii')

    ### Let's give a look at the input data 
    figure(figsize=[18,6])
    subplot(1,2,1); imshow(spect._measurement[:,:,10],cmap='gray');
    subplot(1,2,2); imshow(spect._attenuation[64,:,:],cmap='gray');
    raw_input("Press Enter to start SPECT reconstruction...") 

    ### 4) Define the geometry of the reconstruction volume 
    spect.set_pixel_size(4.8,4.8)
    spect.set_radius(200.0)
    spect.set_psf(fwhm0_mm=5.0, depth_dependence=0.0001)



    ### 5) Reconstruction !
    
    activity = spect.estimate_activity(iterations=10, subset_size=16, subset_mode='random', method='EM') 

    figure(figsize=[18,6])
    subplot(131); imshow(activity.data[:,:,64],cmap='gray',vmax=180.0); 
    subplot(132); imshow(activity.data[:,80,:],cmap='gray',vmax=180.0); 
    subplot(133); imshow(activity.data[64,:,:],cmap='gray',vmax=180.0)

    ### 6) Save the tomogram 
    #activity.save_to_file('../data/spect/spect_01_reconstruction.nii')

    raw_input("Done! Press Enter to quit...") 

if __name__ == "__main__": 
    demo_SPECT_reconstruction("phantom")
