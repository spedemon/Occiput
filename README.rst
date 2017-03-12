============================
Occiput - Tomographic Vision
============================

Tomographic reconstruction software for PET, PET-MRI and SPECT in 2D, 3D (volumetric) and 4D (spatio-temporal) in Python. 

The software provides high-speed reconstruction using Graphics Processing Units (GPU). 
Note: an NVidia CUDA compatible GPU is required.  

``Occiput`` can be utilized with arbitrary scanner geometries. It provides production quality image reconstruction 
with standard algorithms (such as MLEM and OSEM) and implements advanced algorithms for motion correction, 
kinetic imaging and for multi-modal reconstruction. 

The source code contains Jupyter notebooks with examples. 

A Python package implementing the interface to the ``Siemens Biograph mMR`` PET-MRI scanner 
is available upon request and following authorization from Siemens. Please email us at occiput.reconstruction@gmail.com 


Installation 
============

Linux, Windows (not tested recently), MacOS
-------------------------------------------

Pre-requisites: ``Occiput`` requires ``NVidia GPU Drivers`, ``NVidia CUDA`` and the ``NiftyRec`` GPU accelerated tomographic ray-tracing library. 

1. `Install NVidia GPU Drivers and CUDA <https://developer.nvidia.com/cuda-downloads>`_

2. `Install NiftyRec libraries <http://niftyrec.scienceontheweb.net>`_ 
    
3. Make sure that CUDA libraries and NiftyRec libraries are in the system path: 

 - Linux: 
 
    export LD_LIBRARY_PATH:$LD_LIBRARY_PATH:/path_to_cuda_libraries:/path_to_niftyrec_libraries
    
 - MacOS: 

    export DYLD_LIBRARY_PATH:$DYLD_LIBRARY_PATH:/path_to_cuda_libraries:/path_to_niftyrec_libraries

 - Windows: 

    setx path "%path%;c:/path_to_cuda_libraries:/path_to_niftyrec_libraries;"

4. Install ``Occiput``: 

    python setup.py build install 


Getting started
===============
Examples and demos of the features of Occiput are in the /occiput/notebooks folder. 
To get started, install ``Python Jupyter`` and open the scripts in 
`/occiput/notebooks <https://github.com/spedemon/occiput/tree/master/occiput/notebooks>`_


Website
=======
For more information see `occiput.io  <http://www.occiput.io/>`_. 




