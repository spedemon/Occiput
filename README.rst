============================
Occiput - Tomographic Vision
============================

Tomographic reconstruction software for PET, PET-MRI and SPECT in 2D, 3D (volumetric) and 4D (spatio-temporal) in Python. 

The software provides high-speed reconstruction using Graphics Processing Units (GPU). 
Note: an NVidia CUDA-compatible GPU is required.  

``Occiput`` can be utilized with arbitrary scanner geometries. It can be utilized for abstract tomographic 
reconstruction experiments to develop new algorithms and explore new system geometries, or to connect to real-world scanners, 
providing production quality image reconstruction with standard (MLEM, OSEM, Ordinary Poisson OSEM) and advanced algorithms. 

``Occiput`` implements algorithms for motion correction (direct motion estimation), kinetic imaging, multi-modal reconstruction, respiratory and cardiac gated imaging. 
The source code contains Jupyter notebooks with examples. 

Occiput has a plugin module designed to facilitate the creation of interfaces for real-world PET and SPECT scanners. 
A Python extension package ``Occiput_Interface_Biograph_mMR``, implementing the interface plugin for the Siemens Biograph mMR PET-MRI scanner 
is available upon request and following authorization from Siemens. Notebooks containing Biograph_mMR in the title can 
only be executed after installing the extension package. 
Please email us at occiput.reconstruction@gmail.com 


Installation 
============

Linux, Windows (not tested recently), MacOS
-------------------------------------------

Pre-requisites: Occiput requires ``NVidia GPU Drivers``, ``NVidia CUDA`` and the ``NiftyRec`` GPU accelerated tomographic ray-tracing library. 

1. `Install NVidia GPU Drivers and CUDA <https://developer.nvidia.com/cuda-downloads>`_

2. `Install NiftyRec libraries <http://niftyrec.scienceontheweb.net>`_ - build the latest version using CMake
    
3. Make sure that CUDA libraries and NiftyRec libraries are in the system path: 

 - Linux: 
 
    export LD_LIBRARY_PATH:$LD_LIBRARY_PATH:/path_to_cuda_libraries:/path_to_niftyrec_libraries
    
 - MacOS: 

    export DYLD_LIBRARY_PATH:$DYLD_LIBRARY_PATH:/path_to_cuda_libraries:/path_to_niftyrec_libraries

 - Windows: 

    setx path "%path%;c:/path_to_cuda_libraries:/path_to_niftyrec_libraries;"

4. Install ``Occiput``: 

    git clone https://github.com/spedemon/occiput.git 

    python setup.py build install 


Getting started
===============
Examples and demos of the features of Occiput are in the /occiput/notebooks folder. 
To get started, install ``Python Jupyter`` and open the scripts in 
`/occiput/notebooks <https://github.com/spedemon/occiput/tree/master/occiput/notebooks>`_. The 
notebook `/occiput/notebooks/DOCUMENTATION.ipynb <https://github.com/spedemon/occiput/tree/master/occiput/notebooks/DOCUMENTATION.ipynb>`_ contains 
an index and short description of the notebooks. 

Website
=======
For more information see `occiput.io  <http://tomographylab.scienceontheweb.net/>`_. 




