# -*- coding: utf-8 -*-
# occiput  
# Harvard University, Martinos Center for Biomedical Imaging 
# Aalto University, Department of Computer Science

"""Conversions to and from occiput ImageND n-dimensional image. 
Conversion to and from Nifti, Nipy, ndarray. 
"""


import occiput as __occiput
import nibabel as __nibabel
import warnings as __warnings
with __warnings.catch_warnings():
    __warnings.simplefilter("ignore")
    import nipy as __nipy
    from nipy.io.nifti_ref import nifti2nipy as nifti_to_nipy



def nipy_to_occiput(img): 
    """Convert nipy image to Occiput ImageND image.
    
    Args: 
        img (nipy): nipy image. 
    
    Returns: 
        ImageND: occiput image. """
    if img.ndim == 3: 
        img2 = __occiput.Core.Image3D(data=img.get_data(),affine=img.affine,space="world") 
    else: 
        img2 = __occiput.Core.ImageND(data=img.get_data(),affine=img.affine,space="world")   
    return img2 


def nifti_to_occiput(nif): 
    """Convert Nifti image to occiput ImageND image. 
    
    Args: 
        nif (Nifti): Nifti image. 
        
    Returns: 
        ImageND: occiput image. 
    """
    nip = nifti_to_nipy(nif)
    return nipy_to_occiput(nip)


def occiput_to_nifti(occ): 
    """Conver occiput ImageND to Nifti image. 
    
    Args: 
        occ (ImageND): occiput ImageND image. 
        
    Returns: 
        Nifti: Nifti image. 
    """
    nii = __nibabel.nifti1.Nifti1Image(occ.data,occ.affine.data) 
    return nii
    

def occiput_from_array(array):
    """Numpy ndarray to occiput ImageND image. 
    
    Args: 
        array (ndarray): numpy.ndarray. 

    Returns: 
        ImageND: occiput ImageND image.  
    """
    if array.ndim == 3: 
        im = __occiput.Core.Image3D(data=array, space="world") 
    else: 
        raise("Currently only conversion of 3D arrays is supported. ")
    return im 
