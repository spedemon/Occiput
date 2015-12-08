# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# Dec. 2014, Boston, MA 
# April. 2014, Boston, MA


import occiput as __occiput
import nibabel as __nibabel
import warnings as __warnings
with __warnings.catch_warnings():
    __warnings.simplefilter("ignore")
    import nipy as __nipy
    from nipy.io.nifti_ref import nifti2nipy as nifti_to_nipy



def nipy_to_occiput(nib): 
    if nib.ndim == 3: 
        im = __occiput.Core.Image3D(data=nib.get_data(), affine=nib.affine, space="world") 
    else: 
        im = __occiput.Core.ImageND(data=nib.get_data(), affine=nib.affine, space="world")   
    return im 



def nifti_to_occiput(nif): 
    nip = nifti_to_nipy(nif)
    return nipy_to_occiput(nip)



def occiput_to_nifti(occ): 
    nii = __nibabel.nifti1.Nifti1Image(occ.data,occ.affine.data) 
    return nii
    


def occiput_from_array(array):
    if array.ndim == 3: 
        im = __occiput.Core.Image3D(data=array, space="world") 
    else: 
        raise("Currently only conversion of 3D arrays is supported. ")
    return im 