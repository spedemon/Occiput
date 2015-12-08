
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2014, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA
# April 2014
# Nov. 2015 


__all__ = ["guess_file_type_by_name"] 


def guess_file_type_by_name(filename): 
    if filename.endswith("h5"): 
        return "h5"
    elif filename.endswith(".l.hdr"): 
        return "interfile_listmode_header"
    if filename.endswith(".l"): 
        return "interfile_listmode_data"
    if filename.endswith(".h33") or filename.endswith(".s.hdr"): 
        return "interfile_projection_header"
    elif filename.endswith(".a") or filename.endswith(".s"):
        return "interfile_projection_data"
    elif filename.endswith(".v.hdr"): 
        return "interfile_volume_header"
    elif filename.endswith(".v"): 
        return "interfile_volume_data"
    elif filename.endswith(".nii.gz"): 
        return "nifti_compressed"
    elif filename.endswith(".nii"): 
        return "nifti"
    elif filename.endswith(".mat"): 
        return "mat"
    else: 
        return "unknown"
        




