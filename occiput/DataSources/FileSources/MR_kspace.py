
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2014, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA
# March 2015


from occiput.DataSources.FileSources.Volume import import_nifti

def import_kspace(filename): 
    return import_nifti(filename) 

