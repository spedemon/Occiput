# -*- coding: utf-8 -*-
# occiput  
# Harvard University, Martinos Center for Biomedical Imaging 
# Aalto University, Department of Computer Science


def import_listmode(filename):
    return None


def convert_listmode_dicom_to_interfile(dicom_file, interfile_file_name_base): 
    print "Converting listmode dicom file to listmode interfile. "

    import dicom 
    d = dicom.read_file(dicom_file)
    
    header = d.get(('0029', '1010'))
    data   = d.get(('7fe1', '1010'))
    
    f = open(interfile_file_name_base+".l.hdr","w")
    f.write(header.value)
    f.close() 
    
    f = open(interfile_file_name_base+".l","wb")
    f.write(data.value)
    f.close() 
    print "Done. "