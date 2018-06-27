# -*- coding: utf-8 -*-
# occiput  
# Harvard University, Martinos Center for Biomedical Imaging 
# Aalto University, Department of Computer Science



class FileNotFound(Exception): 
    def __init__(self,msg,filename): 
        self.msg = str(msg) 
        self.filename = str(filename)
    def __str__(self): 
        return "Cannot find file '%s' (%s)."%(self.filename, self.msg)

class UnknownParameter(Exception): 
    def __init__(self,msg): 
        self.msg = str(msg) 
    def __str__(self): 
        return "Unkwnown parameter: %s"%(self.msg)

class UnexpectedParameter(Exception): 
    def __init__(self,msg): 
        self.msg = str(msg) 
    def __str__(self): 
        return "Unexpected parameter: %s"%(self.msg)



