#!/usr/bin/python
# -*- coding: UTF-8 -*-
          
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA


from occiput.Visualization.Visualization import ProgressBar
from occiput.DataSources.FileSources.LookupTable import load_freesurfer_lut_file

import subprocess
import logging
import sys
import inspect
import os


class Downloader_HTTP(): 
    def __init__(self): 
        self._filename = 'unknown' 
        self._output = None 
        self._progress_bar = ProgressBar() 
        self._progress_bar.set_percentage(0) 
        self._verbose = False 

    def set_verbose(self, verbose): 
        self._verbose = verbose 

    def _set_percentage(self,percentage): 
        self._progress_bar.set_percentage(percentage) 

    def download(self, url, output=None, overwrite=False): 
        # Check if output exists: 
        if not overwrite: 
            if output!=None: 
                if os.path.exists(output): 
                    if self._verbose: 
                        print "File",output,"exists, not dowloading." 
                    self._set_percentage(100) 
                    return output 
        self._set_percentage(0) 
        self._output = output 
        if output==None: 
            args = ['wget',url,'--continue']
        else: 
            args = ['wget',url,'--continue','-O',output]
        try:
            pipe = subprocess.Popen(args, bufsize = 0,
                shell = False,
                stdout = None, # no redirection, child use parent's stdout
                stderr = subprocess.PIPE) # redirection stderr, create a new pipe, from which later we will read
        except Exception as e:
            #inspect.stack()[1][3] will get caller function name
            logging.error(inspect.stack()[1][3] + ' error: ' + str(e))
            return False
        while 1:
            s = pipe.stderr.readline()
            if s:
                p = self._strip_percentage(s)
                if p!=None: 
                    self._set_percentage(p) 
                if self._verbose:  
                    print s
                name = self._strip_filename(s) 
                if name!=None: 
                    self._set_filename(name)
                #sys.stdout.write(s)
            if pipe.returncode is None:
                code = pipe.poll()
            else:
                break
        if not 0 == pipe.returncode:
            self._set_percentage(0) 
            return False
        self._set_percentage(100)
        if output: 
            return output
        else: 
            return self._filename

    def _strip_percentage(self,s):
        p = s[s.find('%')-2:s.find('%')].strip(' ') 
        try: 
            percentage = int(p)
        except: 
            return None
        else: 
            return percentage
 
    def _strip_filename(self,s): 
        if s.find("Saving to")!=-1: 
            name = s.strip("Saving to").strip("\n").strip("“").strip("'").strip("`")
            if self._verbose: 
                print "Detected name: ",name
            return name 

    def _set_filename(self,name): 
        self._filename = name 

class Dropbox(Downloader_HTTP): 
    pass 
    

def download_Dropbox(url, output=None, overwrite=False, verbose=False): 
    D = Dropbox()
    D.set_verbose(verbose)
    return D.download(url, output, overwrite) 

