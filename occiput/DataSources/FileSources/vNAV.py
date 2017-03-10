
# occiput 
# Stefano Pedemonte 
# April 2014 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA 
# Nov. 2015 

__all__ = ['vNAV_MPRage','load_vnav_mprage'] 

from PIL import Image 
from occiput.Core import Image3D, Transform_Affine 
import numpy 
import nibabel 
import dicom 
import matplotlib.pyplot as plt
import os 
import copy 
from occiput.Core import transformations as tr


BOX_MIN      = [-50.0,-50.0,-50.0]
BOX_MAX      = [50.0,50.0,50.0]
THRESHOLD_MM = 1.1
LINE_COLOR   = "#2f8dff"

def quaternion_to_rotation(q,axes='sxyz'): 
    return tr.euler_from_quaternion(q,axes) 

def angle_axis_to_quaternion(angle_rad,axis):
    f = numpy.sin(0.5*angle_rad)
    quaternion = numpy.asarray( [numpy.cos(0.5*angle_rad),f*axis[0],f*axis[1],f*axis[2]] )
    return quaternion

def angle_axis_to_rotation(angle_rad,axis): 
    quaternion = angle_axis_to_quaternion(angle_rad,axis)
    rotation = quaternion_to_rotation(quaternion) 
    return rotation

def affine_from_quaternion(q,axes='sxyz'): 
    pass 

def quaternion_from_matrix(affine):
    return tx.quaternion_from_matrix(affine)

def rad_to_deg(rad): 
    return numpy.asarray(rad)*180.0/numpy.pi
    
def deg_to_rad(deg): 
    return numpy.asarray(deg)/180.0*numpy.pi




class vNAV_MPRage(): 
    def __init__(self,path=None, from_dicom_comments=True, files_start_with=None, files_end_with='.dcm'): 
        self._n_time_points     = 0 
        self._duration          = [] 
        self._motion            = [] 
        if path  is not None: 
            self.load_data_files(path, from_dicom_comments, files_start_with, files_end_with) 

    def get_volume_dicom(self, index): 
        dcm_file = self._paths[index] 
        f = dicom.read_file(dcm_file) 
        return f 

    def get_motion_quaternion(self, index): 
        motion_affine = self.get_motion_affine(index) 
        return quaternion_from_matrix(motion_affine)

    def get_motion_affine(self, index): 
        return self._motion[index] 

    def get_n_time_points(self): 
        return self._n_time_points 

    def get_duration(self,index): 
        return self.duration[index] 

    def load_data_files(self, path, from_dicom_comments=True, files_start_with=None, files_end_with=None, exclude_files_end_with=['.dat','.txt','.py','.pyc','.nii','.gz','.png','.jpg','.jpeg','.eps','.hdr','.l'] ):
        """Load vNAV dicom files from given path and extract motion information. 
        If from_dicom_comments==True, use the information stored in the dicom comments. 
        If from_dicom_comments==False, use the information stored in the dicom MOCO field. 
        As of April 2014, these two express motion in different coordinate systems. 
        The dicom comments report absolute motion in scanner coordinates (origin is the magnetic iso-center of the scanner).
        The dicom MOCO field is currently for proprietary use (??). """ 

        self._n_time_points     = 0 
        self._duration          = [] 
        self._motion            = [] 

        self._paths = []
        self._tx = []; self._ty = []; self._tz = []; 
        self._rx = []; self._ry = []; self._rz = []; 
        self._tx_comm = []; self._ty_comm = []; self._tz_comm = [];
        self._rx_comm = []; self._ry_comm = []; self._rz_comm = [];
        self._q0_comm = []; self._q1_comm = []; self._q2_comm = []; self._q3_comm = []; 
        self._a0_comm = []; self._a1_comm = []; self._a2_comm = []; self._a3_comm = []; 
        N=0 

        # pick the first dicom file found in path 
        files = os.listdir(path)
        files.sort() 

        # CASE 1: there exist files named with .dcm extension
        for file_name in files: 
            file_valid = True
            if files_start_with is not None: 
                if not file_name.startswith(files_start_with): 
                    file_valid = False
            if files_end_with is not None: 
                if not file_name.endswith(files_end_with): 
                    file_valid = False 
            for s in exclude_files_end_with: 
                if file_name.endswith(s): 
                    file_valid = False           
            if file_valid: 
                full_path = path+os.sep+file_name
                # read moco information from files 
                self._paths.append(full_path)
                try: 
                    f = dicom.read_file(full_path) 
                except: 
                    print "Could not read file ",full_path
                    return 
                t = f.get(0x00191025).value 
                r = f.get(0x00191026).value 
                self._tx.append(t[0]); self._ty.append(t[1]); self._tz.append(t[2]); 
                self._rx.append(r[0]); self._ry.append(r[1]); self._rz.append(r[2]); 
                motion_dicom_moco = []
                
                # extract moco information stored in the dicom comment field
                if from_dicom_comments: 
                    s = f.get(0x00204000).value 
                    if N: 
                        a = numpy.float32(s.split(' ')[1:5])
                        t = numpy.float32(s.split(' ')[6:9])
                        freq = numpy.float32(s.split(' ')[10])
                        r = angle_axis_to_rotation(a[0],a[1:4]) 
                    else: 
                        t = numpy.float32([0,0,0])
                        r = numpy.float32([0,0,0]) 
                        a = numpy.float32([0,1,0,0])  #FIXME: is this right?
                
                    q = angle_axis_to_quaternion(a.copy()[0],a.copy()[1:4]) 
                    self._a0_comm.append(a[0]); self._a1_comm.append(a[1]); self._a2_comm.append(a[2]); self._a3_comm.append(a[3]); 
                    self._tx_comm.append(t[0]); self._ty_comm.append(t[1]); self._tz_comm.append(t[2])
                    self._q0_comm.append(q[0]); self._q1_comm.append(q[1]); self._q2_comm.append(q[2]); self._q3_comm.append(q[3]); 
                    self._rx_comm.append(r[0]); self._ry_comm.append(r[1]); self._rz_comm.append(r[2]); 


                    tra_mat = tr.translation_matrix(t) 
                    rot_mat = tr.quaternion_matrix(q)
                    motion_dicom_comments = numpy.dot(tra_mat,rot_mat) 

                #xaxis, yaxis, zaxis = [1, 0, 0], [0, 1, 0], [0, 0, 1]
                #Rx = tr.rotation_matrix(r[0], xaxis)
                #Ry = tr.rotation_matrix(r[1], yaxis)
                #Rz = tr.rotation_matrix(r[2], zaxis)
                #rot_mat = tr.concatenate_matrices(Rx, Ry, Rz)
                #rot_mat = Ry.copy()
                #motion_dicom_comments = numpy.dot(tra_mat,rot_mat) 
                #motion_dicom_comments = rot_mat.copy()
                
                N += 1 
                if from_dicom_comments: 
                    self._motion.append(motion_dicom_comments) 
                else: 
                    self._motion.append(motion_dicom_moco) 
                acquisition_number = f.get(0x00200012).value 
                creation_time      = f.get(0x00080013).value
#                print "Acquisition number: ", acquisition_number
#                print "Creation time:      ",creation_time
        self._n_time_points = N

    def _draw_rectangle(self, axis, x, y, alpha=0.2, ec="gray", fc="gray"): #"CornflowerBlue"
        axis.add_patch( plt.Rectangle( (x[0],y[0]) , x[1]-x[0],y[1]-y[0], alpha=alpha, ec=ec, fc=fc, visible=True) )
        #plt.draw() 

    def _draw_rectangles(self, axis, windows, range_y): 
        for ii in range(len(windows)):
            tt = windows[ii]
            yy = (range_y[0],range_y[1])
            self._draw_rectangle(axis,tt,yy)

    def _draw_line(self, axis, x, range_y, color="#ff8d8d", linestyle="dashed", label=""): 
        axis.vlines(x, range_y[0], range_y[1], colors=color, linestyles=linestyle, label=label, visible=True)

    def _draw_events(self, axis, time, range_y, events, color="#ff8d8d", linestyle="dashed", label=""): 
        for t_index in time: 
            if t_index: 
                if events[t_index-1]: 
                    #print "Drawing line: ", t_index, range_y, color, linestyle
                    self._draw_line(axis, t_index, range_y, color, linestyle, label)

    def plot_motion(self, save_to_file=None, display_dicom_comments=True, display_dicom_moco=False, range_mm=(-8,8), range_deg=(-7,7), extract_events_threshold=THRESHOLD_MM, method='box', box_min=BOX_MIN, box_max=BOX_MAX, min_duration=1, line_color=LINE_COLOR, figsize=[10,5]):
        t = range(len(self._tx))

        # make windows:
        if extract_events_threshold is not None: 
            windows = []
            events  = self.extract_motion_events(method, extract_events_threshold, box_min, box_max)
        t_index_start = 0
        for t_index in t:  
            if t_index:     #this excludes the possibility of a motion event at time 0 
                if events[t_index-1]: 
                    t_index_end = t_index-1
                    if t_index_end-t_index_start > min_duration: 
                        windows.append( (t_index_start, t_index_end) )  #end window with frame before a motion event 
                    t_index_start = t_index+1                 #start window with frame after a motion event 
        windows.append( (t_index_start, t[-1]) ) 

        if display_dicom_moco:      
            fig1 = plt.figure(1,figsize=figsize)

            ax1 = fig1.add_subplot(311)
            ax1.plot(t, self._tx, line_color)
            ax1.grid(True)
            ax1.set_ylim( range_mm )
            ax1.set_ylabel('TX [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 

            ax1 = fig1.add_subplot(312)
            ax1.plot(t, self._ty, line_color)
            ax1.grid(True)
            ax1.set_ylim( range_mm )
            ax1.set_ylabel('TY [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 
    
            ax1 = fig1.add_subplot(313)
            ax1.plot(t, self._tz, line_color)
            ax1.grid(True)
            ax1.set_ylim( range_mm )
            ax1.set_ylabel('TZ [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 

#            if save_to_file is not None: 
#                plt.savefig(save_to_file)

            fig2 = plt.figure(2,figsize=figsize)
    
            ax1 = fig2.add_subplot(311)
            ax1.plot(t, rad_to_deg(self._rx), line_color)
            ax1.grid(True)
            ax1.set_ylim( range_deg )
            ax1.set_ylabel('RX [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
            
            ax1 = fig2.add_subplot(312)
            ax1.plot(t, rad_to_deg(self._ry), line_color)
            ax1.grid(True)
            ax1.set_ylim( range_deg )
            ax1.set_ylabel('RY [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
                                
            ax1 = fig2.add_subplot(313)
            ax1.plot(t, rad_to_deg(self._rz), line_color)
            ax1.grid(True)
            ax1.set_ylim( range_deg )
            ax1.set_ylabel('RZ [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
                                
#            if save_to_file is not None: 
#                plt.savefig(save_to_file)
        else: 
            fig1=None
            fig2=None

        if display_dicom_comments: 
            fig3 = plt.figure(3,figsize=figsize)

            #mr = numpy.min([self._rx_comm,self._ry_comm,self._rz_comm])
            #Mr = numpy.max([self._rx_comm,self._ry_comm,self._rz_comm])
            #mt = numpy.min([self._tx_comm,self._ty_comm,self._tz_comm])
            #Mt = numpy.max([self._tx_comm,self._ty_comm,self._tz_comm])
            mr = range_deg[0]
            Mr = range_deg[1]
            mt = range_mm[0]
            Mt = range_mm[1]

            ax1 = fig3.add_subplot(311)
            ax1.set_title("Rotation vs. vNAV frame number")
            ax1.plot(t, rad_to_deg(self._rx_comm), line_color)
            ax1.grid(True)
            ax1.set_ylim( (mr,Mr) )
            ax1.set_ylabel('RX [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            #print windows
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
    
            ax1 = fig3.add_subplot(312)
#            ax1.set_title("Rotation comments Y")
            ax1.plot(t, rad_to_deg(self._ry_comm), line_color)
            ax1.grid(True)
            ax1.set_ylim( (mr,Mr) )
            ax1.set_ylabel('RY [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
                
            ax1 = fig3.add_subplot(313)
#            ax1.set_title("Rotation comments Z")
            ax1.plot(t, rad_to_deg(self._rz_comm), line_color)
            ax1.grid(True)
            ax1.set_ylim( (mr,Mr) )
            ax1.set_ylabel('RZ [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_deg) 
                self._draw_events(ax1, t, range_deg, events) 
        
#            if save_to_file is not None: 
#                plt.savefig(save_to_file)

            fig4 = plt.figure(4,figsize=figsize)

            ax1 = fig4.add_subplot(311)
            ax1.set_title("Translation vs. vNAV frame number")
            ax1.plot(t, self._tx_comm, line_color)
            ax1.grid(True)
            ax1.set_ylim( (mt,Mt) )
            ax1.set_ylabel('TX [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 
                
            ax1 = fig4.add_subplot(312)
#            ax1.set_title("Translation comments Y")
            ax1.plot(t, self._ty_comm, line_color)
            ax1.grid(True)
            ax1.set_ylim( (mt,Mt) )
            ax1.set_ylabel('TY [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 
                
            ax1 = fig4.add_subplot(313)
#            ax1.set_title("Translation comments Z")
            ax1.plot(t, self._tz_comm, line_color)
            ax1.grid(True)
            ax1.set_ylim( (mt,Mt) )
            ax1.set_ylabel('TZ [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,range_mm) 
                self._draw_events(ax1, t, range_mm, events) 
                
#            if save_to_file is not None: 
#                plt.savefig(save_to_file)
        else: 
            fig3=None
            fig4=None
            
        plt.show()
        return fig1,fig2,fig3,fig4


    def get_mean_displacement(self, index, method='box', box_min=BOX_MIN, box_max=BOX_MAX): 
        mat = self.get_motion_affine(index)
        mat = Transform_Affine(mat) 
        if method=='box': 
            b = box_min
            B = box_max
            corners = numpy.asarray([ [b[0],b[1],b[2],1], [B[0],b[1],b[2],1], [b[0],B[1],b[2],1], [B[0],B[1],b[2],1], [b[0],b[1],B[2],1], [B[0],b[1],B[2],1], [b[0],B[1],B[2],1], [B[0],B[1],B[2],1] ]).transpose()
            corners_t = mat.left_multiply(corners)
            corners_t[3,:]=0
            corners[3,:]  =0
            #print "CORNERS: ", corners_t
            #dist = numpy.sqrt(((corners-corners_t)**2).sum(0))
            dist = (corners-corners_t).sum(0)
            mean_displ = numpy.mean(dist)
        else: 
            raise "Method to compute mean displacement is unknown. "
        return mean_displ
 
    def get_mean_displacement_variation(self, index, method='box', box_min=BOX_MIN, box_max=BOX_MAX):
        mat  = self.get_motion_affine(index)
        if index > 0: 
            mat0 = self.get_motion_affine(index-1)
        else: 
            mat0 = tr.identity_matrix()
        mat = Transform_Affine(mat) 
        mat0 = Transform_Affine(mat0)
        if method=='box': 
            b = box_min
            B = box_max
            corners = numpy.asarray([ [b[0],b[1],b[2],1], [B[0],b[1],b[2],1], [b[0],B[1],b[2],1], [B[0],B[1],b[2],1], [b[0],b[1],B[2],1], [B[0],b[1],B[2],1], [b[0],B[1],B[2],1], [B[0],B[1],B[2],1] ]).transpose()
            corners_t = mat.left_multiply(corners)
            corners_t[3,:]=0
            corners_t0 = mat0.left_multiply(corners)
            corners_t0[3,:]=0
            dist = numpy.sqrt(((corners_t-corners_t0)**2).sum(0))
            #dist = (corners-corners_t).sum(0)
            mean_displ = numpy.mean(dist)
        else: 
            raise "Method to compute mean displacement is unknown. "
        return mean_displ

    def get_mean_displacement_variation_since_time(self, index_new, index_old, method='box', box_min=BOX_MIN, box_max=BOX_MAX):
        mat  = self.get_motion_affine(index_new)
        mat0 = self.get_motion_affine(index_old)
        mat = Transform_Affine(mat) 
        mat0 = Transform_Affine(mat0)
        if method=='box': 
            b = box_min
            B = box_max
            corners = numpy.asarray([ [b[0],b[1],b[2],1], [B[0],b[1],b[2],1], [b[0],B[1],b[2],1], [B[0],B[1],b[2],1], [b[0],b[1],B[2],1], [B[0],b[1],B[2],1], [b[0],B[1],B[2],1], [B[0],B[1],B[2],1] ]).transpose()
            corners_t = mat.left_multiply(corners)
            corners_t[3,:]=0
            corners_t0 = mat0.left_multiply(corners)
            corners_t0[3,:]=0
            dist = numpy.sqrt(((corners_t-corners_t0)**2).sum(0))
            #dist = (corners-corners_t).sum(0)
            mean_displ = numpy.mean(dist)
        else: 
            raise "Method to compute mean displacement is unknown. "
        return mean_displ

        
    def extract_motion_events(self, method='box', threshold=THRESHOLD_MM, box_min=BOX_MIN, box_max=BOX_MAX):
        t = range(self.get_n_time_points())
        is_event   = numpy.zeros(len(t)-1)
        t_index_old = 0
        for t_index in t[1:]: 
#             mean_displ = self.get_mean_displacement_variation(t_index, method, box_min, box_max ) 
#             if numpy.sqrt((mean_displ)**2) >= threshold: 

             mean_displ = self.get_mean_displacement_variation_since_time(t_index, t_index_old, method, box_min, box_max )
             if numpy.sqrt((mean_displ)**2) >= threshold: 
                 t_index_old = numpy.copy(t_index)
                 is_event[t_index-1] = 1
             else: 
                 is_event[t_index-1] = 0 
        return is_event 


    def plot_mean_displacement(self, method='box', box_min=BOX_MIN,box_max=BOX_MAX, save_to_file=None, plot_zero=False, extract_events_threshold=THRESHOLD_MM, plot_range=[None,None], line_color=LINE_COLOR, min_duration=1, figsize=[10,5] ): 
        t = range(self.get_n_time_points())
        mean_displ = numpy.zeros(len(t))
        mean_displ_var = numpy.zeros(len(t))
        mean_displ_var_since_event = numpy.zeros(len(t))
        
        if extract_events_threshold is not None: 
            events = self.extract_motion_events(method, extract_events_threshold, box_min, box_max)

        t_index_old = 0
        for t_index in t: 
             mean_displ[t_index] = self.get_mean_displacement(t_index, method, box_min, box_max ) 
             mean_displ_var[t_index] = self.get_mean_displacement_variation(t_index, method, box_min, box_max ) 
             mean_displ_var_since_event[t_index] = self.get_mean_displacement_variation_since_time(t_index, t_index_old, method, box_min, box_max ) 
             if t_index: 
                 if events[t_index-1] == 1: 
                     t_index_old = t_index-1
        if not plot_zero: 
            t = t[1:]
            mean_displ = mean_displ[1:]
            mean_displ_var = mean_displ_var[1:]
            mean_displ_var_since_event = mean_displ_var_since_event[1:]

        # make windows:
        if extract_events_threshold is not None: 
            windows = []
            events  = self.extract_motion_events(method, extract_events_threshold, box_min, box_max)
        t_index_start = 0
        for t_index in t:  
            if t_index:     #this excludes the possibility of a motion event at time 0 
                if events[t_index-1]: 
                    t_index_end = t_index-1
                    if t_index_end-t_index_start > min_duration: 
                        windows.append( (t_index_start, t_index_end) )  #end window with frame before a motion event 
                    t_index_start = t_index+1                 #start window with frame after a motion event 
        windows.append( (t_index_start, t[-1]) ) 
        
#        mean_displ[numpy.where(mean_displ==0)]=-1000
#        mean_displ_var[numpy.where(mean_displ_var==0)]=-1000
#        mean_displ_var_since_event[numpy.where(mean_displ_var_since_event==0)]=-1000
        
        fig = plt.figure(5,figsize=figsize)

        ax1 = fig.add_subplot(311)
        ax1.set_title("Times frames - vNAV")
        ax1.plot(t,mean_displ, line_color)
        ax1.grid(True)
        if plot_range[0] is None: 
            pr0=copy.copy(mean_displ.min())
        else: 
            pr0=copy.copy(plot_range[0])
        if plot_range[1] is None: 
            pr1=copy.copy(mean_displ.max())
        else: 
            pr1=copy.copy(plot_range[1])
        ax1.set_ylim( [pr0,pr1] )
        ax1.set_ylabel('disp [mm]')
        if extract_events_threshold is not None: 
            ax1.hold(1)
            E = events*mean_displ
            E[numpy.where(E<0.1)]=-1000
            ax1.plot(t,E,'r.')
            #ax1.plot(t,0.5*(events*(mean_displ.max()-mean_displ.min())+mean_displ.min()),'r.')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')
        self._draw_rectangles(ax1,windows,[pr0,pr1]) 
        self._draw_events(ax1, t, [pr0,pr1], events) 


        ax1 = fig.add_subplot(312)
        #ax1.set_title("Mean displacement delta ")
        ax1.plot(t,mean_displ_var, line_color)
        ax1.grid(True)
        if plot_range[0] is None: 
            pr0=copy.copy(mean_displ_var.min())
        else: 
            pr0=copy.copy(plot_range[0])
        if plot_range[1] is None: 
            pr1=copy.copy(mean_displ_var.max())
        else: 
            pr1=copy.copy(plot_range[1])
        ax1.set_ylim( [pr0,pr1] )
        ax1.set_ylabel('delta [mm]')
        if extract_events_threshold is not None: 
            ax1.hold(1)
            E = events*mean_displ_var
            E[numpy.where(E<0.1)]=-1000
            ax1.plot(t,E,'r.')
            #ax1.plot(t,0.5*(events*(mean_displ.max()-mean_displ.min())+mean_displ.min()),'r.')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')
        self._draw_rectangles(ax1,windows,[pr0,pr1]) 
        self._draw_events(ax1, t, [pr0,pr1], events)             

        ax1 = fig.add_subplot(313)
        #ax1.set_title("Mean displacement event")
        ax1.plot(t,mean_displ_var_since_event, line_color)
        ax1.grid(True)
        if plot_range[0] is None: 
            pr0=copy.copy(mean_displ_var_since_event.min())
        else: 
            pr0=copy.copy(plot_range[0])
        if plot_range[1] is None: 
            pr1=copy.copy(mean_displ_var_since_event.max())
        else: 
            pr1=copy.copy(plot_range[1])
        ax1.set_ylim( [pr0,pr1] )
        ax1.set_ylabel('event [mm]')
        if extract_events_threshold is not None: 
            ax1.hold(1)
            E = events*mean_displ_var_since_event
            E[numpy.where(E<0.1)]=-1000
            ax1.plot(t,E,'r.')
            #ax1.plot(t,0.5*(events*(mean_displ.max()-mean_displ.min())+mean_displ.min()),'r.')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')
        self._draw_rectangles(ax1,windows,[pr0,pr1]) 
        self._draw_events(ax1, t, [pr0,pr1], events) 

        if save_to_file is not None: 
            plt.savefig(save_to_file)
        plt.show() 
        return fig

    def plot_quaternion(self, save_to_file=None, plot_range=[None,None], line_color=LINE_COLOR, figsize=[10,5]): 
        t = range(self.get_n_time_points())[1:]
        s = rad_to_deg(numpy.asarray(self._a0_comm))[1:]
        
        fig = plt.figure(6,figsize=figsize)
        ax1 = fig.add_subplot(211)
        ax1.set_title("Rotation agnle [deg] vs. vNAV frame number")
        ax1.plot(t,s, line_color)
        ax1.grid(True)
        if plot_range[0] is None: 
            plot_range[0]=s.min()
        if plot_range[1] is None: 
            plot_range[1]=s.max()
        ax1.set_ylim( plot_range )
        ax1.set_ylabel('Rotation angle [deg]')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')

        arc = numpy.zeros(self.get_n_time_points())
        v0  = numpy.asarray( self._a1_comm[0],self._a2_comm[0],self._a3_comm[0] )
        for t in range(self.get_n_time_points()): 
            vt  = numpy.asarray( self._a1_comm[t],self._a2_comm[t],self._a3_comm[t] )
            arc[t] = numpy.dot( numpy.transpose(v0), vt )
        ax1 = fig.add_subplot(212)
        ax1.set_title("Arc vs. vNAV frame number")
        t = range(self.get_n_time_points())
        ax1.plot(t,arc, line_color)
        ax1.grid(True)
        plot_range[0]=arc.min()-0.0001
        plot_range[1]=arc.max()+0.0001
        ax1.set_ylim( plot_range )
        ax1.set_ylabel('Arc [steradians]')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')
    
        if save_to_file is not None: 
            plt.savefig(save_to_file)
        plt.show() 
        return fig


    def _repr_html_(self): 
        self.plot_mean_displacement() 

    
def load_vnav_mprage(path, from_dicom_comments=True, files_start_with=None, files_end_with=None): 
    return vNAV_MPRage(path, from_dicom_comments, files_start_with, files_end_with)
