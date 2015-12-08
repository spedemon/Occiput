
# occiput 
# March 2015 
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA 


__all__ = ['load_motion_sensor_data'] 

from occiput.Core import Transform_Affine 
import numpy 
from numpy import float32, int32 
import matplotlib.pyplot as plt
import copy 
from occiput.Core import transformations as tr

BOX_MIN      = [-50.0,-50.0,-50.0]
BOX_MAX      = [50.0,50.0,50.0]
THRESHOLD_MM = 5.0
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





class Motion_Sensor: 
    def __init__(self,filename=None, channels=['B']): 
        self._reset()
        if filename is not None: 
            self.load_from_file(filename, channels) 

    def get_motion_quaternion(self, index): 
        motion_affine = self.get_motion_affine(index) 
        return quaternion_from_matrix(motion_affine)

    def get_motion_affine(self, index): 
        return self._motion[index] 

    def _reset(self):
        self._motion = []
        self._n_time_points = 0 
        self._tx = []; self._ty = []; self._tz = []; 
        self._rx = []; self._ry = []; self._rz = []; 
        self._q0 = []; self._q1 = []; self._q2 = []; self._q3 = [];

    def get_n_time_points(self):
        return self._n_time_points

    def load_from_file(self,filename, channels=['B']): 
        f = open(filename)
        self._reset() 
        while 1: 
            l = f.readline() 
            if l=='': 
                break 
            d = l.split() 
            channel = d[0][1]
            
            do_process = False 
            if channels is None: 
                do_process = True
            elif channels == []: 
                do_process = True 
            elif channels == '': 
                do_process = True 
            else: 
                if channel in channels: 
                    do_process = True
                else: 
                    do_process = False

            if do_process: 
                data    = d[1]
                x  = float32(d[2])
                y  = float32(d[3])
                z  = float32(d[4])
                q0 = float32(d[5])
                q1 = float32(d[6])
                q2 = float32(d[7])
                q3 = float32(d[8])
                q = [q0,q1,q2,q3]
                #q = tr.quaternion_conjugate(q)
                status = int32(d[11])
                uncertainty = float32(d[12]) 
                self._n_time_points += 1
                tra_mat = tr.translation_matrix([x,y,z]) 
                rot_mat = tr.quaternion_matrix(q)
                self._motion.append( numpy.dot(tra_mat,rot_mat) )
                rotation = quaternion_to_rotation(q)
                self._tx.append(x)
                self._ty.append(y)
                self._tz.append(z)
                self._rx.append(rotation[0])
                self._ry.append(rotation[1])
                self._rz.append(rotation[2])
                self._q0.append(q[0])
                self._q1.append(q[1])
                self._q2.append(q[2])
                self._q3.append(q[3])

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
              if t_index<self.get_n_time_points():
                if events[t_index-1]: 
                    #print "Drawing line: ", t_index, range_y, color, linestyle
                    self._draw_line(axis, t_index, range_y, color, linestyle, label)

    def plot_motion(self, save_to_file=None, extract_events_threshold=THRESHOLD_MM, method='box', box_min=BOX_MIN, box_max=BOX_MAX, min_duration=10, line_color=LINE_COLOR):
        t = range(len(self._tx))

        # make windows:
        if extract_events_threshold is not None: 
            windows = []
            events  = self.extract_motion_events(method, extract_events_threshold, box_min, box_max)
        t_index_start = 0
        for t_index in t:  
            if t_index:     #this excludes the possibility of a motion event at time 0 
              if t_index<self.get_n_time_points():
                if events[t_index-1]: 
                    t_index_end = t_index-1
                    if t_index_end-t_index_start > min_duration: 
                        windows.append( (t_index_start, t_index_end) )  #end window with frame before a motion event 
                    t_index_start = t_index+1                 #start window with frame after a motion event 
        windows.append( (t_index_start, t[-1]) ) 

        if 1:      
            fig1 = plt.figure(1,figsize=(17.5,4),dpi=200)

            ax1 = fig1.add_subplot(321)
            ax1.plot(t, self._tx, line_color)
            ax1.grid(True)
            pr0 = float32(self._tx).min()
            pr1 = float32(self._tx).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('TX [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 

            ax1 = fig1.add_subplot(323)
            ax1.plot(t, self._ty, line_color)
            ax1.grid(True)
            pr0 = float32(self._ty).min()
            pr1 = float32(self._ty).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('TY [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 
    
            ax1 = fig1.add_subplot(325)
            ax1.plot(t, self._tz, line_color)
            ax1.grid(True)
            pr0 = float32(self._tz).min()
            pr1 = float32(self._tz).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('TZ [mm]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 

#            if save_to_file is not None: 
#                plt.savefig(save_to_file)

            fig2 = fig1 #plt.figure(2)
    
            ax1 = fig2.add_subplot(322)
            ax1.plot(t, rad_to_deg(self._rx), line_color)
            ax1.grid(True)
            pr0 = float32(rad_to_deg(self._rx)).min()
            pr1 = float32(rad_to_deg(self._rx)).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('RX [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 
            
            ax1 = fig2.add_subplot(324)
            ax1.plot(t, rad_to_deg(self._ry), line_color)
            ax1.grid(True)
            pr0 = float32(rad_to_deg(self._ry)).min()
            pr1 = float32(rad_to_deg(self._ry)).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('RY [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 
                                
            ax1 = fig2.add_subplot(326)
            ax1.plot(t, rad_to_deg(self._rz), line_color)
            ax1.grid(True)
            pr0 = float32(rad_to_deg(self._rz)).min()
            pr1 = float32(rad_to_deg(self._rz)).max()
            d = pr1-pr0
            pr = [pr0-d/4.0, pr1+d/4.0]
            ax1.set_ylim( pr )
            ax1.set_ylabel('RZ [deg]')
            #for label in ax1.get_xticklabels():
            #    label.set_color('r')
            if extract_events_threshold  is not None:
                self._draw_rectangles(ax1,windows,pr) 
                self._draw_events(ax1, t, pr, events) 
                                
#            if save_to_file is not None: 
#                plt.savefig(save_to_file)
            
        plt.show()
        #return fig1,fig2


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

        
    def extract_motion_events(self, method='box', threshold=THRESHOLD_MM, box_min=BOX_MIN, box_max=BOX_MAX, prune_distance=50):
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
        if prune_distance >=2: 
            last_event = 0
            for i in range(len(is_event)): 
                if is_event[i]: 
                    if (i-last_event) <= prune_distance:
                        is_event[i]=0
                    else:     
                        last_event = i
        return is_event 


    def plot_mean_displacement(self, method='box', box_min=BOX_MIN,box_max=BOX_MAX, save_to_file=None, plot_zero=False, extract_events_threshold=THRESHOLD_MM, plot_range=[None,None], line_color=LINE_COLOR, min_duration=10 ): 
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
        
        fig = plt.figure(5,figsize=(8,4),dpi=200)

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
        #return fig

    def plot_quaternion(self, save_to_file=None, line_color=LINE_COLOR): 
        t = range(self.get_n_time_points())[1:]
        s = rad_to_deg(numpy.asarray(self._q0))[1:]
        
        fig = plt.figure(6,figsize=(8,4),dpi=200)
        ax1 = fig.add_subplot(211)
        ax1.set_title("Rotation agnle [deg] vs. vNAV frame number")
        ax1.plot(t,s, line_color)
        ax1.grid(True)
        pr0 = s.min()
        pr1 = s.max() 
        d = pr1-pr0
        pr = [pr0-d/4.0, pr1+d/4.0]
        ax1.set_ylim( pr )
        ax1.set_ylabel('Rotation angle [deg]')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')

        arc = numpy.zeros(self.get_n_time_points()) 
        v0  = numpy.asarray( self._q1[0],self._q2[0],self._q3[0] ) 
        for t in range(self.get_n_time_points()): 
            vt  = numpy.asarray( self._q1[t],self._q2[t],self._q3[t] ) 
            arc[t] = numpy.dot( numpy.transpose(v0), vt ) 
        ax1 = fig.add_subplot(212) 
        ax1.set_title("Arc vs. vNAV frame number") 
        t = range(self.get_n_time_points()) 
        ax1.plot(t,arc, line_color)
        ax1.grid(True)
        pr0 = arc.min()
        pr1 = arc.max() 
        d = pr1-pr0
        pr = [pr0-d/4.0, pr1+d/4.0]
        ax1.set_ylim( pr )
        ax1.set_ylabel('Arc [steradians]')
        #for label in ax1.get_xticklabels():
        #    label.set_color('r')
    
        if save_to_file is not None: 
            plt.savefig(save_to_file)
        plt.show() 
        #return fig

    def _repr_html_(self): 
        self.plot_mean_displacement() 



def load_motion_sensor_data(filename, channels=['B',]): 
    return Motion_Sensor(filename, channels)
    


