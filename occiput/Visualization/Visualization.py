
# occiput 
# Stefano Pedemonte 
# Oct 2015
# Harvard University, Martinos Center for Biomedical Imaging 
# Boston, MA, USA 


from PIL import Image
import numpy
import uuid

from occiput.global_settings import is_gpu_enabled
from DisplayNode import DisplayNode
import Colors as C
from IPython.display import HTML, display, Javascript
from ipynb import is_in_ipynb



class InstallationError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

        
    
    
    
class ProgressBar(): 
    def __init__(self, height='6px', width='100%%', background_color=C.LIGHT_BLUE, foreground_color=C.BLUE): 
        self._percentage = 0.0
        self._divid = str(uuid.uuid4()) 
        self.visible = False
        if is_in_ipynb(): 
            self.set_display_mode("ipynb")
        else: 
            self.set_display_mode("text")
        self._pb = HTML(
            """
            <div style="border: 1px solid white; width:%s; height:%s; background-color:%s">
                <div id="%s" style="background-color:%s; width:0%%; height:%s"> </div>
            </div> 
            """ % ( width, height, background_color, self._divid, foreground_color, height))

    def show(self): 
        if self.mode == "ipynb": 
            display(self._pb)
        self.visible = True 
        
    def set_display_mode(self, mode="ipynb"): 
        self.mode=mode 

    def set_percentage(self,percentage):
        if not self.visible: 
            self.show()
        if percentage < 0.0:
            percentage = 0.0
        if percentage > 100.0:
            percentage = 100.0 
        percentage = int(percentage)
        self._percentage = percentage
        if self.mode == "ipynb": 
            display(Javascript("$('div#%s').width('%i%%')" % (self._divid, percentage)))
        else: 
            print "%2.1f / 100"%percentage

    def get_percentage(self):
        return self._percentage







class MultipleVolumes(): 
    def __init__(self,volumes,axis=0,shrink=256,rotate=90,subsample_slices=None,scales=None,open_browser=None): 
        self.volumes = volumes 
        self._axis = axis
        self._shrink = shrink
        self._rotate = rotate
        self._subsample_slices = subsample_slices
        self._scales = scales
        self._open_browser = open_browser
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=C.LIGHT_GRAY, foreground_color=C.GRAY)

    def get_data(self,volume_index): 
        volume = self.volumes[volume_index] 
        if isinstance(volume,numpy.ndarray): 
            return volume
        else:
            # Volume type
            return volume.data
        
    def get_shape(self,volume_index): 
        return self.volumes[volume_index].shape
        
    def export_image(self,volume_index,slice_index,axis=0,normalise=True,scale=None,shrink=None,rotate=0,global_scale=True): 
        #FIXME: handle 4D, 5D, ..
        
        M = self.get_data(volume_index).max()
        m = self.get_data(volume_index).min()
        
        if axis == 0:
            a = numpy.float32(self.get_data(volume_index)[slice_index,:,:].reshape((self.get_shape(volume_index)[1],self.get_shape(volume_index)[2]))) 
        elif axis == 1: 
            a = numpy.float32(self.get_data(volume_index)[:,slice_index,:].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[2])))
        else: 
            a = numpy.float32(self.get_data(volume_index)[:,:,slice_index].reshape((self.get_shape(volume_index)[0],self.get_shape(volume_index)[1])))

        if m>=0: 
            if normalise: 
                if not global_scale: 
                    if scale is None:
                        a = a * 255/(a.max()+1e-9) 
                    else: 
                        a = a * scale *255/(a.max()+1e-9) 
                else: 
                    if scale is None:
                        a = a * 255/(M+1e-9) 
                    else: 
                        a = a * scale *255/(M+1e-9) 
            im = Image.fromarray(a).convert("RGB") 
        else:
            if normalise: 
                if not global_scale:
                    if scale is None:
                        a = a * 512/(a.max()-a.min()+1e-9) 
                    else: 
                        a = a * scale * 512/(a.max()-a.min()+1e-9) 
                else: 
                    if scale is None:
                        a = a * 512/(M-m+1e-9) 
                    else: 
                        a = a * scale * 512/(M-m+1e-9)                     
            blue = a
            red  = a.copy()
            red[numpy.where(red>=0)]  = 0 
            red  = - red
            blue[numpy.where(blue<0)] = 0 
            green = numpy.zeros(red.shape)
            rgb = numpy.zeros((red.shape[0],red.shape[1],3),dtype=numpy.uint8)
            rgb[:,:,0]=red
            rgb[:,:,1]=green
            rgb[:,:,2]=blue
            im = Image.fromarray(rgb,mode='RGB')
        if shrink  is not None: 
            # scale in order to make the largest dimension equal to 'shrink' 
            shrink = int(shrink) 
            (h,w) = im.size
#FIXME: downsample the volume (with the GPU, if available) before converting to images, it will save a lot of time, conversion to Image and to RGB is slow
            if h>shrink or w>shrink: 
                if h>w: 
                    im = im.resize( (shrink, int(shrink*w/(1.0*h))) )
                else: 
                    im = im.resize( (int(shrink*h/(1.0*w)),shrink) )
        if rotate is not None: 
            im = im.rotate(rotate)
        return im
        

    # display
    def display_in_browser(self,axis=None,shrink=None,rotate=None,subsample_slices=None,scales=None): 
        self.display(axis,shrink,rotate,subsample_slices,scales,open_browser=True)

    def display(self,axis=None,shrink=None,rotate=None,subsample_slices=None,scales=None,open_browser=None):
        if axis is None:
            axis = self._axis
        if shrink  is None:
            shrink = self._shrink
        if rotate  is None: 
            rotate = self._rotate
        if  subsample_slices  is None: 
            subsample_slices = self._subsample_slices
        if subsample_slices is None: 
            subsample_slices=1
        if scales  is None: 
            scales = self._scales 
        if open_browser  is None: 
            open_browser = self._open_browser
        if open_browser  is None: 
            open_browser = False 
        D = DisplayNode()
        images = []
        n=0
        
        N_slices = 0
        for j in range(len(self.volumes)): 
            N_slices += self.get_shape(j)[axis] 
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=C.LIGHT_GRAY, foreground_color=C.GRAY)
        
        for j in range(len(self.volumes)): 
            if scales  is None: 
                scale = 255/(self.get_data(j).max()+1e-9)
            else: 
                scale = scales[j]*255/(self.get_data(j).max()+1e-9)
            images_inner = []
            for i in range(0,self.get_shape(j)[axis],subsample_slices): 
                im = self.export_image(j,i,axis,normalise=True,scale=scale,shrink=shrink,rotate=rotate)
                images_inner.append(im) 
                n+=1
                self._progress_bar.set_percentage(n*100.0/N_slices*subsample_slices) 
            images.append(images_inner) 
        return D.display('tipix', images, open_browser)   

    def _repr_html_(self): 
        return self.display()._repr_html_()
        






class MultipleVolumesNiftyPy(): 
    def __init__(self,volumes,axis=0,open_browser=None): 
        self.volumes = volumes 
        self._axis = axis
        self._open_browser = open_browser
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=C.LIGHT_GRAY, foreground_color=C.GRAY)

    def get_data(self,volume_index): 
        volume = self.volumes[volume_index] 
        if isinstance(volume,numpy.ndarray): 
            return volume
        else:
            # Image3D
            return volume.data
        
    def get_shape(self,volume_index): 
        return self.volumes[volume_index].shape
        
    def _resample_volume(self): 
        pass 

    # display
    def display_in_browser(self,axis=None,max_size=200): 
        self.display(axis,max_size,open_browser=True)

    def display(self,axis=None, max_size=256, open_browser=None):
        if axis is None:
            axis = self._axis
        if open_browser  is None: 
            open_browser = self._open_browser
        if open_browser  is None: 
            open_browser = False 
        D = DisplayNode() 
        
        self._progress_bar = ProgressBar(height='6px', width='100%%', background_color=C.LIGHT_GRAY, foreground_color=C.GRAY)
        
        volume = self.volumes[0]    #FIXME: optionally use other grid
        # make grid of regularly-spaced points 
        box_min = volume.get_world_grid_min() 
        box_max = volume.get_world_grid_max() 
        span = box_max-box_min 
        max_span = span.max() 
        n_points = numpy.uint32(span / max_span * max_size) 
        grid = volume.get_world_grid(n_points) 
        n=0
        images = [] 
        for index in range(len(self.volumes)): 
                volume = self.volumes[index] 
                resampled = volume.compute_resample_on_grid(grid).data
                self._resampled = resampled
                sequence = []
                for slice_index in range(n_points[axis]): 
                    if axis == 0:
                        a = numpy.float32(resampled[slice_index,:,:].reshape((resampled.shape[1],resampled.shape[2]))) 
                    elif axis == 1: 
                        a = numpy.float32(resampled[:,slice_index,:].reshape((resampled.shape[0],resampled.shape[2])))
                    else: 
                        a = numpy.float32(resampled[:,:,slice_index].reshape((resampled.shape[0],resampled.shape[1])))
                    lookup_table = volume.get_lookup_table()
                    im = self.__array_to_im( a, lookup_table )
                    sequence.append(im.rotate(90))  #FIXME: make optional 
                    n+=1
                    self._progress_bar.set_percentage(n*100.0/(len(self.volumes)*max_size)) 
                images.append(sequence) 
        if len(images)==1: 
            return D.display('tipix', images[0], open_browser) 
        else:  
            return D.display('tipix', images, open_browser)   

    def __array_to_im(self, a, lookup_table): 
        if lookup_table  is not None: 
            red,green,blue,alpha = lookup_table.convert_ndarray_to_rgba(a)
            rgb = numpy.zeros((a.shape[0],a.shape[1],3),dtype=numpy.uint8)
            rgb[:,:,0]=red
            rgb[:,:,1]=green
            rgb[:,:,2]=blue
            im = Image.fromarray(rgb,mode='RGB') 
        else: 
            im = Image.fromarray(a).convert("RGB") 
        return im 
        
    def _repr_html_(self): 
        return self.display()._repr_html_()








def deg_to_rad(deg):
    return deg*numpy.pi/180.0

def rad_to_deg(rad):
    return rad*180.0/numpy.pi




try: 
    from NiftyPy.NiftyRec import SPECT_project_parallelholes as projection
except: 
    has_NiftyPy = False
    print "Please install NiftyPy"
else: 
    has_NiftyPy = True





class SPECT_Projection(): 
    """SPECT projection object. """
    def __init__(self, data): 
        self.data = data        
    
    def get_data(self):
        """Returns the raw projection data (note that is can be accessed also as self.data ). """
        return data

    def save_to_file(self, filename): 
        h5f = h5py.File(filename, 'w')
        h5f.create_dataset('data', data=self.data)
        #h5f.create_dataset('size_x', data=size_x)
        #h5f.create_dataset('size_y', data=size_y)
        h5f.close()

    def get_integral(self): 
        return self.data.sum() 

    def to_image(self,data,index=0,scale=None,absolute_scale=False): 
        from PIL import Image
        a = numpy.float32(data[:,:,index].reshape((data.shape[0],data.shape[1])))   
        if scale is None:
            a = 255.0*(a)/(a.max()+1e-12)
        else: 
            if absolute_scale: 
                a = scale*(a) 
            else: 
                a = scale*255.0*(a)/(a.max()+1e-12)
        return Image.fromarray(a).convert("RGB")
        
    def display_in_browser(self,axial=True,azimuthal=False,index=0,scale=None): 
        self.display(axial=axial,azimuthal=azimuthal,index=index,scale=scale,open_browser=True)

    def display(self,scale=None,open_browser=False): 
        data = self.data
        d = DisplayNode()
        images = []
        progress_bar = ProgressBar(height='6px', width='100%%', background_color=C.LIGHT_GRAY, foreground_color=C.GRAY) 
        if scale is not None:
            scale = scale*255.0/(data.max()+1e-12)
        else: 
            scale = 255.0/(data.max()+1e-12) 
        N_projections = self.data.shape[2]
        N_x = self.data.shape[0]
        N_y = self.data.shape[1]
        print "SPECT Projection   [N_projections: %d   N_x: %d   N_y: %d]"%(N_projections,N_x,N_y)
        for i in range( N_projections ): 
                images.append( self.to_image(data,i,scale=scale,absolute_scale=True) ) 
                progress_bar.set_percentage(i*100.0/N_projections)                         
        progress_bar.set_percentage(100.0)
        return d.display('tipix', images, open_browser)
        
    def _repr_html_(self): 
        return self.display()._repr_html_()
        

class VolumeRenderer(): 
    def __init__(self, volume, theta_min_deg=0.0, theta_max_deg=360.0, n_positions=180, truncate_negative=False, psf=None, attenuation=None, scale=1.0): 
        self.volume  = volume 
        self.psf = psf
        self.attenuation = attenuation
        self.use_gpu = is_gpu_enabled() 
        self.theta_min_deg = theta_min_deg 
        self.theta_max_deg = theta_max_deg 
        self.n_positions = n_positions
        self.truncate_negative = truncate_negative
        self.scale=scale

    def __make_cameras(self,axis,direction): 
        #self.cameras = numpy.float32(numpy.linspace(deg_to_rad(self.theta_min_deg),deg_to_rad(self.theta_max_deg),self.n_positions).reshape((self.n_positions,1))) 
        self.cameras = numpy.zeros((self.n_positions,3),dtype=numpy.float32)
        self.cameras[:,0] = self.cameras[:,0] + deg_to_rad(axis[0])
        self.cameras[:,1] = self.cameras[:,0] + deg_to_rad(axis[1])
        self.cameras[:,2] = self.cameras[:,0] + deg_to_rad(axis[2])
        self.cameras[:,direction] = numpy.linspace(deg_to_rad(self.theta_min_deg),deg_to_rad(self.theta_max_deg),self.n_positions)
        
    def render(self, axis=[90,0,0], direction=0, max_n_points=None): 
        if hasattr(self.volume,'compute_resample_on_grid'):   #i.e. if it is a Image3 
            volume = self.volume.copy()
            # make grid of regularly-spaced points
            if max_n_points is None: 
                max_n_points = 256
            box_min = volume.get_world_grid_min() 
            box_max = volume.get_world_grid_max() 
            span = box_max-box_min 
            max_span = span.max() 
            n_points = numpy.uint32(span / max_span * max_n_points) 
            grid = volume.get_world_grid(n_points) 
            volume.compute_resample_on_grid(grid)
            volume = numpy.float32(volume.data)
        else: 
            volume = self.volume
        self.__make_cameras(axis, direction)
        if has_NiftyPy: 
            proj_data = projection(volume, self.cameras, self.attenuation, self.psf, 0.0, 0.0, self.use_gpu, self.truncate_negative)
        else: 
            raise InstallationError("NiftyPy not installed, please install to execute render(). ")
        self.__proj = SPECT_Projection(proj_data)
        return self.__proj  #FIXME: memoize projection (use new style objects - properties)

    def display(self): 
        r = self.render()
        return r.display(scale=self.scale) 

    def _repr_html_(self):
        return self.display()





H_graph = """
<div id="__id__"></div>
<style type="text/css">
path.link {
  fill: none;
  stroke: #666;
  stroke-width: 1.5px;
}

marker#t0 {
  fill: green;
}

path.link.t0 {
  stroke: green;
}

path.link.t2 {
  stroke-dasharray: 0,2 1;
}

circle {
  fill: #ccc;
  stroke: #333;
  stroke-width: 1.5px;
}

text {
  font: 10px sans-serif;
  pointer-events: none;
}

text.shadow {
  stroke: #fff;
  stroke-width: 3px;
  stroke-opacity: .8;
}
</style>
"""

J_graph = """  
<script type="text/Javascript">
require.config({paths: {d3: "http://d3js.org/d3.v3.min"}});
require(["d3"], function(d3) {

var width = 800;
var height = 500;

var graph = __graph_data__; 

var nodes = graph['nodes']; 
var links = graph['links']; 

var nodes2 = {};
var links2 = []; 

for (var i=0; i<nodes.length; i++) {
    node = nodes[i]; 
    nodes2[node['name']] = node; 
};

for (var i=0; i<links.length; i++) {
    links2[i] = {'source':nodes2[links[i]['source']], 'target':nodes2[links[i]['target']], 'type':links[i]['type'],}; 
};


var force = d3.layout.force()
    .nodes(d3.values(nodes2))
    .links(links2)
    .size([width, height])
    .linkDistance(60)
    .charge(-300)
    .on("tick", tick)
    .start();

var svg = d3.select("#__id__").append("svg:svg")
    .attr("width", width)
    .attr("height", height);

// Per-type markers, as they don't inherit styles.
svg.append("svg:defs").selectAll("marker")
    .data(["t0", "t1", "t2"])
  .enter().append("svg:marker")
    .attr("id", String)
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 20)
    .attr("refY", -1.5)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
  .append("svg:path")
    .attr("d", "M0,-5L10,0L0,5");

var path = svg.append("svg:g").selectAll("path")
    .data(force.links())
  .enter().append("svg:path")
    .attr("class", function(d) { return "link " + d.type; })
    .attr("marker-end", function(d) { return "url(#" + d.type + ")"; });

colors  = ['#fff','#ccc','#ee4','#dddddd','#eeeeee','#ffffc0']; 
strokes = ['#333','#333','#333','#dddddd','#dddddd','#dddddd']; 

var circle = svg.append("svg:g").selectAll("circle")
    .data(force.nodes())
  .enter().append("svg:circle")
    .attr("r", 10)
    .call(force.drag)
    .style('fill', function(d){return colors[d.type];})
    .style('stroke', function(d){return strokes[d.type];});

var text = svg.append("svg:g").selectAll("g")
    .data(force.nodes())
  .enter().append("svg:g");

// A copy of the text with a thick white stroke for legibility.
text.append("svg:text")
    .attr("x", 12)
    .attr("y", ".31em")
    .attr("class", "shadow")
    .text(function(d) { return d.name; });

text.append("svg:text")
    .attr("x", 12)
    .attr("y", ".31em")
    .text(function(d) { return d.name; });


// Use elliptical arc path segments to doubly-encode directionality.
function tick() {
  path.attr("d", function(d) {
    var dx = d.target.x - d.source.x,
        dy = d.target.y - d.source.y,
        dr = Math.sqrt(dx * dx + dy * dy);
    return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;
  });

  circle.attr("transform", function(d) {
    return "translate(" + d.x + "," + d.y + ")";
  });

  text.attr("transform", function(d) {
    return "translate(" + d.x + "," + d.y + ")";
  });
}

});
</script>
"""


import json 

class Graph(): 
    def __init__(self, graph): 
        self.set_graph(graph)

    def set_graph(self,g):
        self.graph = g

    def get_html(self): 
        div_id = "graph_"+str(uuid.uuid4())
        H = H_graph.replace("__id__",div_id)
        J = J_graph.replace("__id__",div_id)
        J = J.replace("__graph_data__",json.dumps(self.graph))
        return H+J

    def _repr_html_(self): 
        graph = HTML( self.get_html() )
        return graph._repr_html_() 

    


## ipy_table 
#FIXME: perhaps move this somewhere else

try: 
    import ipy_table 
    has_ipy_table = True
except: 
    print "Please install ipy_table (e.g. 'easy_install ipy_table') to enable ipython notebook tables. "
    ipy_table = None 
    has_ipy_table = False




## svg_write
#FIXME: perhaps move this somewhere else

try: 
    import svgwrite
    has_svgwrite = True
except: 
    print "Please install svgwrite (e.g. 'easy_install svgwrite') to enable svg visualisations. "
    svgwrite = None
    has_svgwrite = False




        

    
