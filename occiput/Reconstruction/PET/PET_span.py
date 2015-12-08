
# occiput.io 
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# 2013 - 2015
# Boston, MA


# Occiput.io does not like the concepts of sinogram-span and meshing. 
# These concepts have been utilized in the past in order to reduce the computational 
# complexity of PET reconstruction. However sinogram meshing makes the emission data 
# complicated and the complicates the reconstruction software. 
# Occiput.io can do without. It uses the GPU acceleration to avoid meshing altogether and 
# it also provides a technique similar to meshing to further speed-up the reconstruction. 
# This is full automated, it does not require the concept of span and Michelogram.  
# For compatibility, however, occiput can deal with the conventional meshing technique. 
# The code for conventional meshing is here. 



from occiput.Visualization import svgwrite, has_svgwrite
from numpy import int32, ones 


class Michelogram: 
    def __init__(self, n_rings=64, span=11, max_ring_difference=60): 
        self.n_rings = n_rings
        self.span = span
        self.max_ring_difference = max_ring_difference 
        self.n_sinograms = None 
        self.n_planes    = None 
        self._svg_string = None 
        self._make()
        self._make_svg()

    def _make(self): 
        if not self._verify_max_ring_difference(): 
            print "Michelogram: Maximum ring difference is not consistent with span. "
        n_sinograms_per_side = (self.max_ring_difference+1-self.span/2)/self.span 
        n_sinograms = n_sinograms_per_side * 2 + 1
        n_planes = zeros([n_sinograms_per_side+1,],dtype=int32)
        n_planes[0] = 2*n_rings-1
        n_planes[1] = n_planes[0] - self.span
        for i in range(n_sinograms_per_side-1): 
            n_planes[2+i]=n_planes[1+i]-2*span
        self.n_sinograms = n_sinograms
        self.n_planes    = n_planes

        self.michelogram_sinogram = -1*ones([n_rings,n_rings],order = 'F',dtype=int32)
        self.michelogram_plane    = -1*ones([n_rings,n_rings],order = 'F',dtype=int32)

        # -2- right sinograms  
        sinogram_index = 0
        for X in range(0,max_ring_diff,span): 
            if sinogram_index == 0: 
                plane_index_offset = -n_rings_small
            else: 
                plane_index_offset = 0
            self._fill_michelogram(X, -n_rings_small, sinogram_index, plane_index_offset)
            sinogram_index += 2

        # -3- left sinograms 
        sinogram_index = 1
        for Y in range(n_rings_large,max_ring_diff,span): 
            plane_index_offset = 0
            self._fill_michelogram(0, Y, sinogram_index, plane_index_offset)
            sinogram_index += 2
            
    def _verify_max_ring_difference(self): 
        return True
        #verify if max_ring_difference is consistent with span and n_rings (give warnings if not): 

    def _fill_michelogram(self, x_start, y_start, sinogram_index, plane_index_offset): 
        n_rings = self.n_rings
        n_rings_small = self.span/2
        n_rings_large = n_rings_small+1 
        for X in range(n_rings + n_rings_small): 
            # long line
            plane = 2*X + plane_index_offset
            for line_index in range(n_rings_large): 
                x = x_start + X - line_index
                y = y_start + X + line_index
                if x>=0 and x<n_rings and y>=0 and y<n_rings: 
                    self.michelogram_sinogram[x,y] = sinogram_index
                    self.michelogram_plane[x,y] = plane
            # short line
            plane = 2*X + 1 + plane_index_offset
            for line_index in range(n_rings_small): 
                x = x_start + X - line_index
                y = y_start + X + 1 + line_index
                if x>=0 and x<n_rings and y>=0 and y<n_rings: 
                    self.michelogram_sinogram[x,y] = sinogram_index
                    self.michelogram_plane[x,y] = plane 

    def _make_svg(self): 
        if not has_svgwrite: 
            self._svg_string = None 
            return self._svg_string 
        w = '100%'
        h = '100%'
        dwg = svgwrite.Drawing('michelogram.svg',size=(w,h), profile='full', debug=True)
        dwg.viewbox(width=100, height=100)

        # blank table 
        rect = dwg.add(dwg.rect(insert=(0, 0), size=(100, 100), rx=0.0, ry=0.0))
        rect.fill('grey',opacity=0.1).stroke('black',width=0.1,opacity=0.001)
        #for ... 
        #line = dwg.add(dwg.polyline([(10, 10), (10, 100), (100, 100), (100, 10), (10, 10)],stroke='black', fill='none'))

        # fill table with Michelogram data 

        self._svg_string = dwg.tostring() 
        return self._svg_string

    def _repr_svg_(self): 
        return self._svg_string    

    def display_image(self): 
        import pylab 
        pylab.imshow(self.michelogram_plane,interpolation='nearest',cmap='gray') 


