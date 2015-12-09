
# occiput 
# Stefano Pedemonte
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Sep. 2015, Boston 


__all__ = ["Michelogram"] 


from occiput.Visualization import ipy_table, has_ipy_table, svgwrite, has_svgwrite 
from PIL import ImageDraw
from PIL import Image as PIL 
from numpy import isscalar, linspace, int32, uint32, ones, zeros, pi, sqrt, float32, float64, where, ndarray, nan
from numpy import inf, asarray, concatenate, fromfile, maximum, exp, asfortranarray, fliplr, transpose 


DEFAULT_SPAN = 11 


class Michelogram: 
    def __init__(self, n_rings=64, span=DEFAULT_SPAN, max_ring_difference=60): 
        span = abs(int(span))
        max_ring_difference = abs(int(max_ring_difference))
        n_rings = abs(int(n_rings))
        if (span/2)*2 == span: 
            span = span+1
            print "Span must be odd, using span %d"%span
        if not self._verify_max_ring_difference(span, n_rings, max_ring_difference): 
            max_ring_difference_actual = self._largest_compatible_max_ring_difference(span, n_rings, max_ring_difference)
            compatible_max_ring_diff = self._compatible_max_ring_differences(span, n_rings)
            if compatible_max_ring_diff is []: 
                print "Michelogram: No 'maximum_ring_difference' exists compatible with given 'span' and 'n_rings'. "
            elif max_ring_difference_actual is None: 
                print "Michelogram: 'maximum_ring_difference' %d too small, using %d. "%(max_ring_difference,compatible_max_ring_diff[0])
                max_ring_difference = compatible_max_ring_diff[0]
            else: 
                print "Michelogram: 'maximum_ring_difference' %d is not compatible with 'span' %d. Compatible values are "%(max_ring_difference,span)+str(compatible_max_ring_diff)+". \nUsing %d"%max_ring_difference_actual
                max_ring_difference = max_ring_difference_actual
        self.n_rings = n_rings
        self.span = span 
        self.max_ring_difference = max_ring_difference 
        self.n_segments  = None 
        self.segments_sizes = None 
        self._svg_string = None 
        self._make()
        self._make_svg()

    def _make(self): 
        n_rings = self.n_rings
        n_rings_small = self.span/2
        n_rings_large = n_rings_small+1 
        if self.span>1: 
            n_segments_per_side = (self.max_ring_difference+1-self.span/2)/self.span 
        else: 
            n_segments_per_side = self.max_ring_difference
        n_segments = n_segments_per_side * 2 + 1
        segments_sizes = zeros([n_segments_per_side+1,],dtype=int32)
        if self.span>1: 
            segments_sizes[0] = 2*n_rings-1 
            if n_segments_per_side>0: 
                segments_sizes[1] = segments_sizes[0] - self.span - 1
                for i in range(n_segments_per_side-1): 
                    segments_sizes[2+i]=segments_sizes[1+i]-2*self.span
        else: 
            segments_sizes[0] = n_rings
            for i in range(n_segments_per_side): 
                segments_sizes[i+1]=segments_sizes[i]-1           
        self.n_segments = n_segments
        self.segments_sizes    = segments_sizes

        self.michelogram_sinogram = -1*ones([n_rings,n_rings],order = 'F',dtype=int32)
        self.michelogram_plane    = -1*ones([n_rings,n_rings],order = 'F',dtype=int32)

        if self.span > 1: 
            # -1- right sinograms  
            sinogram_index = int(n_segments)/2
            for X in range(0,self.max_ring_difference,self.span): 
                if sinogram_index == int(n_segments)/2: 
                    plane_index_offset = -n_rings_small
                else: 
                    plane_index_offset = 0
                self._fill_michelogram(X, -n_rings_small, sinogram_index, plane_index_offset)
                sinogram_index += 1

            # -2- left sinograms 
            sinogram_index = int(n_segments)/2-1
            for Y in range(n_rings_large,self.max_ring_difference,self.span): 
                plane_index_offset = 0
                self._fill_michelogram(0, Y, sinogram_index, plane_index_offset)
                sinogram_index -= 1
        else: 
            # -1- right sinograms 
            for seg in range((self.n_segments+1)/2): 
                for y in range(self.segments_sizes[seg]): 
                    x = y+seg
                    self.michelogram_sinogram[x,y] = seg + (self.n_segments+1)/2 - 1
                    self.michelogram_plane[x,y] = y
            # -2- left sinograms 
            for seg in range((self.n_segments+1)/2): 
                for x in range(self.segments_sizes[seg]): 
                    y = x+seg
                    self.michelogram_sinogram[x,y] = (self.n_segments+1)/2 - seg - 1
                    self.michelogram_plane[x,y] = x

    def _verify_max_ring_difference(self, span, n_rings, max_ring_difference): 
        return max_ring_difference in self._compatible_max_ring_differences(span, n_rings)

    def _compatible_max_ring_differences(self, span, n_rings): 
        c = [span/2]
        while True: 
            cnew = c[-1]+span
            if cnew < n_rings:
                c.append(cnew)
            else: 
                break
        return c 

    def _largest_compatible_max_ring_difference(self, span, n_rings, max_ring_difference):
        c = self._compatible_max_ring_differences(span, n_rings)
        for i in range(len(c)):
            if c[len(c)-1-i] < max_ring_difference: 
                return c[len(c)-1-i]
        return None
    
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

    def display_michelogram(self, cmap='gray', figsize=[16,8]): 
        import matplotlib.pyplot as plt
        plt.figure(figsize=figsize)
        plt.subplot(1,2,1); 
        plt.imshow(self.michelogram_plane,interpolation='nearest')
        plt.set_cmap(cmap)
        plt.subplot(1,2,2); 
        plt.imshow(self.michelogram_sinogram,interpolation='nearest')
        plt.set_cmap(cmap)

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





