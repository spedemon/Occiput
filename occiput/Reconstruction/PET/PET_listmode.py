
# occiput 
# Stefano Pedemonte
# Aalto University, School of Science, Helsinki
# Oct 2013, Helsinki 
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Dec. 2013, Boston
# Martinos Center for Biomedical Imaging, Harvard University/MGH, Boston
# Jan. 2015, Boston
# Feb. 2015, Helsinki
# Nov. 2015, Boston 


# Interface to list-mode data. In order to cope with various data formats of the list-mode data and 
# to enable fast decoding, occiput uses a simple plug-in system. The routines 
# to decode the list-mode data are defined in a C function that is loaded as a plugin. 
# The simple plugin system handles the compilation of the code. 



from occiput.Reconstruction.PET.PET_Projection import PET_Projection, PET_Projection_Sparsity
from occiput.Core.NiftyPy_wrap import PET_listmode_get_measurement_static, PET_listmode_get_measurement

__all__ = ["Listmode_Loader"]


class ListmodeWrapper(): 
    def open_stream(self, stream): 
        PET_listmode_open_stream(stream)
    
    def close_stream(self): 
        PET_listmode_close_stream()
    
    def load_scanner_plugin(self, plugin): 
        scanner_library = plugin.get_library_path()
        PET_listmode_load_scanner_plugin(scanner_library) 

    def load_stream_plugin(self, plugin): 
        stream_library = plugin.get_library_path()
        PET_listmode_load_stream_plugin(stream_library)     
    
    def bin_events(self, n_events, Naxial, Nazim, Nu, Nv, time_bins): 
        PET_listmode_bin(n_events, Naxial, Nazim, Nu, Nv, time_bins)

    def bin_events_intrinsic(self, n_events): 
        PET_listmode_bin_intrinsic(n_events)

    def get_static_measurement(self): 
        R = PET_listmode_get_static()
        return R['counts'], R['offsets'], R['locations'], R['time_start'], R['time_end'] 

    def get_dynamic_measurement(self): 
        R = PET_listmode_get_dynamic()
        return R['counts'], R['offsets'], R['locations'], R['time_start'], R['time_end'] 



class ListmodeLoader(): 
    def __init__(self): 
        self._LM = ListmodeWrapper()
    
    def set_streamer(self, plugin): 
        # call C function that loads dynamically the open(), close() and read() functions from the streamer dynamic library
        self._LM.load_stream_plugin( plugin )
    
    def set_scanner(self, plugin): 
        # call C function that loads dynamically the index_to_lor() function of the scanner plugin 
        self._LM.load_scanner( plugin ) 
    
    def open_stream(self, stream): 
        # call the C function that opens a list-mode stream 
        self._LM.open_stream( stream )
    
    def close_stream(self): 
        # call the C function that closes the list-mode stream 
        self._LM.close_stream()

    def make_projection_given_binning(self, binning, n_events=-1, time_bins=None):  
        # call the C function that parses the list-mode data and bins it into (dynamic) projection data 
        self._LM.bin_events(n_events, binning.N_axial, binning.N_azimuthal, binning.N_u, binning.N_v, time_bins)
        counts,offsets,locations,time_start,time_end = self._LM.get_static_measurement()
        time_bins = int32(linspace(time_start, time_end, 2))
        return PET_Projection( binning, counts, offsets, locations, time_bins) 

    def make_projection_intrinsic_geometry(self, n_events=-1, time_bins=None, N_azimuthal=-1): 
        # call the C function that parses the list-mode data and bins it into (dynamic) projection data 
        self._LM.bin_events_intrinsic(n_events, binning.N_axial, binning.N_azimuthal, binning.N_u, binning.N_v, time_bins)
        counts,offsets,locations,time_start,time_end = self._LM.get_static_measurement()
        time_bins = int32(linspace(time_start, time_end, 2))
        return PET_Projection( binning, counts, offsets, locations, time_bins) 


