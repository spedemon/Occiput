

class BaseCollimatorSPECT(): 
    def __init__(self): 
        self._name = "Generic Collimator" 
        self._collimator_type = "Generic Collimator" 
        self._manufacturer = "No manufacturer"
        self._p_simulation_method    = "analytical" 
        self._p_hole_diameter_mm     = 1.0 
        self._p_thickness_mm         = 40.0 
        self._p_gap_scintillatpr_mm  = 0.0 
        self._p_septum_mm            = 0.3 
        self._p_hole_shape           = "hexagonal" 
        self._need_to_recompute_psf = True 
        self._previous_energy_request = 0
        self._psf = None

    def get_hole_diameter(self):
        return self._p_hole_diameter_mm

    def set_hole_diameter(self,diameter_mm):
        if not np.isscalar(diameter_mm): 
            raise BadParameter('Expected a scalar value.')
        self._p_hole_diameter = diameter_mm
        self._need_to_recompute_psf = True

    def get_thickness(self):
        return self._p_thickness_mm

    def set_thickness(self,thickness_mm):
        if not np.isscalar(thickness_mm): 
            raise BadParameter('Expected a scalar value.')
        self._p_thickness_mm = thickness_mm
        self._need_to_recompute_psf = True

    def get_gap_scintillator(self):
        return self._p_gap_scintillator_mm

    def set_gap_scintillator(self,gap_scintillator_mm):
        if not np.isscalar(gap_scintillator_mm): 
            raise BadParameter('Expected a scalar value.')
        self._p_gap_scintillator_mm = gap_scintillator_mm
        self._need_to_recompute_psf = True

    def get_septum(self):
        return self._p_septum_mm

    def set_septum(self,septum_mm):
        if not np.isscalar(septum_mm): 
            raise BadParameter('Expected a scalar value.')
        self._p_septum_mm = septum_mm
        self._need_to_recompute_psf = True

    def get_hole_shape(self):
        return self._p_hole_shape

    def set_hole_shape(self,hole_shape):
        if str(hole_shape) == "hexagonal" or str(hole_shape) == "circular": 
            self._p_hole_shape = str(hole_shape)
        else: 
            raise BadParameter('Unknown hole shape') 
        self._need_to_recompute_psf = True

    def get_simulation_method(self): 
        return self._p_simulation_method 

    def set_simulation_method(self,simulation_method): 
        if str(simulation_method) == "analytical": 
            self._p_simulation_method = str(simulation_method)
        else: 
            raise BadParameter('Unknown simulation method') 
        self._need_to_recompute_psf = True

    def get_parameters(self):
        parameters = {}
        dic = self.__dict__
        for k in dic.keys(): 
            if k.startswith('_p_'):
                parameters[k[3:]]=dic[k]        
        return parameters 

    def get_psf(self,energy_kev): 
        if self._need_to_recompute_psf or energy_kev!=self._previous_energy_request: 
            self._compute_psf(energy_kev)
            self._previous_energy_request = energy_kev
        return self._psf 

    def _compute_psf(self,energy_kev): 
        print_medium_verbose("Recomputing PSF collimator ..")
        self._need_to_recompute_psf = False
        self._psf = 0



class LEHR(BaseCollimatorSPECT):
    def __init__(self): 
        BaseCollimatorSPECT.__init__(self)
        self._name = "Generic LEHR Collimator"
        self._collimator_type = "Low Energy High Resolution - LEHR"
        self._manufacturer = "No manufacturer" 

    def get_psf(self,energy_kev): 
        if self._need_to_recompute_psf or energy_kev!=self._previous_energy_request: 
            self._compute_psf(energy_kev)
            self._previous_energy_request = energy_kev
        return self._psf 

    def _compute_psf(self,energy_kev): 
        print_medium_verbose("Recomputing PSF collimator ..")
        self._need_to_recompute_psf = False
        self._psf = 0


class HELR(BaseCollimatorSPECT):
    def __init__(self): 
        BaseCollimatorSPECT.__init__(self)
        self._name = "Generic HELR Collimator"
        self._collimator_type = "High Energy Low Resolution - HELR"
        self._manufacturer = "No manufacturer" 

    def get_psf(self,energy_kev): 
        if self._need_to_recompute_psf or energy_kev!=self._previous_energy_request: 
            self._compute_psf(energy_kev)
            self._previous_energy_request = energy_kev
        return self._psf 

    def _compute_psf(self,energy_kev): 
        print 1
        print_low_verbose("Recomputing PSF collimator ..")
        self._need_to_recompute_psf = False
        self._psf = 0



