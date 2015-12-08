

class BaseScintillatorSPECT(): 
    def __init__(self): 
        self._name = "Generic Scintillator"
        self._manufacturer = "No manufacturer"
        self._scintillator_material = "Not specified"
        self._p_thickness_mm = 10.0
        self._need_to_recompute_psf = True
        self._previous_energy_request = 0

    def get_name(self): 
        return self._name

    def get_manufacturer(self): 
        return self._manufacturer

    def get_material(self): 
        return self._scintillator_material

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
        print_medium_verbose("Recomputing PSF..")
        self._need_to_recompute_psf = False
        self._psf = 0


class Ideal(BaseScintillatorSPECT): 
    def __init__(self): 
        BaseScintillatorSPECT.__init__(self)
        self._name = "Ideal Scintillator"

    def get_psf(self,energy_kev): 
        if self._need_to_recompute_psf or energy_kev!=self._previous_energy_request: 
            self._compute_psf(energy_kev)
            self._previous_energy_request = energy_kev
        return self._psf 

    def _compute_psf(self,energy_kev): 
        print_medium_verbose("Recomputing PSF scintillator..")
        self._need_to_recompute_psf = False
        self._psf = 0




