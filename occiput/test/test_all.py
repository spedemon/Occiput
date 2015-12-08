
# occiput
# Stefano Pedemonte
# Harvard University / MGH
# Boston, Jan 2015


from .. import occiput

import unittest


class Test_Reconstruction_PET(unittest.TestCase): 
    """Sequence of tests for tomographic reconstruction - Positron Emission Tomography """   
    def setUp(self):
        pass 

    def test_projection_wrapper(self):
        """Test the Python wrapper for the projection algorithm. """
        number = 0.1
        descriptor = [  {'name':'input',  'type':'int', 'value':number},
                        {'name':'output', 'type':'int', 'value':None },  ]
        r = c_python.call_c_function( self.lib.echo, descriptor ) 
        self.assertTrue(r.output == number)

    def test_projection(self):
        """Test the projection algorithm. """
        A = 0
        B = 0
        self.assertTrue(A == B)


class Test_Reconstruction_SPECT(unittest.TestCase): 
    def setUp(self):
        pass 
        
class Test_Reconstruction_CT(unittest.TestCase): 
    def setUp(self): 
        pass 



if __name__ == '__main__':
    unittest.main()
    

