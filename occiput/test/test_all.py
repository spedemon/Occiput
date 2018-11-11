# -*- coding: utf-8 -*-
# occiput
# Harvard University, Martinos Center for Biomedical Imaging
# Aalto University, Department of Computer Science

from .. import occiput
import unittest


class Test_Reconstruction_PET(unittest.TestCase):
    """Sequence of tests for tomographic reconstruction - Positron Emission Tomography."""

    def setUp(self):
        pass

    def test_projection_wrapper(self):
        """Test the Python wrapper for the projection algorithm. """
        number = 0.1
        descriptor = [
            {"name": "input", "type": "int", "value": number},
            {"name": "output", "type": "int", "value": None},
        ]
        r = c_python.call_c_function(self.lib.echo, descriptor)
        self.assertTrue(r.output == number)


if __name__ == "__main__":
    unittest.main()
