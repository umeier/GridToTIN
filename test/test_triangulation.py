import os
import unittest

from triangulation import Triangulation


class TestTriangulation(unittest.TestCase):
    def test_raster(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/dgm1.img')
        tri = Triangulation(path)
        tri.insert_next()