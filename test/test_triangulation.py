import os
import unittest

from grid2tin.triangulation import Triangulation


class TestTriangulationRaster(unittest.TestCase):
    def setUp(self):
        self.path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/dgm5.tif')

    def test_limits(self):
        tri = Triangulation(self.path)
        repeat = True
        vertex_limit = 4000
        error_limit = 1.0
        vertex_count = None
        while repeat:
            error, vertex_count = tri.insert_next()
            if vertex_count >= vertex_limit or error <= error_limit:
                repeat = False
        self.assertLessEqual(vertex_count, vertex_limit)
