import os
import unittest

import numpy as np
import rasterio

from grid2tin.quadedge import Vertex
from grid2tin.triangulation import Triangulation


class TestTriangulationRaster(unittest.TestCase):
    def setUp(self):
        self.path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/dgm5.tif')

    def test_triag_file_gap(self):
        tri = Triangulation(self.path)
        self.do_triangulation(tri)

    def test_triag_file_no_gap(self):
        tri = Triangulation(self.path, minimum_gap=0)
        self.do_triangulation(tri)

    def test_triag_array_gap(self):
        """
        Grid example from rasterio doku
        """
        x = np.linspace(-4.0, 4.0, 240)
        y = np.linspace(-3.0, 3.0, 180)
        X, Y = np.meshgrid(x, y)
        Z1 = np.exp(-2 * np.log(2) * ((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 1 ** 2)
        Z2 = np.exp(-3 * np.log(2) * ((X + 0.5) ** 2 + (Y + 0.5) ** 2) / 2.5 ** 2)
        Z = 10.0 * (Z2 - Z1)
        tri = Triangulation(Z)
        self.do_triangulation(tri)

    def test_triag_force_points_edges(self):
        tri = Triangulation(self.path, minimum_gap=0)
        x_vals = np.linspace(tri.min_x, tri.max_x, 10)
        y_vals = np.linspace(tri.min_y, tri.max_y, 10)
        for x in x_vals:
            tri.insert_point(Vertex(int(x), tri.min_y))
            tri.insert_point(Vertex(int(x), tri.max_y))
        for y in y_vals:
            tri.insert_point(Vertex(tri.min_x, int(y)))
            tri.insert_point(Vertex(tri.max_x, int(y)))
        self.do_triangulation(tri)

    def do_triangulation(self, tri):
        repeat = True
        vertex_limit = 100
        error_limit = 5.0
        vertex_count = None
        while repeat:
            error, vertex_count = tri.insert_next()
            if vertex_count >= vertex_limit or error <= error_limit:
                repeat = False
        self.assertLessEqual(vertex_count, vertex_limit)
