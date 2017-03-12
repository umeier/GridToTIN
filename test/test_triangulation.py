import os
import unittest

import numpy as np

from grid2tin.quadedge import Vertex
from grid2tin.triangulation import Triangulation


class TestTriangulationRaster(unittest.TestCase):
    def setUp(self):
        self.path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/dgm5.tif')
        self.grid = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], np.int32)

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
        mx, my = np.meshgrid(x, y)
        z1 = np.exp(-2 * np.log(2) * ((mx - 0.5) ** 2 + (my - 0.5) ** 2) / 1 ** 2)
        z2 = np.exp(-3 * np.log(2) * ((mx + 0.5) ** 2 + (my + 0.5) ** 2) / 2.5 ** 2)
        z = 10.0 * (z2 - z1)
        tri = Triangulation(z)
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

    def test_insert_point_out_of_grid(self):
        tri = Triangulation(self.path, minimum_gap=0)
        with self.assertRaises(IndexError):
            tri.insert_point(Vertex(tri.max_y + 1, tri.max_y + 1))

    def test_write_obj_minimal_grid(self):
        tri = Triangulation(self.grid, minimum_gap=0)
        tri.write_obj('grid.obj')

    def test_write_obj(self):
        tri = Triangulation(self.path, minimum_gap=0)
        self.do_triangulation(tri, limit=4000)
        tri.write_obj('dgm5.obj')

    def do_triangulation(self, tri, limit=100):
        repeat = True
        vertex_limit = limit
        error_limit = 5.0
        vertex_count = None
        while repeat:
            error, vertex_count = tri.insert_next()
            if vertex_count >= vertex_limit or error <= error_limit:
                repeat = False
        self.assertLessEqual(vertex_count, vertex_limit)
