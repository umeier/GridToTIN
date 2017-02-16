import unittest

from grid2tin.quadedge import Vertex, QuadEdge


class TestVertex(unittest.TestCase):
    def setUp(self):
        self.vex_a = Vertex(2, 1)
        self.vex_b = Vertex(6, 5)
        self.vex_c = Vertex(2, 3)
        self.vex_d = Vertex(5, 2)
        self.vex_e = Vertex(4, 3)
        self.vex_f = Vertex(4, 7)
        self.vex_g = Vertex(3, 3)
        self.vex_h = Vertex(2, 6)
        self.vex_i = Vertex(3, 10)
        self.vex_j = Vertex(7, 7)
        self.q = QuadEdge(self.vex_a, self.vex_b)
        self.e = self.q.base

    def test_vertex(self):
        self.assertEqual(self.vex_a.x, 2)
        self.assertEqual(self.vex_a.y, 1)
        self.assertEqual(self.vex_a.z, 0)

    def test_on_edge(self):
        self.assertTrue(self.vex_e.on_edge(self.e))
        self.assertFalse(self.vex_c.on_edge(self.e))

    def test_left_right(self):
        self.assertTrue(self.vex_c.left_of(self.e))
        self.assertFalse(self.vex_d.left_of(self.e))
        self.assertTrue(self.vex_d.right_of(self.e))
        self.assertFalse(self.vex_c.right_of(self.e))

    def test_in_triangle(self):
        self.assertTrue(self.vex_g.in_triangle(self.vex_a, self.vex_b, self.vex_c))
        self.assertFalse(self.vex_f.in_triangle(self.vex_a, self.vex_b, self.vex_c))

    def test_in_circle(self):
        self.assertTrue(self.vex_f.in_circle(self.vex_h, self.vex_j, self.vex_i))
        self.assertFalse(self.vex_g.in_circle(self.vex_h, self.vex_j, self.vex_i))

    def test_encroaches(self):
        self.assertTrue(self.vex_c.encroaches(self.e))
        self.assertFalse(self.vex_h.encroaches(self.e))
