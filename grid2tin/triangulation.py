# Implementation of Garland and Heckbert's sequential greedy insertion
# algorithm for terrain approximation: http://mgarland.org/software/terra.html


import logging
from math import ceil

import numpy as np
import pyximport
import rasterio
from affine import Affine

from .heap import Heap
from .quadedge import Vertex, splice, connect, swap, make_edge, Triangle

pyximport.install()
# noinspection PyPep8
from .calculation import scan_triangle_line

logging.basicConfig(level=logging.WARN)


class Triangulation:
    def __init__(self, dem, minimum_gap=5):
        if isinstance(dem, np.ndarray):
            self.dem = dem
            self.affine = None
        elif isinstance(dem, str):
            with rasterio.Env():
                with rasterio.open(dem) as src:
                    rawdata = src.read()
                    self.affine = src.affine
                    self.dem = np.array(rawdata.squeeze(), dtype=float)

        self.minimum_gap = minimum_gap

        min_x = 0
        min_y = 0
        max_x = self.dem.shape[1] - 1
        max_y = self.dem.shape[0] - 1

        self.heap = Heap()
        self.triangle_list = []

        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y

        self.available = np.ones_like(self.dem)

        self.vertex_dict = dict()
        self.edge_dict = dict()

        v0 = Vertex(min_x, min_y, self.dem[0, 0])
        v1 = Vertex(max_x, min_y, self.dem[0, -1])
        v2 = Vertex(max_x, max_y, self.dem[-1, -1])
        v3 = Vertex(min_x, max_y, self.dem[-1, 0])

        self.vertex_dict[0] = v0
        self.vertex_dict[1] = v1
        self.vertex_dict[2] = v2
        self.vertex_dict[3] = v3
        self.next_vertex_id = 4

        self.next_edge_id = 0
        # Boundary rectangle
        q0 = make_edge(v0, v1)
        q1 = make_edge(v2, v3)
        q2 = make_edge(v3, v0)
        q3 = make_edge(v1, v2)
        # Diagonal
        q4 = make_edge(v1, v3)

        splice(q0.sym, q4)
        splice(q4.sym, q2)
        splice(q2.sym, q0)
        splice(q0.sym, q3)
        splice(q3.sym, q1)
        splice(q1.sym, q4.sym)

        self.add_edge(q0)
        self.add_edge(q1)
        self.add_edge(q2)
        self.add_edge(q3)

        # Mark area around border edges as unavailable
        for e in self.edges:
            self.mark_availability(e.origin,
                                   e.destination,
                                   radius=self.minimum_gap,
                                   value=0)
        # Mark area on border edges themselves as available
        for e in self.edges:
            self.mark_availability(e.origin,
                                   e.destination,
                                   radius=0,
                                   value=1)
        # Mark area around border vertices as unavailable
        for v in self.vertices:
            self.mark_availability(v, radius=minimum_gap, value=0)

        self.add_edge(q4)

        self.base = q0

        self.history = Triangle(self.base, anchor=False, id_=-1)
        self.history.children = [Triangle(q4),
                                 Triangle(q4.sym)]

        for triangle in self.history.children:
            self.triangle_list.append(triangle)
            self.scan_triangle(triangle)
            triangle.id = self.heap.insert(triangle.candidate_error,
                                           (triangle.candidate, triangle))

    @property
    def vertices(self):
        return [self.vertex_dict[key]
                for key in self.vertex_dict]

    @property
    def edges(self):
        return [self.edge_dict[key]
                for key in self.edge_dict]

    @property
    def triangles(self):
        return list(filter(lambda x: x.id != -1, self.triangle_list))

    # Point location using the history graph. Much faster than the walking
    # method
    def search(self, v):
        triangle = None
        current_triangle = self.history

        while len(current_triangle.children) > 0:
            for triangle in current_triangle.children:
                if v.in_triangle(triangle.vertices[0],
                                 triangle.vertices[1],
                                 triangle.vertices[2]):
                    current_triangle = triangle
                    break
            if not current_triangle == triangle:
                # None of the children contained the point, point is not in
                # triangulation
                logging.debug("Point {} not in triangulation, no edge fount.".format(v))
                return None
        return current_triangle.anchor

    def add_edge(self, e):
        self.edge_dict[self.next_edge_id] = e
        e.id = e.sym.id = self.next_edge_id
        self.next_edge_id += 1

    def delete_edge(self, e):
        splice(e, e.o_prev)
        splice(e.sym, e.sym.o_prev)
        del self.edge_dict[e.id]
        del e

    def insert_site(self, v, e=None):
        """
        Insert a new site in the 2D triangulation, maintaining the Delaunay
        criterion
        :param v: Vertex to be inserted into the triangulation
        :param e: Optional: edge of the triangle that contains v, to speed up
        the process
        :return:
        """
        deleted_triangles = []
        created_triangles = []
        boundary_edge = None

        # Get elevation from map if none was provided in vertex
        if v.z == 0:
            v.z = self.dem[v.y, v.x]

        if e is None:
            e = self.search(v)
        else:
            assert (not v.right_of(e) and
                    not v.right_of(e.l_next) and
                    not v.right_of(e.l_prev)), \
                'Edge %s is not an edge of the ' \
                'triangle containing edge %s' % (e, v)
        if v == e.origin or v == e.destination:
            return deleted_triangles, created_triangles
        elif v.on_edge(e):
            if not e.o_prev.destination.right_of(e):
                parents = [e.triangle]
                boundary_edge = e
            else:
                parents = [e.triangle, e.sym.triangle]
                e = e.o_prev
                self.delete_edge(e.o_next)
        else:
            parents = [e.triangle]

        # Add point to triangulation
        self.vertex_dict[self.next_vertex_id] = v
        self.next_vertex_id += 1

        # Create first spoke from origin of base to new site
        spoke = make_edge(e.origin, v)
        self.add_edge(spoke)

        splice(spoke, e)
        starting_spoke = spoke

        # Create second spoke from destination of base to new site
        spoke = connect(e, spoke.sym)
        self.add_edge(spoke)

        e = spoke.o_prev
        while e.l_next != starting_spoke:
            spoke = connect(e, spoke.sym)
            self.add_edge(spoke)
            e = spoke.o_prev

        if boundary_edge is not None:
            # We might be deleting the base edge of the mesh, so use an edge
            # that is guaranteed to survive
            self.base = e
            self.delete_edge(boundary_edge)

        # Create (potentially ephemeral) triangles for all the spokes and assign
        # them as children to the parent triangles
        current_spoke = starting_spoke
        while True:
            current_spoke = current_spoke.d_next
            if current_spoke.o_next.destination.left_of(current_spoke):
                child = Triangle(current_spoke)
                created_triangles.append(child)
                for parent in parents:
                    parent.children.append(child)

            if current_spoke == starting_spoke:
                break

        for parent in parents:
            parent.anchor = None
        deleted_triangles.extend(parents)

        while True:
            t = e.o_prev
            if t.destination.right_of(e) and v.in_circle(e.origin,
                                                         t.destination,
                                                         e.destination):
                # Delaunay criterion is violated, swap an edge to fix it
                # This deletes two triangles and creates two new ones
                parents = [e.triangle, e.sym.triangle]
                swap(e)
                deleted_triangles.extend(parents)

                children = [Triangle(e),
                            Triangle(e.sym)]
                created_triangles.extend(children)
                for parent in parents:
                    parent.children.extend(children)
                    parent.anchor = None

                e = e.o_prev
            elif e.o_next == starting_spoke:
                break
            else:
                e = e.o_next.l_prev

        created_triangles = list(set(created_triangles) -
                                 set(deleted_triangles))

        self.mark_availability(v, radius=self.minimum_gap, value=0)
        return created_triangles, deleted_triangles

    def scan_triangle(self, t, interpolation_map=None, only_return_points=False):
        v0, v1, v2 = t.vertices

        # Sort vertices in ascending order
        if v0.y > v1.y:
            v0, v1 = v1, v0
        if v0.y > v2.y:
            v0, v2 = v2, v0
        if v1.y > v2.y:
            v1, v2 = v2, v1

        # Check if base of triangle is flat
        if v1.y == v0.y:
            dx0 = 0.0
        else:
            dx0 = (v1.x - v0.x) / (v1.y - v0.y)

        dx1 = (v2.x - v0.x) / (v2.y - v0.y)

        y_start = v0.y
        y_end = v1.y
        x_a = v0.x
        x_b = v0.x

        points = []
        # If the base of the triangle is flat, this loop won't be executed
        for y in range(y_start, y_end):
            points.extend(scan_triangle_line(self.available, self.dem, t, y,
                                             x_a, x_b,
                                             interpolation_map,
                                             only_return_points))
            x_a += dx0
            x_b += dx1

        # We have reached the second half of the triangle
        # Check if top of triangle is flat
        if v2.y == v1.y:
            dx0 = 0.0
        else:
            dx0 = (v2.x - v1.x) / (v2.y - v1.y)

        y_start = v1.y
        y_end = v2.y
        x_a = v1.x

        # If the top of the triangle is flat, this loop will be executed once
        for y in range(y_start, y_end + 1):
            points.extend(scan_triangle_line(self.available, self.dem, t, y,
                                             x_a, x_b,
                                             interpolation_map,
                                             only_return_points))
            x_a += dx0
            x_b += dx1

        if only_return_points:
            return points

    def circle_points(self, center, radius):
        circle_points = []
        y_start = (max(round(center.y - radius), self.min_y))
        y_end = (min(round(radius + 1 + center.y), self.max_y + 1))

        for y in range(y_start, y_end):
            x_max = max(radius ** 2 - (y - center.y) ** 2, 0) ** 0.5
            x_start = int(max(round(center.x - x_max), self.min_x))
            x_end = int(min(round(x_max + 1 + center.x), self.max_x))
            for x in range(x_start, x_end + 1):
                circle_points.append((x, y))

        return circle_points

    @staticmethod
    def segment_points(s0, s1):
        """
        Find all points along the selection segment of an edge
        :param s0: starting vertex of segment
        :param s1: ending vertex of segment
        :return: list of 2D coordinates along segment
        """
        segment_points = []
        a = s1 - s0
        d = ceil(a.norm)

        step = 1 / d

        for i in range(d):
            v = s0 + i * step * a
            # noinspection PyUnresolvedReferences
            x = int(round(v.x))
            # noinspection PyUnresolvedReferences
            y = int(round(v.y))
            segment_points.append((x, y))
        return segment_points

    def mark_availability(self, v0, v1=None, radius=0, value=0):
        if v1 is not None:
            segment_points = self.segment_points(v0, v1)
        else:
            segment_points = [(v0.x, v0.y)]

        for s in segment_points:
            cp = self.circle_points(Vertex(s[0], s[1]), radius)
            for p in cp:
                self.available[p[1], p[0]] = value

    def insert_point(self, v, e=None):
        """
        Insert a new vertex into the triangulation and scan the newly created
        triangles for the error
        :param v: Vertex to be inserted
        :param e: Optional edge for starting the triangle search
        :return:
        """
        new, deleted = self.insert_site(v, e)

        for triangle in deleted:
            if not triangle.id == -1:
                self.heap.delete(triangle.id)
                triangle.id = -1

        for triangle in new:
            triangle.id = -2

        for triangle in new:
            self.scan_triangle(triangle)
            triangle.id = self.heap.insert(triangle.candidate_error,
                                           (triangle.candidate, triangle))
        self.mark_availability(v, radius=self.minimum_gap, value=0)
        self.triangle_list.extend(new)

    def insert_next(self):
        """
        Pop the candidate with the greatest error and insert it into the
        triangulation
        :return:
        """
        error, (candidate, triangle) = self.heap.pop()
        triangle.id = -1  # Mark it as removed from the heap

        new, deleted = self.insert_site(candidate)

        for triangle in deleted:
            if not triangle.id == -1:
                self.heap.delete(triangle.id)
                triangle.id = -1

        for triangle in new:
            self.scan_triangle(triangle)
            triangle.id = self.heap.insert(triangle.candidate_error,
                                           (triangle.candidate, triangle))
        self.triangle_list.extend(new)
        return error, len(self.vertex_dict)

    def interpolated_map(self):
        """
        The height map resulting from linear interpolation of the triangle
        mesh
        :return:
        """
        interpolated_map = self.dem.copy()
        for triangle in self.triangles:
            self.scan_triangle(triangle, interpolated_map)
        return interpolated_map

    def error_map(self):
        """
        The difference map between the original height map and the interpolated
        height map
        :return:
        """
        error_map = self.dem.copy() - self.interpolated_map()
        return error_map

    def write_obj(self, filename):
        if not self.affine:
            self.affine = Affine.identity()
        coordinates = np.array([self.affine * v.pos[:2] + tuple([v.pos[2]]) for v in self.vertices])
        v_fmt = "v %.3f %.3f %.3f"
        triangles = [[self.vertices.index(vertex) for vertex in triangle.vertices][::-1] for triangle in self.triangles]

        texture_coordinates = coordinates[:, :2].copy()
        texture_coordinates -= texture_coordinates.min(axis=0)
        texture_coordinates /= texture_coordinates.ptp(axis=0)

        with open(filename, 'wb') as outfile:
            np.savetxt(outfile, coordinates, fmt=v_fmt)
            np.savetxt(outfile, texture_coordinates, fmt="vt %.3f %.3f")
            np.savetxt(outfile,
                       np.dstack([triangles, triangles]).reshape(-1, 6) + 1,
                       fmt="f %i/%i/ %i/%i/ %i/%i/")
