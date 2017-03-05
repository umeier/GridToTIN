from math import ceil, floor

def calc_interpolation(double a, double b, double c, int x, int y):
    # for call from external
    return a * x + b * y + c

cdef double interpolation(double a, double b, double c, int x, int y):
    # for call from internal
    return a * x + b * y + c

def scan_triangle_line(available, dem, t, int y, x_a, x_b, interpolation_map=None,
                       only_return_points=False):
    """
    This is the most time consuming part of the triangulation.
    """

    cdef int x_start
    cdef int x_end
    cdef double a
    cdef double b
    cdef double c
    cdef int x

    x_start = int(ceil(min(x_a, x_b)))
    x_end = int(floor(max(x_a, x_b)))

    a = t.a
    b = t.b
    c = t.c

    points = []
    for x in range(x_start, x_end + 1):
        if only_return_points:
            points.append((x, y))
        else:
            z_map = dem[y, x]
            error = abs(z_map - interpolation(a, b, c, x, y))
            if interpolation_map is not None:
                interpolation_map[y, x] = interpolation(a, b, c, x, y)
            if error > t.candidate_error and available[y, x] == 1:
                t.candidate_error = error
                t.candidate.pos = (x, y, z_map)
    return points
