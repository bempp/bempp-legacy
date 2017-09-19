"""Implement various fractal shapes."""

from bempp.api.shapes.shapes import __generate_grid_from_geo_string

# 3D fractals

#pylint: disable=invalid-name
#pylint: disable=too-many-locals
#pylint: disable=too-many-statements
#pylint: disable=too-many-nested-blocks
#pylint: disable=too-many-branches

def sierpinski_pyramid(h=0.1, level=2):
    """
    Returns a Sierpinski pyramid.

    Parameters
    ----------
    h : float
        Size of the elements.
    level : int
        Refinement level.

    """
    from numpy import sqrt
    from numpy import array
    if level < 1:
        raise ValueError("level must be 1 or larger")
    geo = "lc = " + str(h) + ";\n"
    points = [array((.5, 1 / (2 * sqrt(3)), sqrt(2. / 3.))),
              array((0, 0, 0)), array((1, 0, 0)), array((.5, sqrt(3) / 2, 0))]
    tetras = [[0, 1, 2, 3]]
    for _ in range(level):
        new_tetras = []
        vertex_map = {}
        for t in tetras:
            for i in t:
                new_t = []
                for j in t:
                    if i == j:
                        new_t.append(i)
                    else:
                        a = min(i, j)
                        b = max(i, j)
                        if b not in vertex_map or a not in vertex_map[b]:
                            if b not in vertex_map:
                                vertex_map[b] = {}
                            vertex_map[b][a] = len(points)
                            points.append((points[a] + points[b]) / 2.)
                        new_t.append(vertex_map[b][a])
                new_tetras.append(new_t)
        tetras = new_tetras

    for i, p in enumerate(points):
        geo += "Point(" + str(1 + i) + \
            ") = {" + str(p[0]) + "," + str(p[1]) + "," + str(p[2]) + ",lc};\n"
    for n, t in enumerate(tetras):
        geo += "Line(" + str(1 + 6 * n) + \
            ") = {" + str(1 + t[0]) + "," + str(1 + t[1]) + "};\n"
        geo += "Line(" + str(1 + 6 * n + 1) + \
            ") = {" + str(1 + t[0]) + "," + str(1 + t[2]) + "};\n"
        geo += "Line(" + str(1 + 6 * n + 2) + \
            ") = {" + str(1 + t[0]) + "," + str(1 + t[3]) + "};\n"
        geo += "Line(" + str(1 + 6 * n + 3) + \
            ") = {" + str(1 + t[1]) + "," + str(1 + t[2]) + "};\n"
        geo += "Line(" + str(1 + 6 * n + 4) + \
            ") = {" + str(1 + t[2]) + "," + str(1 + t[3]) + "};\n"
        geo += "Line(" + str(1 + 6 * n + 5) + \
            ") = {" + str(1 + t[3]) + "," + str(1 + t[1]) + "};\n"
        geo += ("Line Loop(" + str(1 + 4 * n) + ") = { " + str(1 + 6 * n) +
                ", " + str(1 + 6 * n + 3) + ",-" + str(1 + 6 * n + 1) + "};\n")
        geo += ("Line Loop(" + str(1 + 4 * n + 1) + ") = { " + str(
            1 + 6 * n + 1) + ", " + str(1 + 6 * n + 4) + ",-" +
                str(1 + 6 * n + 2) + "};\n")
        geo += ("Line Loop(" + str(1 + 4 * n + 2) + ") = { " + str(
            1 + 6 * n + 2) + ", " + str(1 + 6 * n + 5) + ",-" +
                str(1 + 6 * n) + "};\n")
        geo += ("Line Loop(" + str(1 + 4 * n + 3) + ") = {-" + str(
            1 + 6 * n + 4) + ",-" + str(1 + 6 * n + 3) + ",-" +
                str(1 + 6 * n + 5) + "};\n")
        for i in range(4):
            geo += "Plane Surface(" + str(1 + 4 * n + i) + \
                ") = {" + str(1 + 4 * n + i) + "};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)


def menger_sponge(h=0.1, level=2):
    """
    Returns a Menger sponge.

    Parameters
    ----------
    h : float
        Size of the elements.
    level : int
        Refinement level.

    """
    from numpy import array
    from itertools import product
    geo = "lc = " + str(h) + ";\n"

    N = 3**level
    def get_point(x, y, z):
        return x*(N+1)**2 + y*(N+1) + z

    def coords(x, y, z):
        out = ""
        for a in (x,y,z):
            out += str(a*1./N) + ","
        out += "lc"
        return out

    def is_solid(x, y, z):
        if x < 0 or y < 0 or z < 0:
            return False
        if x >= N or y >= N or z >= N:
            return False
        for i in range(level):
            for a,b in [(x,y),(x,z),(y,z)]:
                pow = 3**i
                if (a//pow)%3==1 and (b//pow)%3==1:
                    return False
        return True


    class EdgeGenerator(object):
        def __init__(self):
            self.geo_edges = ""
            self.edges = {}
            self.edge_n = 1

        def edge(self, x1,y1,z1, x2,y2,z2):
            start = get_point(x1,y1,z1)
            end = get_point(x2,y2,z2)
            if start not in self.edges:
                self.edges[start] = {}
            if end not in self.edges:
                self.edges[end] = {}
            if end not in self.edges[start]:
                self.geo_edges += "Line("+str(self.edge_n)+")="
                self.geo_edges += "{"+str(start)+","+str(end)+"};\n"
                self.edges[start][end] = str(self.edge_n)
                self.edges[end][start] = str(-self.edge_n)
                self.edge_n += 1
            return self.edges[start][end]

    geo_faces = ""
    geo_surfaces = ""
    loop_n = 1

    eg = EdgeGenerator()

    for x in range(-1,N+1):
        out = ""
        for y in range(-1,N+1):
            if is_solid(x,y,0):
                out += "X"
            else:
                out += "-"
        print(out)

    for x,y,z in product(range(N), repeat=3):
        if is_solid(x, y, z):
            if not is_solid(x-1, y, z):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x,y,z,     x,y+1,z)+","
                geo_faces += eg.edge(x,y+1,z,   x,y+1,z+1)+","
                geo_faces += eg.edge(x,y+1,z+1, x,y,z+1)+","
                geo_faces += eg.edge(x,y,z+1,   x,y,z)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1
            if not is_solid(x+1, y, z):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x+1,y,z,     x+1,y,z+1)+","
                geo_faces += eg.edge(x+1,y,z+1,   x+1,y+1,z+1)+","
                geo_faces += eg.edge(x+1,y+1,z+1, x+1,y+1,z)+","
                geo_faces += eg.edge(x+1,y+1,z,   x+1,y,z)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1
            if not is_solid(x, y-1, z):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x,y,z,     x,y,z+1)+","
                geo_faces += eg.edge(x,y,z+1,   x+1,y,z+1)+","
                geo_faces += eg.edge(x+1,y,z+1, x+1,y,z)+","
                geo_faces += eg.edge(x+1,y,z,   x,y,z)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1
            if not is_solid(x, y+1, z):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x,y+1,z,     x+1,y+1,z)+","
                geo_faces += eg.edge(x+1,y+1,z,   x+1,y+1,z+1)+","
                geo_faces += eg.edge(x+1,y+1,z+1, x,y+1,z+1)+","
                geo_faces += eg.edge(x,y+1,z+1,   x,y+1,z)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1
            if not is_solid(x, y, z-1):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x,y,z,     x+1,y,z)+","
                geo_faces += eg.edge(x+1,y,z,   x+1,y+1,z)+","
                geo_faces += eg.edge(x+1,y+1,z, x,y+1,z)+","
                geo_faces += eg.edge(x,y+1,z,   x,y,z)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1
            if not is_solid(x, y, z+1):
                geo_faces += "Line Loop("+str(loop_n)+")={"
                geo_faces += eg.edge(x,y,z+1,     x,y+1,z+1)+","
                geo_faces += eg.edge(x,y+1,z+1,   x+1,y+1,z+1)+","
                geo_faces += eg.edge(x+1,y+1,z+1, x+1,y,z+1)+","
                geo_faces += eg.edge(x+1,y,z+1,   x,y,z+1)
                geo_faces += "};\n"
                geo_surfaces += "Plane Surface("+str(loop_n)+")="
                geo_surfaces += "{"+str(loop_n)+"};\n"
                loop_n += 1

    for x, y, z in product(range(N+1),repeat=3):
        if get_point(x, y, z) in eg.edges:
            geo += "Point(" + str(get_point(x, y, z)) + ")="
            geo += "{" + coords(x,y,z) + "};\n"
    geo += eg.geo_edges + "\n"
    geo += geo_faces + "\n"
    geo += geo_surfaces + "\n"

    geo += "\nMesh.Algorithm = 6;"
    with open("/home/matt/python/bempp-test/mesh.geo","w") as f:
        f.write(geo)
    return __generate_grid_from_geo_string(geo)


# 2D fractals

def sierpinski_triangle(h=.1, level=2):
    """
    Returns a Sierpinski triangle.

    Parameters
    ----------
    h : float
        Size of the elements.
    level : int
        Refinement level.

    """
    from numpy import sqrt, array
    if level < 1:
        raise ValueError("level must be 1 or larger")
    h = min(.5**(level + 1), h)
    geo = "lc = " + str(h) + ";\n"
    points = [array((0, 0)), array((1, 0)), array((.5, sqrt(3) / 2))]
    tris = [[0, 1, 2]]
    for _ in range(level):
        new_tris = []
        for t in tris:
            new_p = []
            for i in range(3):
                new_p.append(len(points))
                points.append(
                    (points[
                        t[0]] +
                     points[
                         t[1]] +
                     points[
                         t[2]] -
                     points[
                         t[i]]) /
                    2.)
            new_tris.append([t[0], new_p[2], new_p[1]])
            new_tris.append([new_p[2], t[1], new_p[0]])
            new_tris.append([new_p[1], new_p[0], t[2]])
        tris = new_tris
    for i, p in enumerate(points):
        geo += "Point(" + str(1 + i) + \
            ") = {" + str(p[0]) + "," + str(p[1]) + ",0,lc};\n"
    for n, t in enumerate(tris):
        geo += "Line(" + str(1 + 3 * n) + \
            ") = {" + str(1 + t[0]) + "," + str(1 + t[1]) + "};\n"
        geo += "Line(" + str(1 + 3 * n + 1) + \
            ") = {" + str(1 + t[1]) + "," + str(1 + t[2]) + "};\n"
        geo += "Line(" + str(1 + 3 * n + 2) + \
            ") = {" + str(1 + t[2]) + "," + str(1 + t[0]) + "};\n"
        geo += "Line Loop(" + str(1 + n) + ") = {" + str(1 + 3 * n) + "," + str(
            1 + 3 * n + 1) + "," + str(1 + 3 * n + 2) + "};\n"
        geo += "Plane Surface(" + str(1 + n) + ") = { " + str(1 + n) + "};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)


def sierpinski_carpet(h=.1, level=2):
    """
    Returns a Sierpinski carpet.

    Parameters
    ----------
    h : float
        Size of the elements.
    level : int
        Refinement level.

    """
    from numpy import array
    from itertools import product
    geo = "lc = " + str(h) + ";\n"
    if level < 1:
        raise ValueError("level must be 1 or larger")
    points = [array((0, 0)), array((1, 0)), array((1, 1)), array((0, 1))]
    squs = [(0, 1, 2, 3)]
    for i in range(level):
        new_squs = []
        for c in squs:
            n = []
            for x, y in product((0, 1. / 3, 2. / 3, 1), repeat=2):
                if x not in [0, 1] or y not in [0, 1]:
                    n.append(len(points))
                    points.append(points[c[0]]
                                  + x * (points[c[1]] - points[c[0]])
                                  + y * (points[c[3]] - points[c[0]])
                                 )
            new_squs.append((n[1], n[4], n[5], c[3]))
            new_squs.append((n[0], n[3], n[4], n[1]))
            new_squs.append((c[0], n[2], n[3], n[0]))
            new_squs.append((n[2], n[6], n[7], n[3]))
            new_squs.append((n[6], c[1], n[10], n[7]))
            new_squs.append((n[7], n[10], n[11], n[8]))
            new_squs.append((n[8], n[11], c[2], n[9]))
            new_squs.append((n[4], n[8], n[9], n[5]))
        squs = new_squs
    for i, p in enumerate(points):
        geo += "Point(" + str(1 + i) + \
            ") = {" + str(p[0]) + "," + str(p[1]) + ",0,lc};\n"
    line_n = 1
    all_edges = {}
    for n, t in enumerate(squs):
        edges = []
        for a in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            v1 = 1 + t[a[0]]
            v2 = 1 + t[a[1]]
            if v1 in all_edges and v2 in all_edges[v1]:
                edges.append(all_edges[v1][v2])
            elif v2 in all_edges and v1 in all_edges[v2]:
                edges.append(-all_edges[v2][v1])
            else:
                geo += "Line(" + str(line_n) + \
                    ") = {" + str(v1) + "," + str(v2) + "};\n"
                if v1 not in all_edges:
                    all_edges[v1] = {}
                all_edges[v1][v2] = line_n
                edges.append(line_n)
                line_n += 1
        geo += "Line Loop(" + str(n) + ") = { " + str(edges[0]) + ", " + str(
            edges[1]) + ", " + str(edges[2]) + ", " + str(edges[3]) + "};\n"
        geo += "Plane Surface(" + str(n) + ") = {" + str(n) + "};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)


def koch_snowflake(h=.1, level=2):
    """
    Returns a Koch snowflake.

    Parameters
    ----------
    h : float
        Size of the elements.
    level : int
        Refinement level.

    """
    from numpy import sqrt, array
    if level < 1:
        raise ValueError("level must be 1 or larger")
    h = min(.5**(level + 1), h)
    geo = "lc = " + str(h) + ";\n"
    points = [array((0, 0)), array((1, 0)), array((.5, sqrt(3) / 2))]
    lines = [(0, 1), (1, 2), (2, 0)]

    for l in range(level):
        new_lines = []
        for line in lines:
            new_lines.append((line[0], len(points)))
            new_lines.append((len(points), len(points) + 1))
            new_lines.append((len(points) + 1, len(points) + 2))
            new_lines.append((len(points) + 2, line[1]))

            edge = points[line[1]] - points[line[0]]
            normal = sqrt(3) / 6 * array((-edge[1], edge[0]))
            points.append(points[line[0]] + edge / 3.)
            points.append(points[line[0]] + (edge) / 2. - normal)
            points.append(points[line[0]] + 2 * edge / 3.)

        lines = new_lines

    for i, p in enumerate(points):
        geo += "Point(" + str(1 + i) + \
            ") = {" + str(p[0]) + "," + str(p[1]) + ",0,lc};\n"
    for i, l in enumerate(lines):
        geo += "Line(" + str(1 + i) + \
            ") = {" + str(l[0] + 1) + "," + str(l[1] + 1) + "};\n"
    geo += "Line Loop(1) = { " + ",".join([str(i + 1)
                                           for i in range(len(lines))]) + "};\n"
    geo += "Plane Surface(1) = {1};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)
