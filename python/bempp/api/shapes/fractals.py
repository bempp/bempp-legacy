from .shapes import __generate_grid_from_geo_string

def _make_word(geo,h):
    return __generate_grid_from_geo_string(intro+geo+outro)

# 3D fractals

def sierpinski_pyramid(h=0.1, level=2):
    from numpy import sqrt, array
    if level < 1:
        raise ValueError("level must be 1 or larger")
    geo = "lc = "+str(h)+";\n"
    points = [array((.5,1/(2*sqrt(3)),sqrt(2./3.))),array((0,0,0)),array((1,0,0)),array((.5,sqrt(3)/2,0))]
    tetras = [[0,1,2,3]]
    for l in range(level):
        new_tetras = []
        vertex_map = {}
        for t in tetras:
            for i in t:
                new_t = []
                for j in t:
                    if i==j:
                        new_t.append(i)
                    else:
                        a = min(i,j)
                        b = max(i,j)
                        if b not in vertex_map or a not in vertex_map[b]:
                            if b not in vertex_map:
                                vertex_map[b] = {}
                            vertex_map[b][a] = len(points)
                            points.append((points[a]+points[b])/2.)
                        new_t.append(vertex_map[b][a])
                new_tetras.append(new_t)
        tetras = new_tetras
        
    for i,p in enumerate(points):
        geo += "Point("+str(1+i)+") = {"+str(p[0])+","+str(p[1])+","+str(p[2])+",lc};\n"
    for n,t in enumerate(tetras):
        geo += "Line("+str(1+6*n  )+") = {"+str(1+t[0])+","+str(1+t[1])+"};\n"
        geo += "Line("+str(1+6*n+1)+") = {"+str(1+t[0])+","+str(1+t[2])+"};\n"
        geo += "Line("+str(1+6*n+2)+") = {"+str(1+t[0])+","+str(1+t[3])+"};\n"
        geo += "Line("+str(1+6*n+3)+") = {"+str(1+t[1])+","+str(1+t[2])+"};\n"
        geo += "Line("+str(1+6*n+4)+") = {"+str(1+t[2])+","+str(1+t[3])+"};\n"
        geo += "Line("+str(1+6*n+5)+") = {"+str(1+t[3])+","+str(1+t[1])+"};\n"
        geo += "Line Loop("+str(1+4*n  )+") = { "+str(1+6*n  )+", "+str(1+6*n+3)+",-"+str(1+6*n+1)+"};\n"
        geo += "Line Loop("+str(1+4*n+1)+") = { "+str(1+6*n+1)+", "+str(1+6*n+4)+",-"+str(1+6*n+2)+"};\n"
        geo += "Line Loop("+str(1+4*n+2)+") = { "+str(1+6*n+2)+", "+str(1+6*n+5)+",-"+str(1+6*n  )+"};\n"
        geo += "Line Loop("+str(1+4*n+3)+") = {-"+str(1+6*n+4)+",-"+str(1+6*n+3)+",-"+str(1+6*n+5)+"};\n"
        for i in range(4):
            geo += "Plane Surface("+str(1+4*n+i)+") = {"+str(1+4*n+i)+"};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)

def menger_sponge(h=0.1, level=2):
    from numpy import array
    from itertools import product
    geo = "lc = "+str(h)+";\n"
    if level < 1:
        raise ValueError("level must be 1 or larger")
    points = [array((0,0,1)),array((1,0,1)),array((1,1,1)),array((0,1,1)),array((0,0,0)),array((1,0,0)),array((1,1,0)),array((0,1,0))]
    cubes = [(0,1,2,3,4,5,6,7,(True,True,True,True,True,True))]
    for i in range(level):
        new_cubes = []
        for c in cubes:
            n = []
            for x,y,z in product((0,1./3,2./3,1),repeat=3):
                if x not in [0,1] or y not in [0,1] or z not in [0,1]:
                    n.append(len(points))
                    points.append(points[c[4]]
                                  + x*(points[c[5]]-points[c[4]])
                                  + y*(points[c[7]]-points[c[4]])
                                  + z*(points[c[0]]-points[c[4]])
                                 )
            new_cubes.append((c[0],n[15],n[19],n[5],n[1],n[14],n[18],n[4],     (c[8][0],c[8][1],False,False,c[8][4],False)))
            new_cubes.append((n[31],c[1],n[49],n[35],n[30],n[45],n[48],n[34],  (c[8][0],c[8][1],c[8][2],False,False,False)))
            new_cubes.append((n[39],n[53],c[2],n[43],n[38],n[52],n[55],n[42],  (c[8][0],False,c[8][2],c[8][3],False,False)))
            new_cubes.append((n[9],n[23],n[27],c[3],n[8],n[22],n[26],n[11],    (c[8][0],False,False,c[8][3],c[8][4],False)))

            new_cubes.append((n[5],n[19],n[23],n[9],n[4],n[18],n[22],n[8],     (c[8][0],False,True,False,c[8][4],True)))
            new_cubes.append((n[15],n[31],n[35],n[19],n[14],n[30],n[34],n[18], (c[8][0],c[8][1],False,True,False,True)))
            new_cubes.append((n[35],n[49],n[53],n[39],n[34],n[48],n[52],n[38], (c[8][0],False,c[8][2],False,True,True)))
            new_cubes.append((n[23],n[39],n[43],n[27],n[22],n[38],n[42],n[26], (c[8][0],True,False,c[8][3],False,True)))

            new_cubes.append((n[1],n[14],n[18],n[4],n[0],n[13],n[17],n[3],     (False,c[8][1],True,True,c[8][4],False)))
            new_cubes.append((n[30],n[45],n[48],n[34],n[29],n[44],n[47],n[33], (False,c[8][1],c[8][2],True,True,False)))
            new_cubes.append((n[38],n[52],n[55],n[42],n[37],n[51],n[54],n[41], (False,True,c[8][2],c[8][3],True,False)))
            new_cubes.append((n[8],n[22],n[26],n[11],n[7],n[21],n[25],n[10],   (False,True,True,c[8][3],c[8][4],False)))

            new_cubes.append((n[0],n[13],n[17],n[3],c[4],n[12],n[16],n[2],     (False,c[8][1],False,False,c[8][4],c[8][5])))
            new_cubes.append((n[29],n[44],n[47],n[33],n[28],c[5],n[46],n[32],  (False,c[8][1],c[8][2],False,False,c[8][5])))
            new_cubes.append((n[37],n[51],n[54],n[41],n[36],n[50],c[6],n[40],  (False,False,c[8][2],c[8][3],False,c[8][5])))
            new_cubes.append((n[7],n[21],n[25],n[10],n[6],n[20],n[24],c[7],    (False,False,False,c[8][3],c[8][4],c[8][5])))

            new_cubes.append((n[3],n[17],n[21],n[7],n[2],n[16],n[20],n[6],     (True,False,True,False,c[8][4],c[8][5])))
            new_cubes.append((n[13],n[29],n[33],n[17],n[12],n[28],n[32],n[16], (True,c[8][1],False,True,False,c[8][5])))
            new_cubes.append((n[33],n[47],n[51],n[37],n[32],n[46],n[50],n[36], (True,False,c[8][2],False,True,c[8][5])))
            new_cubes.append((n[21],n[37],n[41],n[25],n[20],n[36],n[40],n[24], (True,True,False,c[8][3],False,c[8][5])))
        cubes = new_cubes
    for i,p in enumerate(points):
        geo += "Point("+str(1+i)+") = {"+str(p[0])+","+str(p[1])+","+str(p[2])+",lc};\n"
    loop_n = 1
    line_n = 1
    all_edges = {}
    for n,t in enumerate(cubes):
        edges = []
        for a in [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]:
            v1 = 1+t[a[0]]
            v2 = 1+t[a[1]]
            if v1 in all_edges and v2 in all_edges[v1]:
                edges.append(all_edges[v1][v2])
            elif v2 in all_edges and v1 in all_edges[v2]:
                edges.append(-all_edges[v2][v1])
            else:
                geo += "Line("+str(line_n)+") = {"+str(v1)+","+str(v2)+"};\n"
                if v1 not in all_edges:
                    all_edges[v1] = {}
                all_edges[v1][v2] = line_n
                edges.append(line_n)
                line_n += 1

        if t[8][0]:
            geo += "Line Loop("+str(loop_n)+") = { "+str(edges[0])+", "+str(edges[1])+", "+str(edges[2])+", "+str(edges[3])+"};\n"
            loop_n += 1
        if t[8][1]:
            geo += "Line Loop("+str(loop_n)+") = {"+str(-edges[0])+", "+str(edges[4])+", "+str(edges[8])+","+str(-edges[5])+"};\n"
            loop_n += 1
        if t[8][2]:
            geo += "Line Loop("+str(loop_n)+") = {"+str(-edges[1])+", "+str(edges[5])+", "+str(edges[9])+","+str(-edges[6])+"};\n"
            loop_n += 1
        if t[8][3]:
            geo += "Line Loop("+str(loop_n)+") = {"+str(-edges[2])+", "+str(edges[6])+", "+str(edges[10])+","+str(-edges[7])+"};\n"
            loop_n += 1
        if t[8][4]:
            geo += "Line Loop("+str(loop_n)+") = {"+str(-edges[3])+", "+str(edges[7])+", "+str(edges[11])+","+str(-edges[4])+"};\n"
            loop_n += 1
        if t[8][5]:
            geo += "Line Loop("+str(loop_n)+") = {"+str(-edges[11])+","+str(-edges[10])+","+str(-edges[9])+","+str(-edges[8])+"};\n"
            loop_n += 1
    for i in range(1,loop_n):
        geo += "Plane Surface("+str(i)+") = {"+str(i)+"};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)


# 2D fractals

def sierpinski_triangle(h=.1, level=2):
    from numpy import sqrt, array
    if level < 1:
        raise ValueError("level must be 1 or larger")
    h = min(.5**(level+1),h)
    geo = "lc = "+str(h)+";\n"
    points = [array((0,0)),array((1,0)),array((.5,sqrt(3)/2))]
    tris = [[0,1,2]]
    for l in range(level):
        new_tris = []
        vertex_map = {}
        for t in tris:
            new_p = []
            for i in range(3):
                new_p.append(len(points))
                points.append((points[t[0]]+points[t[1]]+points[t[2]]-points[t[i]])/2.)
            new_tris.append([t[0],new_p[2],new_p[1]])
            new_tris.append([new_p[2],t[1],new_p[0]])
            new_tris.append([new_p[1],new_p[0],t[2]])
        tris = new_tris
    for i,p in enumerate(points):
        geo += "Point("+str(1+i)+") = {"+str(p[0])+","+str(p[1])+",0,lc};\n"
    for n,t in enumerate(tris):
        geo += "Line("+str(1+3*n  )+") = {"+str(1+t[0])+","+str(1+t[1])+"};\n"
        geo += "Line("+str(1+3*n+1)+") = {"+str(1+t[1])+","+str(1+t[2])+"};\n"
        geo += "Line("+str(1+3*n+2)+") = {"+str(1+t[2])+","+str(1+t[0])+"};\n"
        geo += "Line Loop("+str(1+n)+") = {"+str(1+3*n)+","+str(1+3*n+1)+","+str(1+3*n+2)+"};\n"
        geo += "Plane Surface("+str(1+n)+") = { "+str(1+n)+"};\n"
    geo += "\nMesh.Algorithm = 6;"
    with open("/home/matt/python/bempp/fractals/tri.geo","w") as f:
        f.write(geo)
    return __generate_grid_from_geo_string(geo)

def sierpinski_carpet(h=.1, level=2):
    from numpy import array
    from itertools import product
    geo = "lc = "+str(h)+";\n"
    if level < 1:
        raise ValueError("level must be 1 or larger")
    points = [array((0,0)),array((1,0)),array((1,1)),array((0,1))]
    squs = [(0,1,2,3)]
    for i in range(level):
        new_squs = []
        for c in squs:
            n = []
            for x,y in product((0,1./3,2./3,1),repeat=2):
                if x not in [0,1] or y not in [0,1]:
                    n.append(len(points))
                    points.append(points[c[0]]
                                  + x*(points[c[1]]-points[c[0]])
                                  + y*(points[c[3]]-points[c[0]])
                                 )
            new_squs.append((n[1],n[4],n[5],c[3]))
            new_squs.append((n[0],n[3],n[4],n[1]))
            new_squs.append((c[0],n[2],n[3],n[0]))
            new_squs.append((n[2],n[6],n[7],n[3]))
            new_squs.append((n[6],c[1],n[10],n[7]))
            new_squs.append((n[7],n[10],n[11],n[8]))
            new_squs.append((n[8],n[11],c[2],n[9]))
            new_squs.append((n[4],n[8],n[9],n[5]))
        squs = new_squs
    for i,p in enumerate(points):
        geo += "Point("+str(1+i)+") = {"+str(p[0])+","+str(p[1])+",0,lc};\n"
    line_n = 1
    all_edges = {}
    for n,t in enumerate(squs):
        edges = []
        for a in [(0,1),(1,2),(2,3),(3,0)]:
            v1 = 1+t[a[0]]
            v2 = 1+t[a[1]]
            if v1 in all_edges and v2 in all_edges[v1]:
                edges.append(all_edges[v1][v2])
            elif v2 in all_edges and v1 in all_edges[v2]:
                edges.append(-all_edges[v2][v1])
            else:
                geo += "Line("+str(line_n)+") = {"+str(v1)+","+str(v2)+"};\n"
                if v1 not in all_edges:
                    all_edges[v1] = {}
                all_edges[v1][v2] = line_n
                edges.append(line_n)
                line_n += 1
        geo += "Line Loop("+str(n)+") = { "+str(edges[0])+", "+str(edges[1])+", "+str(edges[2])+", "+str(edges[3])+"};\n"
        geo += "Plane Surface("+str(n)+") = {"+str(n)+"};\n"
    geo += "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(geo)
