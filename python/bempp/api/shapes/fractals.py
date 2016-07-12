from .shapes import __generate_grid_from_geo_string

def _make_word(geo,h):
    return __generate_grid_from_geo_string(intro+geo+outro)

def sierpinski_pyramid(h=0.1, level=2):
    if level < 1:
        raise ValueError("level must be 1 or larger")
    from numpy import sqrt
    geo = "lc = "+str(h)+";\n"
    points = [(.5,1/(2*sqrt(3)),sqrt(2./3.)),(0,0,0),(1,0,0),(.5,sqrt(3)/2,0)]
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
                            points.append(tuple((p+q)/2. for p,q in zip(points[a],points[b])))
                        new_t.append(vertex_map[b][a])
                new_tetras.append(new_t)
        tetras = new_tetras
        
    for i,p in enumerate(points):
        print p
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
    with open("/home/matt/python/bempp/fractals/test.geo","w") as f:
        f.write(geo)
    return __generate_grid_from_geo_string(geo)
