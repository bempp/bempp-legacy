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



def _b(x=0,y=0):
    """Return gmsh code for a letter B with bottom left corner at (x,y)."""
    global _l_s
    out = "ls = "+str(_l_s)+""";
           Point(ls+1) = {"""+str(x)+","+str(y)+""",1,lc};
           Point(ls+2) = {"""+str(x+1)+","+str(y)+""",1,lc};
           Point(ls+3) = {"""+str(x+1)+","+str(y+1)+""",1,lc};
           Point(ls+4) = {"""+str(x+1.28)+","+str(y+1.96)+""",1,lc};
           Point(ls+5) = {"""+str(x+1)+","+str(y+2.92)+""",1,lc};
           Point(ls+6) = {"""+str(x+1)+","+str(y+3.92)+""",1,lc};
           Point(ls+7) = {"""+str(x)+","+str(y+3.92)+""",1,lc};
           Point(ls+8) = {"""+str(x+.5)+","+str(y+.5)+""",1,lc};
           Point(ls+9) = {"""+str(x+.9)+","+str(y+.5)+""",1,lc};
           Point(ls+10) = {"""+str(x+.9)+","+str(y+1)+""",1,lc};
           Point(ls+11) = {"""+str(x+.9)+","+str(y+1.5)+""",1,lc};
           Point(ls+12) = {"""+str(x+.5)+","+str(y+1.5)+""",1,lc};
           Point(ls+13) = {"""+str(x+.5)+","+str(y+2.42)+""",1,lc};
           Point(ls+14) = {"""+str(x+.9)+","+str(y+2.42)+""",1,lc};
           Point(ls+15) = {"""+str(x+.9)+","+str(y+2.92)+""",1,lc};
           Point(ls+16) = {"""+str(x+.9)+","+str(y+3.42)+""",1,lc};
           Point(ls+17) = {"""+str(x+.5)+","+str(y+3.42)+""",1,lc};
           Point(ls+18) = {"""+str(x)+","+str(y)+""",0,lc};
           Point(ls+19) = {"""+str(x+1)+","+str(y)+""",0,lc};
           Point(ls+20) = {"""+str(x+1)+","+str(y+1)+""",0,lc};
           Point(ls+21) = {"""+str(x+1.28)+","+str(y+1.96)+""",0,lc};
           Point(ls+22) = {"""+str(x+1)+","+str(y+2.92)+""",0,lc};
           Point(ls+23) = {"""+str(x+1)+","+str(y+3.92)+""",0,lc};
           Point(ls+24) = {"""+str(x)+","+str(y+3.92)+""",0,lc};
           Point(ls+25) = {"""+str(x+.5)+","+str(y+.5)+""",0,lc};
           Point(ls+26) = {"""+str(x+.9)+","+str(y+.5)+""",0,lc};
           Point(ls+27) = {"""+str(x+.9)+","+str(y+1)+""",0,lc};
           Point(ls+28) = {"""+str(x+.9)+","+str(y+1.5)+""",0,lc};
           Point(ls+29) = {"""+str(x+.5)+","+str(y+1.5)+""",0,lc};
           Point(ls+30) = {"""+str(x+.5)+","+str(y+2.42)+""",0,lc};
           Point(ls+31) = {"""+str(x+.9)+","+str(y+2.42)+""",0,lc};
           Point(ls+32) = {"""+str(x+.9)+","+str(y+2.92)+""",0,lc};
           Point(ls+33) = {"""+str(x+.9)+","+str(y+3.42)+""",0,lc};
           Point(ls+34) = {"""+str(x+.5)+","+str(y+3.42)+""",0,lc};

           Line(ls+1) = {ls+1,ls+2};
           Circle(ls+2) = {ls+2,ls+3,ls+4};
           Circle(ls+3) = {ls+4,ls+5,ls+6};
           Line(ls+4) = {ls+6,ls+7};
           Line(ls+5) = {ls+7,ls+1};
           Line(ls+6) = {ls+8,ls+9};
           Circle(ls+7) = {ls+9,ls+10,ls+11};
           Line(ls+8) = {ls+11,ls+12};
           Line(ls+9) = {ls+12,ls+8};
           Line(ls+10) = {ls+13,ls+14};
           Circle(ls+11) = {ls+14,ls+15,ls+16};
           Line(ls+12) = {ls+16,ls+17};
           Line(ls+13) = {ls+17,ls+13};
           Line(ls+14) = {ls+18,ls+19};
           Circle(ls+15) = {ls+19,ls+20,ls+21};
           Circle(ls+16) = {ls+21,ls+22,ls+23};
           Line(ls+17) = {ls+23,ls+24};
           Line(ls+18) = {ls+24,ls+18};
           Line(ls+19) = {ls+25,ls+26};
           Circle(ls+20) = {ls+26,ls+27,ls+28};
           Line(ls+21) = {ls+28,ls+29};
           Line(ls+22) = {ls+29,ls+25};
           Line(ls+23) = {ls+30,ls+31};
           Circle(ls+24) = {ls+31,ls+32,ls+33};
           Line(ls+25) = {ls+33,ls+34};
           Line(ls+26) = {ls+34,ls+30};
           Line(ls+27) = {ls+1,ls+18};
           Line(ls+28) = {ls+7,ls+24};
           Line(ls+29) = {ls+8,ls+25};
           Line(ls+30) = {ls+12,ls+29};
           Line(ls+31) = {ls+13,ls+30};
           Line(ls+32) = {ls+17,ls+34};
           Line(ls+33) = {ls+4,ls+21};
           Line(ls+34) = {ls+2,ls+19};
           Line(ls+35) = {ls+6,ls+23};
           Line(ls+36) = {ls+16,ls+33};
           Line(ls+37) = {ls+14,ls+31};
           Line(ls+38) = {ls+11,ls+28};
           Line(ls+39) = {ls+9,ls+26};

           Line Loop(ls+1) = {ls+1,ls+2,ls+3,ls+4,ls+5};
           Line Loop(ls+2) = {ls+6,ls+7,ls+8,ls+9};
           Line Loop(ls+3) = {ls+10,ls+11,ls+12,ls+13};
           Line Loop(ls+4) = {-ls-29,ls+6,ls+39,-ls-19};
           Line Loop(ls+5) = {ls+30,-ls-21,-ls-38,ls+8};
           Line Loop(ls+6) = {-ls-31,ls+10,ls+37,-ls-23};
           Line Loop(ls+7) = {ls+32,-ls-25,-ls-36,ls+12};
           Line Loop(ls+8) = {-ls-37,ls+11,ls+36,-ls-24};
           Line Loop(ls+9) = {-ls-39,ls+7,ls+38,-ls-20};
           Line Loop(ls+10) = {ls+27,ls+14,-ls-34,-ls-1};
           Line Loop(ls+11) = {-ls-28,-ls-4,ls+35,ls+17};
           Line Loop(ls+12) = {ls+34,ls+15,-ls-33,-ls-2};
           Line Loop(ls+13) = {-ls-35,-ls-3,ls+33,ls+16};
           Line Loop(ls+14) = {-ls-18,-ls-17,-ls-16,-ls-15,-ls-14};
           Line Loop(ls+15) = {-ls-22,-ls-21,-ls-20,-ls-19};
           Line Loop(ls+16) = {-ls-26,-ls-25,-ls-24,-ls-23};
           Line Loop(ls+17) = {ls+18,-ls-27,-ls-5,ls+28};
           Line Loop(ls+18) = {ls+9,-ls-30,-ls-22,ls+29};
           Line Loop(ls+19) = {ls+13,ls+31,-ls-26,-ls-32};

           Plane Surface(ls+1) = {ls+1,ls+2,ls+3};
           Plane Surface(ls+2) = {ls+14,ls+15,ls+16};
           Plane Surface(ls+3) = {ls+4};
           Plane Surface(ls+4) = {ls+7};
           Plane Surface(ls+5) = {ls+5};
           Plane Surface(ls+6) = {ls+6};
           Plane Surface(ls+7) = {ls+10};
           Plane Surface(ls+8) = {ls+11};
           Plane Surface(ls+9) = {ls+17};
           Plane Surface(ls+10) = {ls+18};
           Plane Surface(ls+11) = {ls+19};
           Ruled Surface(ls+12) = {ls+8};
           Ruled Surface(ls+13) = {ls+9};
           Ruled Surface(ls+14) = {ls+12};
           Ruled Surface(ls+15) = {ls+13};"""
    _l_s += 39
    return out

_twelve = """Line(ls+1) = {ls+1,ls+2};
           Line(ls+2) = {ls+2,ls+3};
           Line(ls+3) = {ls+3,ls+4};
           Line(ls+4) = {ls+4,ls+5};
           Line(ls+5) = {ls+5,ls+6};
           Line(ls+6) = {ls+6,ls+7};
           Line(ls+7) = {ls+7,ls+8};
           Line(ls+8) = {ls+8,ls+9};
           Line(ls+9) = {ls+9,ls+10};
           Line(ls+10) = {ls+10,ls+11};
           Line(ls+11) = {ls+11,ls+12};
           Line(ls+12) = {ls+12,ls+1};

           Line(ls+13) = {ls+13,ls+14};
           Line(ls+14) = {ls+14,ls+15};
           Line(ls+15) = {ls+15,ls+16};
           Line(ls+16) = {ls+16,ls+17};
           Line(ls+17) = {ls+17,ls+18};
           Line(ls+18) = {ls+18,ls+19};
           Line(ls+19) = {ls+19,ls+20};
           Line(ls+20) = {ls+20,ls+21};
           Line(ls+21) = {ls+21,ls+22};
           Line(ls+22) = {ls+22,ls+23};
           Line(ls+23) = {ls+23,ls+24};
           Line(ls+24) = {ls+24,ls+13};

           Line(ls+25) = {ls+1,ls+13};
           Line(ls+26) = {ls+2,ls+14};
           Line(ls+27) = {ls+3,ls+15};
           Line(ls+28) = {ls+4,ls+16};
           Line(ls+29) = {ls+5,ls+17};
           Line(ls+30) = {ls+6,ls+18};
           Line(ls+31) = {ls+7,ls+19};
           Line(ls+32) = {ls+8,ls+20};
           Line(ls+33) = {ls+9,ls+21};
           Line(ls+34) = {ls+10,ls+22};
           Line(ls+35) = {ls+11,ls+23};
           Line(ls+36) = {ls+12,ls+24};


           Line Loop(ls+1) = {ls+1,ls+2,ls+3,ls+4,ls+5,ls+6,ls+7,ls+8,ls+9,ls+10,ls+11,ls+12};
           Line Loop(ls+2) = {-ls-24,-ls-23,-ls-22,-ls-21,-ls-20,-ls-19,-ls-18,-ls-17,-ls-16,-ls-15,-ls-14,-ls-13};

           Line Loop(ls+3) = {-ls-1,ls+25,ls+13,-ls-26};
           Line Loop(ls+4) = {-ls-2,ls+26,ls+14,-ls-27};
           Line Loop(ls+5) = {-ls-3,ls+27,ls+15,-ls-28};
           Line Loop(ls+6) = {-ls-4,ls+28,ls+16,-ls-29};
           Line Loop(ls+7) = {-ls-5,ls+29,ls+17,-ls-30};
           Line Loop(ls+8) = {-ls-6,ls+30,ls+18,-ls-31};
           Line Loop(ls+9) = {-ls-7,ls+31,ls+19,-ls-32};
           Line Loop(ls+10) = {-ls-8,ls+32,ls+20,-ls-33};
           Line Loop(ls+11) = {-ls-9,ls+33,ls+21,-ls-34};
           Line Loop(ls+12) = {-ls-10,ls+34,ls+22,-ls-35};
           Line Loop(ls+13) = {-ls-11,ls+35,ls+23,-ls-36};
           Line Loop(ls+14) = {-ls-12,ls+36,ls+24,-ls-25};

           Plane Surface(ls+1) = {ls+1};
           Plane Surface(ls+2) = {ls+2};
           Plane Surface(ls+3) = {ls+3};
           Plane Surface(ls+4) = {ls+4};
           Plane Surface(ls+5) = {ls+5};
           Plane Surface(ls+6) = {ls+6};
           Plane Surface(ls+7) = {ls+7};
           Plane Surface(ls+8) = {ls+8};
           Plane Surface(ls+9) = {ls+9};
           Plane Surface(ls+10) = {ls+10};
           Plane Surface(ls+11) = {ls+11};
           Plane Surface(ls+12) = {ls+12};
           Plane Surface(ls+14) = {ls+14};
           Plane Surface(ls+13) = {ls+13};"""


def _e(x=0,y=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    global _l_s
    out = "ls = "+str(_l_s)+""";
           Point(ls+1) = {"""+str(x)+","+str(y)+""",1,lc};
           Point(ls+2) = {"""+str(x+2)+","+str(y)+""",1,lc};
           Point(ls+3) = {"""+str(x+2)+","+str(y+.784)+""",1,lc};
           Point(ls+4) = {"""+str(x+.5)+","+str(y+.784)+""",1,lc};
           Point(ls+5) = {"""+str(x+.5)+","+str(y+1.568)+""",1,lc};
           Point(ls+6) = {"""+str(x+1.5)+","+str(y+1.568)+""",1,lc};
           Point(ls+7) = {"""+str(x+1.5)+","+str(y+2.352)+""",1,lc};
           Point(ls+8) = {"""+str(x+.5)+","+str(y+2.352)+""",1,lc};
           Point(ls+9) = {"""+str(x+.5)+","+str(y+3.136)+""",1,lc};
           Point(ls+10) = {"""+str(x+2)+","+str(y+3.136)+""",1,lc};
           Point(ls+11) = {"""+str(x+2)+","+str(y+3.92)+""",1,lc};
           Point(ls+12) = {"""+str(x)+","+str(y+3.92)+""",1,lc};

           Point(ls+13) = {"""+str(x)+","+str(y)+""",0,lc};
           Point(ls+14) = {"""+str(x+2)+","+str(y)+""",0,lc};
           Point(ls+15) = {"""+str(x+2)+","+str(y+.784)+""",0,lc};
           Point(ls+16) = {"""+str(x+.5)+","+str(y+.784)+""",0,lc};
           Point(ls+17) = {"""+str(x+.5)+","+str(y+1.568)+""",0,lc};
           Point(ls+18) = {"""+str(x+1.5)+","+str(y+1.568)+""",0,lc};
           Point(ls+19) = {"""+str(x+1.5)+","+str(y+2.352)+""",0,lc};
           Point(ls+20) = {"""+str(x+.5)+","+str(y+2.352)+""",0,lc};
           Point(ls+21) = {"""+str(x+.5)+","+str(y+3.136)+""",0,lc};
           Point(ls+22) = {"""+str(x+2)+","+str(y+3.136)+""",0,lc};
           Point(ls+23) = {"""+str(x+2)+","+str(y+3.92)+""",0,lc};
           Point(ls+24) = {"""+str(x)+","+str(y+3.92)+""",0,lc};
           """ + _twelve

    _l_s += 36
    return out

def _m(x=0,y=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    global _l_s
    out = "ls = "+str(_l_s)+""";
           Point(ls+1) = {"""+str(x)+","+str(y)+""",1,lc};
           Point(ls+2) = {"""+str(x+.5)+","+str(y)+""",1,lc};
           Point(ls+3) = {"""+str(x+.5)+","+str(y+3)+""",1,lc};
           Point(ls+4) = {"""+str(x+1.5)+","+str(y+1)+""",1,lc};
           Point(ls+5) = {"""+str(x+2.5)+","+str(y+3)+""",1,lc};
           Point(ls+6) = {"""+str(x+2.5)+","+str(y)+""",1,lc};
           Point(ls+7) = {"""+str(x+3)+","+str(y)+""",1,lc};
           Point(ls+8) = {"""+str(x+3)+","+str(y+3.92)+""",1,lc};
           Point(ls+9) = {"""+str(x+2.3)+","+str(y+3.92)+""",1,lc};
           Point(ls+10) = {"""+str(x+1.5)+","+str(y+2.32)+""",1,lc};
           Point(ls+11) = {"""+str(x+0.9)+","+str(y+3.92)+""",1,lc};
           Point(ls+12) = {"""+str(x)+","+str(y+3.92)+""",1,lc};

           Point(ls+13) = {"""+str(x)+","+str(y)+""",0,lc};
           Point(ls+14) = {"""+str(x+.5)+","+str(y)+""",0,lc};
           Point(ls+15) = {"""+str(x+.5)+","+str(y+3)+""",0,lc};
           Point(ls+16) = {"""+str(x+1.5)+","+str(y+1)+""",0,lc};
           Point(ls+17) = {"""+str(x+2.5)+","+str(y+3)+""",0,lc};
           Point(ls+18) = {"""+str(x+2.5)+","+str(y)+""",0,lc};
           Point(ls+19) = {"""+str(x+3)+","+str(y)+""",0,lc};
           Point(ls+20) = {"""+str(x+3)+","+str(y+3.92)+""",0,lc};
           Point(ls+21) = {"""+str(x+2.3)+","+str(y+3.92)+""",0,lc};
           Point(ls+22) = {"""+str(x+1.5)+","+str(y+2.32)+""",0,lc};
           Point(ls+23) = {"""+str(x+0.9)+","+str(y+3.92)+""",0,lc};
           Point(ls+24) = {"""+str(x)+","+str(y+3.92)+""",0,lc};
           """ + _twelve

    _l_s += 36
    return out

def _plus(x=0,y=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    global _l_s
    out = "ls = "+str(_l_s)+""";
           Point(ls+1) = {"""+str(x)+","+str(y+1.71)+""",1,lc};
           Point(ls+2) = {"""+str(x+.5)+","+str(y+1.71)+""",1,lc};
           Point(ls+3) = {"""+str(x+.5)+","+str(y+1.21)+""",1,lc};
           Point(ls+4) = {"""+str(x+1)+","+str(y+1.21)+""",1,lc};
           Point(ls+5) = {"""+str(x+1)+","+str(y+1.71)+""",1,lc};
           Point(ls+6) = {"""+str(x+1.5)+","+str(y+1.71)+""",1,lc};
           Point(ls+7) = {"""+str(x+1.5)+","+str(2.21)+""",1,lc};
           Point(ls+8) = {"""+str(x+1)+","+str(y+2.21)+""",1,lc};
           Point(ls+9) = {"""+str(x+1)+","+str(y+2.71)+""",1,lc};
           Point(ls+10) = {"""+str(x+.5)+","+str(y+2.71)+""",1,lc};
           Point(ls+11) = {"""+str(x+.5)+","+str(y+2.21)+""",1,lc};
           Point(ls+12) = {"""+str(x)+","+str(y+2.21)+""",1,lc};

           Point(ls+13) = {"""+str(x)+","+str(y+1.71)+""",0,lc};
           Point(ls+14) = {"""+str(x+.5)+","+str(y+1.71)+""",0,lc};
           Point(ls+15) = {"""+str(x+.5)+","+str(y+1.21)+""",0,lc};
           Point(ls+16) = {"""+str(x+1)+","+str(y+1.21)+""",0,lc};
           Point(ls+17) = {"""+str(x+1)+","+str(y+1.71)+""",0,lc};
           Point(ls+18) = {"""+str(x+1.5)+","+str(y+1.71)+""",0,lc};
           Point(ls+19) = {"""+str(x+1.5)+","+str(2.21)+""",0,lc};
           Point(ls+20) = {"""+str(x+1)+","+str(y+2.21)+""",0,lc};
           Point(ls+21) = {"""+str(x+1)+","+str(y+2.71)+""",0,lc};
           Point(ls+22) = {"""+str(x+.5)+","+str(y+2.71)+""",0,lc};
           Point(ls+23) = {"""+str(x+.5)+","+str(y+2.21)+""",0,lc};
           Point(ls+24) = {"""+str(x)+","+str(y+2.21)+""",0,lc};
           """ + _twelve

    _l_s += 36
    return out

