"""Generate the BEM++ Letters."""

#pylint: disable=invalid-name
from bempp.api.shapes.shapes import __generate_grid_from_geo_string


def letters_bempp(h=0.1):
    """Create BEM++ as grid."""
    return _make_word(_b(0, 0, 0) + _e(2.5, 0, 39) + _m(5, 0, 75) +
                      _plus(8.5, 0, 111) + _plus(10.5, 0, 147), h)


def _make_word(geo, h):
    """Create the string for a mesh from components."""
    intro = "lc = " + str(h) + ";\n"
    outro = "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(intro + geo + outro)


def _b(x=0, y=0, l_s=0):
    """Return gmsh code for a letter B with bottom left corner at (x,y)."""
    out = "ls = " + str(l_s) + """;
       Point(ls+1) = {""" + str(x) + "," + str(y) + """,1,lc};
       Point(ls+2) = {""" + str(x + 1) + "," + str(y) + """,1,lc};
       Point(ls+3) = {""" + str(x + 1) + "," + str(y + 1) + """,1,lc};
       Point(ls+4) = {""" + str(x + 1.28) + "," + str(y + 1.96) + """,1,lc};
       Point(ls+5) = {""" + str(x + 1) + "," + str(y + 2.92) + """,1,lc};
       Point(ls+6) = {""" + str(x + 1) + "," + str(y + 3.92) + """,1,lc};
       Point(ls+7) = {""" + str(x) + "," + str(y + 3.92) + """,1,lc};
       Point(ls+8) = {""" + str(x + .5) + "," + str(y + .5) + """,1,lc};
       Point(ls+9) = {""" + str(x + .9) + "," + str(y + .5) + """,1,lc};
       Point(ls+10) = {""" + str(x + .9) + "," + str(y + 1) + """,1,lc};
       Point(ls+11) = {""" + str(x + .9) + "," + str(y + 1.5) + """,1,lc};
       Point(ls+12) = {""" + str(x + .5) + "," + str(y + 1.5) + """,1,lc};
       Point(ls+13) = {""" + str(x + .5) + "," + str(y + 2.42) + """,1,lc};
       Point(ls+14) = {""" + str(x + .9) + "," + str(y + 2.42) + """,1,lc};
       Point(ls+15) = {""" + str(x + .9) + "," + str(y + 2.92) + """,1,lc};
       Point(ls+16) = {""" + str(x + .9) + "," + str(y + 3.42) + """,1,lc};
       Point(ls+17) = {""" + str(x + .5) + "," + str(y + 3.42) + """,1,lc};
       Point(ls+18) = {""" + str(x) + "," + str(y) + """,0,lc};
       Point(ls+19) = {""" + str(x + 1) + "," + str(y) + """,0,lc};
       Point(ls+20) = {""" + str(x + 1) + "," + str(y + 1) + """,0,lc};
       Point(ls+21) = {""" + str(x + 1.28) + "," + str(y + 1.96) + """,0,lc};
       Point(ls+22) = {""" + str(x + 1) + "," + str(y + 2.92) + """,0,lc};
       Point(ls+23) = {""" + str(x + 1) + "," + str(y + 3.92) + """,0,lc};
       Point(ls+24) = {""" + str(x) + "," + str(y + 3.92) + """,0,lc};
       Point(ls+25) = {""" + str(x + .5) + "," + str(y + .5) + """,0,lc};
       Point(ls+26) = {""" + str(x + .9) + "," + str(y + .5) + """,0,lc};
       Point(ls+27) = {""" + str(x + .9) + "," + str(y + 1) + """,0,lc};
       Point(ls+28) = {""" + str(x + .9) + "," + str(y + 1.5) + """,0,lc};
       Point(ls+29) = {""" + str(x + .5) + "," + str(y + 1.5) + """,0,lc};
       Point(ls+30) = {""" + str(x + .5) + "," + str(y + 2.42) + """,0,lc};
       Point(ls+31) = {""" + str(x + .9) + "," + str(y + 2.42) + """,0,lc};
       Point(ls+32) = {""" + str(x + .9) + "," + str(y + 2.92) + """,0,lc};
       Point(ls+33) = {""" + str(x + .9) + "," + str(y + 3.42) + """,0,lc};
       Point(ls+34) = {""" + str(x + .5) + "," + str(y + 3.42) + """,0,lc};

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
       Ruled Surface(ls+15) = {ls+13};

           """
    for i in range(1,16):
        out += "Physical Surface(ls+"+str(i)+") = {ls+"+str(i)+"};\n"
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


def _e(x=0, y=0, l_s=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    out = "ls = " + str(l_s) + """;
       Point(ls+1) = {""" + str(x) + "," + str(y) + """,1,lc};
       Point(ls+2) = {""" + str(x + 2) + "," + str(y) + """,1,lc};
       Point(ls+3) = {""" + str(x + 2) + "," + str(y + .784) + """,1,lc};
       Point(ls+4) = {""" + str(x + .5) + "," + str(y + .784) + """,1,lc};
       Point(ls+5) = {""" + str(x + .5) + "," + str(y + 1.568) + """,1,lc};
       Point(ls+6) = {""" + str(x + 1.5) + "," + str(y + 1.568) + """,1,lc};
       Point(ls+7) = {""" + str(x + 1.5) + "," + str(y + 2.352) + """,1,lc};
       Point(ls+8) = {""" + str(x + .5) + "," + str(y + 2.352) + """,1,lc};
       Point(ls+9) = {""" + str(x + .5) + "," + str(y + 3.136) + """,1,lc};
       Point(ls+10) = {""" + str(x + 2) + "," + str(y + 3.136) + """,1,lc};
       Point(ls+11) = {""" + str(x + 2) + "," + str(y + 3.92) + """,1,lc};
       Point(ls+12) = {""" + str(x) + "," + str(y + 3.92) + """,1,lc};

       Point(ls+13) = {""" + str(x) + "," + str(y) + """,0,lc};
       Point(ls+14) = {""" + str(x + 2) + "," + str(y) + """,0,lc};
       Point(ls+15) = {""" + str(x + 2) + "," + str(y + .784) + """,0,lc};
       Point(ls+16) = {""" + str(x + .5) + "," + str(y + .784) + """,0,lc};
       Point(ls+17) = {""" + str(x + .5) + "," + str(y + 1.568) + """,0,lc};
       Point(ls+18) = {""" + str(x + 1.5) + "," + str(y + 1.568) + """,0,lc};
       Point(ls+19) = {""" + str(x + 1.5) + "," + str(y + 2.352) + """,0,lc};
       Point(ls+20) = {""" + str(x + .5) + "," + str(y + 2.352) + """,0,lc};
       Point(ls+21) = {""" + str(x + .5) + "," + str(y + 3.136) + """,0,lc};
       Point(ls+22) = {""" + str(x + 2) + "," + str(y + 3.136) + """,0,lc};
       Point(ls+23) = {""" + str(x + 2) + "," + str(y + 3.92) + """,0,lc};
       Point(ls+24) = {""" + str(x) + "," + str(y + 3.92) + """,0,lc};
       """ + _twelve

    return out


def _m(x=0, y=0, l_s=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    out = "ls = " + str(l_s) + """;
           Point(ls+1) = {""" + str(x) + "," + str(y) + """,1,lc};
           Point(ls+2) = {""" + str(x + .5) + "," + str(y) + """,1,lc};
           Point(ls+3) = {""" + str(x + .5) + "," + str(y + 3) + """,1,lc};
           Point(ls+4) = {""" + str(x + 1.5) + "," + str(y + 1) + """,1,lc};
           Point(ls+5) = {""" + str(x + 2.5) + "," + str(y + 3) + """,1,lc};
           Point(ls+6) = {""" + str(x + 2.5) + "," + str(y) + """,1,lc};
           Point(ls+7) = {""" + str(x + 3) + "," + str(y) + """,1,lc};
           Point(ls+8) = {""" + str(x + 3) + "," + str(y + 3.92) + """,1,lc};
           Point(ls+9) = {""" + str(x + 2.3) + "," + str(y + 3.92) + """,1,lc};
           Point(ls+10) = {""" + str(x + 1.5) + "," + str(y + 2.32) + """,1,lc};
           Point(ls+11) = {""" + str(x + 0.9) + "," + str(y + 3.92) + """,1,lc};
           Point(ls+12) = {""" + str(x) + "," + str(y + 3.92) + """,1,lc};

           Point(ls+13) = {""" + str(x) + "," + str(y) + """,0,lc};
           Point(ls+14) = {""" + str(x + .5) + "," + str(y) + """,0,lc};
           Point(ls+15) = {""" + str(x + .5) + "," + str(y + 3) + """,0,lc};
           Point(ls+16) = {""" + str(x + 1.5) + "," + str(y + 1) + """,0,lc};
           Point(ls+17) = {""" + str(x + 2.5) + "," + str(y + 3) + """,0,lc};
           Point(ls+18) = {""" + str(x + 2.5) + "," + str(y) + """,0,lc};
           Point(ls+19) = {""" + str(x + 3) + "," + str(y) + """,0,lc};
           Point(ls+20) = {""" + str(x + 3) + "," + str(y + 3.92) + """,0,lc};
           Point(ls+21) = {""" + str(x + 2.3) + "," + str(y + 3.92) + """,0,lc};
           Point(ls+22) = {""" + str(x + 1.5) + "," + str(y + 2.32) + """,0,lc};
           Point(ls+23) = {""" + str(x + 0.9) + "," + str(y + 3.92) + """,0,lc};
           Point(ls+24) = {""" + str(x) + "," + str(y + 3.92) + """,0,lc};
           """ + _twelve

    return out


def _plus(x=0, y=0, l_s=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    out = "ls = " + str(l_s) + """;
           Point(ls+1) = {""" + str(x) + "," + str(y + 1.71) + """,1,lc};
           Point(ls+2) = {""" + str(x + .5) + "," + str(y + 1.71) + """,1,lc};
           Point(ls+3) = {""" + str(x + .5) + "," + str(y + 1.21) + """,1,lc};
           Point(ls+4) = {""" + str(x + 1) + "," + str(y + 1.21) + """,1,lc};
           Point(ls+5) = {""" + str(x + 1) + "," + str(y + 1.71) + """,1,lc};
           Point(ls+6) = {""" + str(x + 1.5) + "," + str(y + 1.71) + """,1,lc};
           Point(ls+7) = {""" + str(x + 1.5) + "," + str(2.21) + """,1,lc};
           Point(ls+8) = {""" + str(x + 1) + "," + str(y + 2.21) + """,1,lc};
           Point(ls+9) = {""" + str(x + 1) + "," + str(y + 2.71) + """,1,lc};
           Point(ls+10) = {""" + str(x + .5) + "," + str(y + 2.71) + """,1,lc};
           Point(ls+11) = {""" + str(x + .5) + "," + str(y + 2.21) + """,1,lc};
           Point(ls+12) = {""" + str(x) + "," + str(y + 2.21) + """,1,lc};

           Point(ls+13) = {""" + str(x) + "," + str(y + 1.71) + """,0,lc};
           Point(ls+14) = {""" + str(x + .5) + "," + str(y + 1.71) + """,0,lc};
           Point(ls+15) = {""" + str(x + .5) + "," + str(y + 1.21) + """,0,lc};
           Point(ls+16) = {""" + str(x + 1) + "," + str(y + 1.21) + """,0,lc};
           Point(ls+17) = {""" + str(x + 1) + "," + str(y + 1.71) + """,0,lc};
           Point(ls+18) = {""" + str(x + 1.5) + "," + str(y + 1.71) + """,0,lc};
           Point(ls+19) = {""" + str(x + 1.5) + "," + str(2.21) + """,0,lc};
           Point(ls+20) = {""" + str(x + 1) + "," + str(y + 2.21) + """,0,lc};
           Point(ls+21) = {""" + str(x + 1) + "," + str(y + 2.71) + """,0,lc};
           Point(ls+22) = {""" + str(x + .5) + "," + str(y + 2.71) + """,0,lc};
           Point(ls+23) = {""" + str(x + .5) + "," + str(y + 2.21) + """,0,lc};
           Point(ls+24) = {""" + str(x) + "," + str(y + 2.21) + """,0,lc};
           """ + _twelve

    return out
