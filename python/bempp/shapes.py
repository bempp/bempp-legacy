# Copyright (C) 2011-2012 by the BEM++ Authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np


gmsh_exe = None

# Solution from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
# to find executable in Path

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def getGmshFile():
    """
    Return a 3-tuple (geo_file,geo_name,msh_name), where
    geo_file is a file descriptor to an empty .geo file, geo_name is
    the corresponding filename and msh_name is the name of the
    Gmsh .msh file that will be generated.

    """
    import os, tempfile
    geo, geo_name = tempfile.mkstemp(suffix='.geo',dir=os.getcwd(), text=True)
    geo_file = os.fdopen(geo,"w")
    msh_name = os.path.splitext(geo_name)[0]+".msh"
    return (geo_file,geo_name,msh_name)

def findGmsh(fpath=None):
    """
    Find the Gmsh executable and return the corresponding path.
    Optionally, the full path to the Gmsh executable can be given
    """

    import sys

    gmshPath = None
    if sys.platform.startswith('linux2'):
        exe_name = 'gmsh'
    elif sys.platform.startswith('darwin'):
        exe_name = '/Applications/Gmsh.app/Contents/MacOS/gmsh'
    else:
        raise Exception("findGmsh: Platform not supported")

    if fpath is not None:
        gmshPath = which(fpath)
        if gmshPath is None:
            gmshPath = which(exe_name)
    else:
        gmshPath = which(exe_name)

    if gmshPath is not None:
        gmsh_exe = gmshPath
        print "Found Gmsh at "+gmsh_exe
    else:
        raise Exception("findGmsh: Could not find Gmsh. Try to set bempp.shapes.gmsh_exe directly to an existing executable.")
    return gmsh_exe

def __generate_grid_from_string(geo_string,grid=True,msh_file=False,parallel=True):
    """Helper routine that implements the grid generation
    """

    import os
    import subprocess

    def msh_from_string(geo_string):
        gmsh_command = findGmsh(gmsh_exe)
        f,geo_name,msh_name = getGmshFile()
        f.write(geo_string)
        f.close()

        fnull = open(os.devnull,'w')
        cmd = gmsh_command+" -2 "+geo_name
        try:
            print "Generating Gmsh grid..."
            subprocess.check_call(cmd,shell=True,stdout=fnull,stderr=fnull)
        except:
            print "The following command failed: "+cmd
            fnull.close()
            raise
        os.remove(geo_name)
        fnull.close()
        return msh_name
        
    if parallel:
        from PyTrilinos.Epetra import PyComm
        comm = PyComm()
        my_rank = comm.MyPID()
    else:
        my_rank = 0
        

    if parallel:
        if my_rank == 0:
            rank0_msh_name = msh_from_string(geo_string)
            nchar = np.array(len(rank0_msh_name),dtype='int')
            comm.Broadcast(nchar,0)
            msh_name_array=np.array(rank0_msh_name)
            comm.Broadcast(msh_name_array,0)
        else:
            nchar=np.array([0],dtype='int')
            comm.Broadcast(nchar,0)
            msh_name_array=''
            msh_name_array=np.array(msh_name_array.zfill(nchar))
            comm.Broadcast(msh_name_array,0)
        msh_name = msh_name_array.item()
    else:
        msh_name = msh_from_string(geo_string)
        
    grid_obj = None
    if grid:
        from bempp.lib import createGridFactory
        grid_obj = createGridFactory().importGmshGrid("triangular",msh_name)
        comm.Barrier()
    if not msh_file:
        if parallel:
            if my_rank==0:
                os.remove(msh_name)
            msh_name = None
        else:
            os.remove(msh_name)
            msh_name = None            
    if grid and not msh_name:
        return grid_obj
    elif msh_file and not grid:
        return msh_name
    elif msh_file and grid:
        return (grid_obj,msh_file)
    else:
        return None


def sphere(radius=1,origin=(0,0,0),h=0.1,grid=True,msh_file=False,parallel=True):
    """
    Return a shpere grid.

    If 'grid=True' and 'msh_file=False' a Bempp grid object is returned. If 'grid=False' and 'msh_file=True'
    a string to a .msh file containing a sphere mesh is returned. If both are true a tuple (sphere,fname) 
    is returned. If this function is run in an MPI context it guarantees that each process returns the same
    grid. Note that all processes need to have access to the same filesystem.

    *Parameters:*
       - radius (real number)
            Radius of the sphere.
       - origin (3-tuple)
            Origin of the sphere.
       - h (real number)
            Approximate element size.
       - grid (True/False)
            If true return an assembled grid object.
       - msh_file (True/False)
            If true return the Gmsh .msh file.
            This file is deleted if msh_file=False.
       - parallel (True/False)
            Run the function as parallelized MPI routine so that
            all processes return the same grid
    """
    import subprocess,os

    sphere_stub = """
    Point(1) = {orig0,orig1,orig2,cl};
    Point(2) = {orig0+rad,orig1,orig2,cl};
    Point(3) = {orig0,orig1+rad,orig2,cl};
    Circle(1) = {2,1,3};
    Point(4) = {orig0-rad,orig1,orig2,cl};
    Point(5) = {orig0,orig1-rad,orig2,cl};
    Circle(2) = {3,1,4};
    Circle(3) = {4,1,5};
    Circle(4) = {5,1,2};
    Point(6) = {orig0,orig1,orig2-rad,cl};
    Point(7) = {orig0,orig1,orig2+rad,cl};
    Circle(5) = {3,1,6};
    Circle(6) = {6,1,5};
    Circle(7) = {5,1,7};
    Circle(8) = {7,1,3};
    Circle(9) = {2,1,7};
    Circle(10) = {7,1,4};
    Circle(11) = {4,1,6};
    Circle(12) = {6,1,2};
    Line Loop(13) = {2,8,-10};
    Ruled Surface(14) = {13};
    Line Loop(15) = {10,3,7};
    Ruled Surface(16) = {15};
    Line Loop(17) = {-8,-9,1};
    Ruled Surface(18) = {17};
    Line Loop(19) = {-11,-2,5};
    Ruled Surface(20) = {19};
    Line Loop(21) = {-5,-12,-1};
    Ruled Surface(22) = {21};
    Line Loop(23) = {-3,11,6};
    Ruled Surface(24) = {23};
    Line Loop(25) = {-7,4,9};
    Ruled Surface(26) = {25};
    Line Loop(27) = {-4,12,-6};
    Ruled Surface(28) = {27};
    Surface Loop(29) = {28,26,16,14,20,24,22,18};
    Volume(30) = {29};
    Mesh.Algorithm = 6;
    """

    sphere_geometry = ("rad = "+str(radius)+";\n"+
            "orig0 = "+str(origin[0])+";\n"+
            "orig1 = "+str(origin[1])+";\n"+
            "orig2 = "+str(origin[2])+";\n"+
            "cl = "+str(h)+";\n"+sphere_stub)

    return __generate_grid_from_string(sphere_geometry,grid,msh_file,parallel)

        
def cube(length=1,origin=(0,0,0),h=0.1,grid=True,msh_file=False):
    """
    Return a shpere grid.

    If 'grid=True' and 'msh_file=False' a Bempp grid object is returned. If 'grid=False' and 'msh_file=True'
    a string to a .msh file containing a cube mesh is returned. If both are true a tuple (cube,fname) 
    is returned. If this function is run in an MPI context it guarantees that each process returns the same
    grid. Note that all processes need to have access to the same filesystem.

    *Parameters:*
       - length (real number)
            Length of the cube.
       - origin (3-tuple)
            Origin of the cube.
       - h (real number)
            Approximate element size.
       - grid (True/False)
            If true return an assembled grid object.
       - msh_file (True/False)
            If true return the Gmsh .msh file.
            This file is deleted if msh_file=False.
       - parallel (True/False)
            Run the function as parallelized MPI routine so that
            all processes return the same grid

    """
    import subprocess,os

    cube_stub = """
    Point(1) = {orig0,orig1,orig2,cl};
    Point(2) = {orig0+l,orig1,orig2,cl};
    Point(3) = {orig0+l,orig1+l,orig2,cl};
    Point(4) = {orig0,orig1+l,orig2,cl};
    Point(5) = {orig0,orig1,orig2+l,cl};
    Point(6) = {orig0+l,orig1,orig2+l,cl};
    Point(7) = {orig0+l,orig1+l,orig2+l,cl};
    Point(8) = {orig0,orig1+l,orig2+l,cl};

    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Line(5) = {1,5};
    Line(6) = {2,6};
    Line(7) = {3,7};
    Line(8) = {4,8};
    Line(9) = {5,6};
    Line(10) = {6,7};
    Line(11) = {7,8};
    Line(12) = {8,5};

    Line Loop(1) = {-1,-4,-3,-2};
    Line Loop(2) = {1,6,-9,-5};
    Line Loop(3) = {2,7,-10,-6};
    Line Loop(4) = {3,8,-11,-7};
    Line Loop(5) = {4,5,-12,-8};
    Line Loop(6) = {9,10,11,12};

    Plane Surface(1) = {1};
    Plane Surface(2) = {2};
    Plane Surface(3) = {3};
    Plane Surface(4) = {4};
    Plane Surface(5) = {5};
    Plane Surface(6) = {6};

    Surface Loop (1) = {1,2,3,4,5,6};

    Volume (1) = {1};

    Mesh.Algorithm = 6;
    """

    cube_geometry = ("l = "+str(length)+";\n"+
            "orig0 = "+str(origin[0])+";\n"+
            "orig1 = "+str(origin[1])+";\n"+
            "orig2 = "+str(origin[2])+";\n"+
            "cl = "+str(h)+";\n"+cube_stub)

    return __generate_grid_from_string(cube_geometry,grid,msh_file,parallel)
