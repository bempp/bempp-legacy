import numpy as np
from .common import FileInterface, Vertex, Element 

def read_version(s):
    tokens = s.split()
    if len(tokens)!=3:
        raise ValueError("File header has unsupported format.")
    try:
        version = float(tokens[0])
    except:
        raise ValueError("Version number not recognized.")
    return version

def read_vertex(s):
    tokens = s.split()
    if len(tokens)!=4:
        raise ValueError("Unsupported format for vertex in string {0}".format(s))
    try:
        index = int(tokens[0])
        x = float(tokens[1])
        y = float(tokens[2])
        z = float(tokens[3])
    except:
        raise ValueError("Vertex value not recognized in string {0}".format(s))

    return Vertex(index,x,y,z)

def read_element(s):
    tokens = s.split()
    try:
        index = int(tokens[0])
        elem_type = int(tokens[1])
    except:
        raise ValueError("Unspported format for element in string {0}".format(s))
    if elem_type!=2: return None
    try:
        n_tags = int(tokens[2])
        phys_id = int(tokens[3])
        v2 = int(tokens[-1])
        v1 = int(tokens[-2])
        v0 = int(tokens[-3])
    except:
        raise ValueError("Unsupported format for element in string {0}".format(s))
    return Element(index,v0,v1,v2,phys_id)

class GmshInterface(FileInterface):

    def __init__(self):
        from collections import OrderedDict

        self._version = None
        self._vertices = {}
        self._elements = {}

    def write(self, file_name):
        pass

    @classmethod
    def read(cls, file_name):

        gmsh_interface = GmshInterface()

        with open(file_name) as f:
            while True:
                line = f.readline()
                if line=='': break
                s = line.rstrip()
                if s=="$MeshFormat":
                    s = f.readline().rstrip()
                    gmsh_interface._version = read_version(s)
                    s = f.readline().rstrip()
                    if not s=="$EndMeshFormat":
                        raise ValueError("Expected $EndMeshFormat but got {0}".format(s))
                    continue
                if s=="$Nodes":
                    s = f.readline().rstrip()
                    try:
                        number_of_vertices = int(s)
                    except:
                        raise ValueError("Expected integer, got {0}".format(s))

                    count = 0
                    s = f.readline().rstrip()
                    while s!="$EndNodes":
                        vertex = read_vertex(s)
                        gmsh_interface._vertices[vertex.index] = vertex.data
                        count += 1
                        s = f.readline().rstrip()
                        if count==number_of_vertices:
                            break
                    if count!=number_of_vertices:
                        raise ValueError("Expected {0} vertices but got {1} vertices.".format(number_of_vertices,count))
                    if s!="$EndNodes":
                        raise ValueError("Expected $EndNodes but got {0}.".format(s))
                if s=="$Elements":
                    s = f.readline().rstrip()
                    try:
                        number_of_elements = int(s)
                    except:
                        raise ValueError("Expected integer, got {0}".format(s))
                    count = 0
                    s = f.readline().rstrip()
                    while s!= "$EndElements":
                        element = read_element(s)
                        count += 1
                        if element is not None:
                            gmsh_interface._elements[element.index] = {'data':element.data, 'domain_index':element.domain_index}
                        s = f.readline().rstrip()
                        if count==number_of_elements:
                            break
                    if count!=number_of_elements:
                        raise ValueError("Expected {0} elements but got {1} elements.".format(number_of_elements,count))
                    if s!="$EndElements":
                        raise ValueError("Expected $EndElements but got {0}.".format(s))
        return gmsh_interface


