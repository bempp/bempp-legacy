import numpy as np
from .common import FileInterfaceImpl, Vertex, Element 

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

class GmshInterface(FileInterfaceImpl):

    def __init__(self):

        super(GmshInterface,self).__init__()
        self._version = None
        self._node_data = None
        self._element_data = None
        self._element_node_data = None

    def default_data_type(self):
        return 'element_node'

    def write(self, file_name):
        with open(file_name,'w') as f:
            self.write_version(f)
            self.write_vertices(f)
            self.write_elements(f)
            if self._node_data is not None:
                self.write_node_data(f)
            if self._element_data is not None:
                self.write_element_data(f)
            if self._element_node_data is not None:
                self.write_element_node_data(f)

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
                        gmsh_interface.vertices[vertex.index] = vertex.data
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
                            gmsh_interface.elements[element.index] = {'data':element.data, 'domain_index':element.domain_index}
                        s = f.readline().rstrip()
                        if count==number_of_elements:
                            break
                    if count!=number_of_elements:
                        raise ValueError("Expected {0} elements but got {1} elements.".format(number_of_elements,count))
                    if s!="$EndElements":
                        raise ValueError("Expected $EndElements but got {0}.".format(s))
        return gmsh_interface

    def write_vertices(self, f):
        n_vertices = len(self.vertices)
        f.write("$Nodes\n")
        f.write(str(n_vertices)+"\n")
        for key in self.vertices:
            f.write(str(key)+" "+str(self.vertices[key][0])+" "+str(self.vertices[key][1])+" "+str(self.vertices[key][2])+"\n")
        f.write("$EndNodes\n")

    def write_elements(self, f):
        n_elements = len(self.elements)
        f.write("$Elements\n")
        f.write(str(n_elements)+"\n")
        for key in self.elements:
            v0, v1, v2 = self.elements[key]['data']
            f.write(str(key)+" "+"2"+" "+"2 "+str(self.elements[key]['domain_index'])+" "+"0 "+str(v0)+" "+str(v1)+" "+str(v2)+"\n")
        f.write("$EndElements\n")

    def write_version(self, f):
        f.write("$MeshFormat\n")
        f.write("2.2 0 8\n")
        f.write("$EndMeshFormat\n")

    def write_node_data(self, f):
        label = self._node_data['label']
        data = self._node_data['data']
        f.write("$NodeData\n")
        f.write("1\n")
        f.write(label+"\n")
        f.write("1\n")
        f.write("0\n")
        f.write("4\n")
        f.write("0\n"+str(len(data.values()[0]))+"\n"+str(len(data))+"\n"+"0\n")
        for key in data:
            f.write(str(key))
            for val in data[key]:
                f.write(" "+str(val))
            f.write("\n")
        f.write("$EndNodeData\n")

    def add_node_data(self, data, label):
        self._node_data = {'label':label,'data':data}

    def write_element_data(self, f):
        label = self._element_data['label']
        data = self._element_data['data']
        f.write("$ElementData\n")
        f.write("1\n")
        f.write(label+"\n")
        f.write("1\n")
        f.write("0\n")
        f.write("4\n")
        f.write("0\n"+str(len(data.values()[0]))+"\n"+str(len(data))+"\n"+"0\n")
        for key in data:
            f.write(str(key))
            for val in data[key]:
                f.write(" "+str(val))
            f.write("\n")
        f.write("$EndElementData\n")

    def add_element_data(self, data, label):
        self._element_data = {'label':label, 'data':data}

    def write_element_node_data(self, f):
        label = self._element_node_data['label']
        data = self._element_node_data['data']
        f.write("$ElementNodeData\n")
        f.write("1\n")
        f.write(label+"\n")
        f.write("1\n")
        f.write("0\n")
        f.write("4\n")
        f.write("0\n"+str(len(data.values()[0][:,0]))+"\n"+str(len(data))+"\n"+"0\n")
        for key in data:
            f.write(str(key)+" 3")
            for i in range(3):
                for val in data[key][:,i]:
                    f.write(" "+str(val))
            f.write("\n")
        f.write("$EndElementNodeData\n")

    def add_element_node_data(self, data, label):
        self._element_node_data = {'label':label, 'data':data}





















