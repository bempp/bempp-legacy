%{
#include "io/gmsh.hpp"
%}


%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (std::vector< std::vector< double > >& ARGOUT_VECVEC_VALUES)
    (PyArrayObject* array=NULL, std::vector< std::vector< double > > std_vecvec)
{
    $1 = &std_vecvec;
}
%typemap(argout)
    (std::vector< std::vector< double > >& ARGOUT_VECVEC_VALUES)
{

    int listSize$argnum = std_vecvec$argnum.size();
    PyObject* outerList$argnum = PyList_New(listSize$argnum);
    for (int i = 0; i < listSize$argnum; ++i ){
        npy_intp dims[1];
        dims[0] = std_vecvec$argnum[i].size();
        array$argnum = reinterpret_cast<PyArrayObject*>(
            PyArray_EMPTY(1, dims, NPY_DOUBLE, NPY_FORTRAN));
        if (!array$argnum)
            SWIG_fail;
        std::copy(std_vecvec$argnum[i].begin(), std_vecvec$argnum[i].end(),
            reinterpret_cast<double*>(array_data(array$argnum)));
        PyList_SET_ITEM(outerList$argnum,i,reinterpret_cast<PyObject*>(array$argnum));
    }
        $result = SWIG_Python_AppendOutput($result,outerList$argnum);
}

%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (std::vector< std::string >& ARGOUT_STRING)
    (PyArrayObject* array=NULL, std::vector< std::string > std_vec)
{
    $1 = &std_vec;
}
%typemap(argout)
    (std::vector< std::string >& ARGOUT_STRING)
{

    int listSize$argnum = std_vec$argnum.size();
    PyObject* outerList$argnum = PyList_New(listSize$argnum);
    for (int i = 0; i < listSize$argnum; ++i )
        PyList_SET_ITEM(outerList$argnum,i,PyString_FromString(std_vec$argnum[i].c_str()));
    $result = SWIG_Python_AppendOutput($result,outerList$argnum);
}

%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (std::vector< std::vector< std::vector<double> > >& ARGOUT_VECVECVEC_VALUES)
    (PyArrayObject* array=NULL, std::vector< std::vector< std::vector<double> > > std_vecvecvec)
{
    $1 = &std_vecvecvec;
}
%typemap(argout)
    (std::vector< std::vector< std::vector<double> > >& ARGOUT_VECVECVEC_VALUES)
{

    int listSize$argnum = std_vecvecvec$argnum.size();
    PyObject* elementList$argnum = PyList_New(listSize$argnum);
    for (int i = 0; i < listSize$argnum; ++i ){
        PyObject* nodeList = PyList_New(std_vecvecvec$argnum[i].size());
        for (int j = 0; j < std_vecvecvec$argnum[i].size(); ++j) {
            npy_intp dims[1];
            dims[0] = std_vecvecvec$argnum[i][j].size();
            array$argnum = reinterpret_cast<PyArrayObject*>(
            PyArray_EMPTY(1, dims, NPY_DOUBLE, NPY_FORTRAN));
            if (!array$argnum)
                SWIG_fail;
            std::copy(std_vecvecvec$argnum[i][j].begin(), std_vecvecvec$argnum[i][j].end(),
                reinterpret_cast<double*>(array_data(array$argnum)));
            PyList_SET_ITEM(nodeList,j,reinterpret_cast<PyObject*>(array$argnum));
         }
         PyList_SET_ITEM(elementList$argnum,i,nodeList);
    }
    $result = SWIG_Python_AppendOutput($result,elementList$argnum);
}


namespace Bempp {

class GmshData {

public:

    GmshData();
    int numberOfNodes() const;
    int numberOfElements() const;
    int numberOfPeriodicEntities() const;
    int numberOfPeriodicNodes() const;
    int numberOfPhysicalNames() const;
    int numberOfNodeDataSets() const;
    int numberOfElementDataSets() const;
    int numberOfElementNodeDataSets() const;
    int numberOfInterpolationSchemeSets() const;

    void addNode(int index, double x, double y, double z);

    void addElement(int index, int elementType, const std::vector<int>& nodes,
                    int physicalEntity, int elementaryEntity,
                    const std::vector<int>& partitions = std::vector<int>());

    void addPeriodicEntity(int dimension, int slaveEntityTag, int masterEntityTag);
    void addPeriodicNode(int slaveNode, int masterNode);

    void addPhysicalName(int dimension, int number, std::string name);

    void addNodeDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags,
                       int numberOfFieldComponents, int capacity=0, int timeStep=0, int partition=0);
    void addElementDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags,
                        int numberOfFieldComponents, int capacity=0, int timeStep=0, int partition=0);
    void addElementNodeDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags,
                        int numberOfFieldComponents,
                        int capacity=0, int timeStep=0, int partition=0);

    void addNodeData(int dataSetIndex, int node, const std::vector<double>& values);
    void addElementData(int dataSetIndex, int element, const std::vector<double>& values);
    void addElementNodeData(int dataSetIndex, int element, const std::vector<std::vector<double> >& values);

    void addInterpolationSchemeSet(std::string name, int topology);
    void addInterpolationMatrix(int dataSetIndex, int nrows, int ncols, const std::vector<double>& values);

    %apply double& OUTPUT { double& x };
    %apply double& OUTPUT { double& y };
    %apply double& OUTPUT { double& z };
    %apply int& OUTPUT { int& physicalEntity };
    %apply int& OUTPUT { int& elementaryEntity };
    %apply int& OUTPUT { int& elementType };
    %apply int& OUTPUT { int& slaveEntityTag };
    %apply int& OUTPUT { int& masterEntityTag };
    %apply int& OUTPUT { int& slaveNode };
    %apply int& OUTPUT { int& masterNode };
    %apply int& OUTPUT { int& dimension };
    %apply int& OUTPUT { int& timeStep };
    %apply int& OUTPUT { int& partition };
    %apply int& OUTPUT { int& numberOfFieldComponents };
    %apply int& OUTPUT { int& topology };

    %apply std::vector<int>& ARGOUT_VEC { std::vector<int>& indices };
    %apply std::vector<int>& ARGOUT_VEC { std::vector<int>& nodeIndices };
    %apply std::vector<int>& ARGOUT_VEC { std::vector<int>& elementIndices };
    %apply std::vector<int>& ARGOUT_VEC { std::vector<int>& nodes };
    %apply std::vector<int>& ARGOUT_VEC { std::vector<int>& partitions };
    %apply std::vector<double>& ARGOUT_VEC { std::vector<double>& realTags };
    %apply std::vector<std::string>& ARGOUT_STRING { std::vector<std::string>& stringTags };
    %apply std::vector<std::vector<double> >& ARGOUT_VECVEC_VALUES { std::vector<std::vector<double> >& values };
    %apply std::vector<std::vector<std::vector<double> > >& ARGOUT_VECVECVEC_VALUES { std::vector<std::vector<std::vector<double> > >& values };

    %pythonappend getNodeDataSet %{
        val = {'stringTags':val[0],'realTags':val[1],'numberOfFieldComponents':val[2],'indices':val[3],'values':val[4],'timeStep':val[5],'partition':val[6]}
    %}
    %pythonappend getElementDataSet %{
        val = {'stringTags':val[0],'realTags':val[1],'numberOfFieldComponents':val[2],'indices':val[3],'values':val[4],'timeStep':val[5],'partition':val[6]}
    %}
    %pythonappend getElementNodeDataSet %{
        val = {'stringTags':val[0],'realTags':val[1],'numberOfFieldComponents':val[2],'indices':val[3],'values':val[4],'timeStep':val[5],'partition':val[6]}
    %}

    %pythonappend getElement %{
        val = {'index':args[0],'type':val[0],'nodes':val[1],'physicalEntity':val[2],'elementaryEntity':val[3]}
    %}

    %pythonappend getNode %{
        import numpy as np
        val = {'index':args[0],'coordinates':np.array(val,dtype='float64')}

    %}

    %pythoncode %{

        @property
        def elements(self):
            """Iterator for mesh elements"""

            indices = self.getElementIndices()
            numberOfIndices = len(indices)
            index = 0
            while index<numberOfIndices:
                yield self.getElement(indices[index])
                index += 1

        @property
        def nodes(self):
            """Iterator for mesh nodes"""

            indices = self.getNodeIndices()
            numberOfIndices = len(indices)
            index = 0
            while index<numberOfIndices:
                yield self.getNode(indices[index])
                index += 1

        @property
        def nodeDataSets(self):
            """Iterator for node datasets"""

            numberOfDataSets = self.numberOfNodeDataSets
            index = 0
            while index<numberOfDataSets:
                yield self.getNodeDataSet(index)
                index += 1

        @property
        def elementDataSets(self):
            """Iterator for element datasets"""

            numberOfDataSets = self.numberOfElementDataSets
            index = 0
            while index<numberOfDataSets:
                yield self.getElementDataSet(index)
                index += 1

        @property
        def elementNodeDataSets(self):
            """Iterator for elementNode datasets"""

            numberOfDataSets = self.numberOfElementNodeDataSets
            index = 0
            while index<numberOfDataSets:
                yield self.getElementNodeDataSet(index)
                index += 1


    %}


    void getNodeIndices(std::vector<int>& indices) const;
    void getElementIndices(std::vector<int>& indices) const;

    void getNode(int index, double &x, double &y, double &z) const;
    void getElement(int index, int& elementType, std::vector<int>& nodes,
                    int& physicalEntity, int& elementaryEntity, std::vector<int>& partitions) const;
    void getPeriodicEntity(int index, int& dimension, int& slaveEntityTag, int& masterEntityTag) const;
    void getPeriodicNode(int index, int& slaveNode, int& masterNode) const;
    void getPhysicalName(int index, int& dimension, int& number, std::string& name) const;
    void getNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags,
                        int& numberOfFieldComponents,
                        std::vector<int>& nodeIndices, std::vector<std::vector<double> >& values,
                        int& timeStep, int& partition) const;
    void getElementDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags,
                        int& numberOfFieldComponents,
                        std::vector<int>& elementIndices, std::vector<std::vector<double> >& values,
                        int& timeStep, int& partition) const;
    void getElementNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags,
                        int& numberOfFieldComponents,
                        std::vector<int>& elementIndices, std::vector<std::vector<std::vector<double> > >& values,
                        int& timeStep, int& partition) const;
    void getInterpolationSchemeSet(int index, std::string& name, int& topology, std::vector<int> &nrows,
                                   std::vector<int> &ncols,
                                   std::vector<std::vector<double> >& values);

    void reserveNumberOfNodes(int n);
    void reserveNumberOfElements(int n);

    void write(const std::string& fname) const;
    static GmshData read(const std::string& fname);


    %clear double& x;
    %clear double& y;
    %clear double& z;
    %clear int& physicalEntity;
    %clear int& elementaryEntity;
    %clear int& elementType;
    %clear int& slaveEntityTag;
    %clear int& masterEntityTag;
    %clear int& slaveNode;
    %clear int& masterNode;
    %clear int& dimension;
    %clear int& timeStep;
    %clear int& partition;
    %clear int& numberOfFieldComponents;
    %clear int& topology;

    %clear std::vector<int>& indices;
    %clear std::vector<int>& nodeIndices;
    %clear std::vector<int>& elementIndices;
    %clear std::vector<int>& nodes;
    %clear std::vector<int>& partitions;
    %clear std::vector<double>& realTags;
    %clear std::vector<std::string>& stringTags;
    %clear std::vector<std::vector<double> >& values;
    %clear std::vector<std::vector<std::vector<double> > >& values;



};

}

namespace Bempp
{

class GmshIo {

public:

    GmshIo();
    GmshIo(GmshData gmshData);
    GmshIo(std::string fname);

    boost::shared_ptr<Grid> grid() const;
    const std::vector<int>& nodePermutation() const;
    const std::vector<int>& elementPermutation() const;
    const std::vector<int>& inverseNodePermutation() const;
    const std::vector<int>& inverseElementPermutation() const;

};

}


