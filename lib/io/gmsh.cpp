// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "gmsh.hpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"
#include "../space/piecewise_constant_scalar_space.hpp"
#include "../space/piecewise_linear_continuous_scalar_space.hpp"
#include "../space/piecewise_linear_discontinuous_scalar_space.hpp"
#include "../common/scalar_traits.hpp"
#include <boost/array.hpp>
#include <armadillo>

#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/entity.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_factory.hpp"
#include "../grid/grid_parameters.hpp"
#include "../grid/geometry.hpp"
#include "../common/complex_aux.hpp"
#include "../common/acc.hpp"






namespace {

typedef std::vector<std::string> StringVector;

StringVector stringTokens(std::string input) {

    StringVector result;
    boost::char_delimiters_separator<char> separator(false," ");
    boost::tokenizer<boost::char_delimiters_separator<char> > tok(input,separator);
    boost::tokenizer<>::iterator it;
    for (it = tok.begin(); it != tok.end(); ++it)
        result.push_back(*it);
    return result;
}

} // namespace

namespace Bempp {

GmshData::GmshData() :
    m_fileType(0),
    m_dataSize(0),
    m_numberOfNodes(0),
    m_numberOfElements(0) {}



int GmshData::numberOfNodes() const {

    return m_numberOfNodes;

}
int GmshData::numberOfElements() const {

    return m_numberOfElements;

}

int GmshData::numberOfPeriodicEntities() const {

    return m_periodicEntities.size();

}

int GmshData::numberOfPeriodicNodes() const {

    return m_periodicNodes.size();
}

int GmshData::numberOfPhysicalNames() const {

    return m_physicalNames.size();

}

int GmshData::numberOfNodeDataSets() const {

    return m_nodeDataSets.size();

}

int GmshData::numberOfElementDataSets() const {

    return m_elementDataSets.size();
}


int GmshData::numberOfElementNodeDataSets() const {

    return m_elementNodeDataSets.size();

}
int GmshData::numberOfInterpolationSchemeSets() const {

    return m_interpolationSchemeSets.size();

}

void GmshData::addNode(int index, double x, double y, double z){

    if (index>=m_nodes.size()) m_nodes.resize(index+1);
    if (!m_nodes[index]) {
        m_nodes[index] = shared_ptr<Node>(new Node());
        ++m_numberOfNodes;
    }

    m_nodes[index]->x = x;
    m_nodes[index]->y = y;
    m_nodes[index]->z = z;

}

void GmshData::addElement(int index, int elementType, const std::vector<int>& nodes,
                int physicalEntity, int elementaryEntity,
                const std::vector<int>& partitions){

    if (index>=m_elements.size()) m_elements.resize(index+1);
    if (!m_elements[index]) {
        m_elements[index] = shared_ptr<Element>(new Element());
        ++m_numberOfElements;
    }

    m_elements[index]->type = elementType;
    m_elements[index]->physicalEntity = physicalEntity;
    m_elements[index]->elementaryEntity = elementaryEntity;
    m_elements[index]->nodes = nodes;
    m_elements[index]->partitions = partitions;


}

void GmshData::addPeriodicEntity(int dimension, int slaveEntityTag, int masterEntityTag){

    PeriodicEntity e;
    e.dimension = dimension;
    e.slaveTag = slaveEntityTag;
    e.masterTag = masterEntityTag;

    m_periodicEntities.push_back(e);

}

void GmshData::addPeriodicNode(int slaveNode, int masterNode){

    PeriodicNode pn;
    pn.masterNode = masterNode;
    pn.slaveNode = slaveNode;
    m_periodicNodes.push_back(pn);

}

void GmshData::addPhysicalName(int dimension, int number, std::string name) {

    PhysicalName p;
    p.dimension = dimension;
    p.number = number;
    p.name = name;
    m_physicalNames.push_back(p);


}

void GmshData::addNodeDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags, int numberOfFieldComponents,
                    int capacity, int timeStep, int partition) {

    m_nodeDataSets.push_back(shared_ptr<NodeDataSet>(new NodeDataSet()));
    NodeDataSet& data = *m_nodeDataSets.back();

    data.stringTags = stringTags;
    data.realTags = realTags;
    data.timeStep = timeStep;
    data.numberOfFieldComponents = numberOfFieldComponents;
    data.partition = partition;
    data.nodeIndices.reserve(capacity);
    data.values.reserve(capacity);

}
void GmshData::addElementDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags,
                    int numberOfFieldComponents, int capacity, int timeStep, int partition) {

    m_elementDataSets.push_back(shared_ptr<ElementDataSet>(new ElementDataSet()));
    ElementDataSet& data = *m_elementDataSets.back();

    data.stringTags = stringTags;
    data.realTags = realTags;
    data.timeStep = timeStep;
    data.numberOfFieldComponents = numberOfFieldComponents;
    data.partition = partition;
    data.elementIndices.reserve(capacity);
    data.values.reserve(capacity);


}
void GmshData::addElementNodeDataSet(const std::vector<std::string>& stringTags, const std::vector<double>& realTags, int numberOfFieldComponents, int capacity,
                    int timeStep,
                    int partition) {

    m_elementNodeDataSets.push_back(shared_ptr<ElementNodeDataSet>(new ElementNodeDataSet()));
    ElementNodeDataSet& data = *m_elementNodeDataSets.back();

    data.stringTags = stringTags;
    data.realTags = realTags;
    data.timeStep = timeStep;
    data.numberOfFieldComponents = numberOfFieldComponents;
    data.partition = partition;
    data.elementIndices.reserve(capacity);
    data.values.reserve(capacity);

}

void GmshData::addNodeData(int dataSetIndex, int node, const std::vector<double>& values) {

    NodeDataSet& nodeDataSet = *m_nodeDataSets.at(dataSetIndex);
    nodeDataSet.nodeIndices.push_back(node);
    nodeDataSet.values.push_back(values);
}
void GmshData::addElementData(int dataSetIndex, int element, const std::vector<double>& values) {

    ElementDataSet& elementDataSet = *m_elementDataSets.at(dataSetIndex);
    elementDataSet.elementIndices.push_back(element);
    elementDataSet.values.push_back(values);

}
void GmshData::addElementNodeData(int dataSetIndex, int element, const std::vector<std::vector<double> >& values) {

    ElementNodeDataSet& elementNodeDataSet = *m_elementNodeDataSets.at(dataSetIndex);
    elementNodeDataSet.elementIndices.push_back(element);
    elementNodeDataSet.values.push_back(values);

}

void GmshData::addInterpolationSchemeSet(std::string name, int topology) {

    m_interpolationSchemeSets.push_back(shared_ptr<InterpolationSchemeSet>(new InterpolationSchemeSet()));
    InterpolationSchemeSet& scheme = *m_interpolationSchemeSets.back();
    scheme.name = name;
    scheme.topology = topology;
}

void GmshData::addInterpolationMatrix(int dataSetIndex, int nrows, int ncols, const std::vector<double>& values) {

    InterpolationSchemeSet& scheme = *m_interpolationSchemeSets.at(dataSetIndex);
    scheme.ncols.push_back(ncols);
    scheme.nrows.push_back(nrows);
    scheme.values.push_back(values);

}

void GmshData::getNodeIndices(std::vector<int>& indices) const {

    indices.clear();
    indices.reserve(m_numberOfNodes);
    for (int i = 0; i < m_nodes.size(); i++ )
        if (m_nodes[i]) indices.push_back(i);


}

void GmshData::getElementIndices(std::vector<int>& indices) const {

    indices.clear();
    indices.reserve(m_numberOfElements);
    for (int i = 0; i < m_elements.size(); i++ )
        if (m_elements[i]) indices.push_back(i);


}


void GmshData::getNode(int index, double& x, double& y, double& z) const {

    if (index >= m_nodes.size()) throw std::runtime_error(
                "GmshData::getNode(): Index does not exist.");
    if (m_nodes[index]) {
        x = m_nodes[index]->x;
        y = m_nodes[index]->y;
        z = m_nodes[index]->z;
    }
    else throw std::runtime_error(
                "GmshData::getNode(): Index does not exist.");

}

void GmshData::getElement(int index, int& elementType, std::vector<int>& nodes,
                int& physicalEntity, int& elementaryEntity, std::vector<int>& partitions) const {

    if (index >= m_elements.size()) throw std::runtime_error(
                "GmshData::getElement(): Index does not exist.");

    if (m_elements[index]) {
        elementType = m_elements[index]->type;
        nodes = m_elements[index]->nodes;
        physicalEntity = m_elements[index]->physicalEntity;
        elementaryEntity = m_elements[index]->elementaryEntity;
        partitions = m_elements[index]->partitions;
    }
    else throw std::runtime_error(
                "GmshData::getElement(): Index does not exist.");

}
void GmshData::getElement(int index, int& elementType, std::vector<int>& nodes,
                int& physicalEntity, int& elementaryEntity) const {

    std::vector<int> partitions;
    getElement(index, elementType, nodes, physicalEntity, elementaryEntity,partitions);

}
void GmshData::getPeriodicEntity(int index, int& dimension, int& slaveEntityTag, int& masterEntityTag) const {

    if (index<m_periodicEntities.size()) {

        dimension = m_periodicEntities[index].dimension;
        slaveEntityTag = m_periodicEntities[index].slaveTag;
        masterEntityTag = m_periodicEntities[index].masterTag;
    }
    else throw std::runtime_error(
                "GmshData::getPeriodicEntity(): Index does not exist.");

}

void GmshData::getPeriodicNode(int index, int& slaveNode, int& masterNode) const {

    if (index < m_periodicNodes.size()) {

        slaveNode = m_periodicNodes[index].slaveNode;
        masterNode = m_periodicNodes[index].masterNode;
    }
    else throw std::runtime_error(
                "GmshData::getPeriodicNode(): Index does not exist.");


}

void GmshData::getPhysicalName(int index, int& dimension, int& number, std::string& name) const {

    if (index < m_physicalNames.size()) {
        dimension = m_physicalNames[index].dimension;
        number = m_physicalNames[index].number;
        name = m_physicalNames[index].name;
    }
    else throw std::runtime_error(
                "GmshData::getPhysicalName(): Index does not exist");

}

void GmshData::getNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& nodeIndices, std::vector<std::vector<double> >& values,
                    int& timeStep,
                    int& partition) const {

    if (index< m_nodeDataSets.size()) {
        stringTags = m_nodeDataSets[index]->stringTags;
        realTags = m_nodeDataSets[index]->realTags;
        numberOfFieldComponents = m_nodeDataSets[index]->numberOfFieldComponents;
        nodeIndices = m_nodeDataSets[index]->nodeIndices;
        values = m_nodeDataSets[index]->values;
        timeStep = m_nodeDataSets[index]->timeStep;
        partition = m_nodeDataSets[index]->partition;
    }
    else throw std::runtime_error(
                "Gmsh::getNodeDataSet(): Index does not exist.");

}
void GmshData::getNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& nodeIndices, std::vector<std::vector<double> >& values) const {

    int timeStep;
    int partition;
    getNodeDataSet(index,stringTags,realTags,numberOfFieldComponents,nodeIndices,values,timeStep,partition);

}
void GmshData::getElementDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& elementIndices, std::vector<std::vector<double> >& values,
                    int& timeStep,
                    int& partition) const {

    if (index< m_elementDataSets.size()) {
        stringTags = m_elementDataSets[index]->stringTags;
        realTags = m_elementDataSets[index]->realTags;
        numberOfFieldComponents = m_elementDataSets[index]->numberOfFieldComponents;
        elementIndices = m_elementDataSets[index]->elementIndices;
        values = m_elementDataSets[index]->values;
        timeStep = m_elementDataSets[index]->timeStep;
        partition = m_elementDataSets[index]->partition;
    }
    else throw std::runtime_error(
                "Gmsh::getElementDataSet(): Index does not exist.");



}
void GmshData::getElementDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& elementIndices, std::vector<std::vector<double> >& values) const {

    int timeStep;
    int partition;
    getElementDataSet(index,stringTags,realTags,numberOfFieldComponents,elementIndices,values,timeStep,partition);

}
void GmshData::getElementNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& elementIndices, std::vector<std::vector<std::vector<double> > >& values,
                    int& timeStep,
                    int& partition) const {

    if (index< m_elementNodeDataSets.size()) {
        stringTags = m_elementNodeDataSets[index]->stringTags;
        realTags = m_elementNodeDataSets[index]->realTags;
        numberOfFieldComponents = m_elementNodeDataSets[index]->numberOfFieldComponents;
        elementIndices = m_elementNodeDataSets[index]->elementIndices;
        values = m_elementNodeDataSets[index]->values;
        timeStep = m_elementNodeDataSets[index]->timeStep;
        partition = m_elementNodeDataSets[index]->partition;
    }
    else throw std::runtime_error(
                "Gmsh::getElementNodeDataSet(): Index does not exist.");


}
void GmshData::getElementNodeDataSet(int index, std::vector<std::string>& stringTags, std::vector<double>& realTags, int& numberOfFieldComponents,
                    std::vector<int>& elementIndices, std::vector<std::vector<std::vector<double > > >& values) const {

    int timeStep;
    int partition;
    getElementNodeDataSet(index,stringTags,realTags,numberOfFieldComponents,elementIndices,
                          values,timeStep,partition);

}
void GmshData::getInterpolationSchemeSet(int index, std::string& name, int& topology, std::vector<int>& nrows,
                               std::vector<int>& ncols,
                               std::vector<std::vector<double> >& values) {

    if (index < m_interpolationSchemeSets.size()) {

        name = m_interpolationSchemeSets[index]->name;
        topology = m_interpolationSchemeSets[index]->topology;
        nrows = m_interpolationSchemeSets[index]->nrows;
        ncols = m_interpolationSchemeSets[index]->ncols;
        values = m_interpolationSchemeSets[index]->values;

    }
    else throw std::runtime_error(
                "Gmsh::getInterpolationSchemeSet(): Index does not exist.");

}

void GmshData::reserveNumberOfNodes(int n) {

    m_nodes.reserve(n+1);

}
void GmshData::reserveNumberOfElements(int n) {

    m_elements.reserve(n+1);
}

void GmshData::write(std::ostream& output) const {


    output << "$MeshFormat" << std::endl;
    output << "2.2" << " " <<
              0 << " " <<
              sizeof(double) << std::endl;
    output << "$EndMeshFormat" << std::endl;

    if (m_numberOfNodes >0){
        std::vector<int> nodeIndices;
        getNodeIndices(nodeIndices);

        output << "$Nodes" << std::endl;
        output << m_numberOfNodes << std::endl;
        for (int i=0; i < m_numberOfNodes;i++){
            output << nodeIndices[i] << " " <<
            boost::lexical_cast<std::string>(m_nodes[nodeIndices[i]]->x) << " " <<
            boost::lexical_cast<std::string>(m_nodes[nodeIndices[i]]->y) << " " <<
            boost::lexical_cast<std::string>(m_nodes[nodeIndices[i]]->z) << std::endl;
        }
        output << "$EndNodes" << std::endl;
    }

    if (m_numberOfElements>0) {
        std::vector<int> elementIndices;
        getElementIndices(elementIndices);
        output << "$Elements" << std::endl;
        output << m_numberOfElements << std::endl;
        for (int i = 0; i< m_numberOfElements; i++) {
            Element& element = *m_elements[elementIndices[i]];
            int ntags;
            if (element.partitions.size())
                ntags = 3+element.partitions.size();
            else
                ntags = 2;
            output << elementIndices[i] << " " <<
                      element.type << " " <<
                      ntags << " " <<
                      element.physicalEntity << " " <<
                      element.elementaryEntity;
            for (int j = 0; j < element.partitions.size(); j++)
                output << " " << element.partitions[j];
            for (int j = 0; j < element.nodes.size(); j++)
                output << " " << element.nodes[j];
            output << std::endl;
        }
        output << "$EndElements" << std::endl;
    }


    if (!(m_periodicEntities.size() == 0) && !(m_periodicNodes.size() == 0) )
    {

        output << "$Periodic" << std::endl;
        output << m_periodicEntities.size() << std::endl;
        for (int i = 0; i < m_periodicEntities.size(); ++i){
            output << m_periodicEntities[i].dimension << " " <<
                      m_periodicEntities[i].slaveTag << " " <<
                      m_periodicEntities[i].masterTag << std::endl;
        }
        output << m_periodicNodes.size() << std::endl;
        for (int i = 0; i < m_periodicNodes.size(); ++i){
            output << m_periodicNodes[i].slaveNode << " " <<
                      m_periodicNodes[i].masterNode << std::endl;
        }
        output << "$EndPeriodic" << std::endl;

    }

    if (m_physicalNames.size() > 0) {
        output << "$PhysicalNames" << std::endl;
        for (int i = 0;i < m_physicalNames.size(); ++i) {
            output << m_physicalNames[i].dimension << " " <<
                      m_physicalNames[i].number << " " <<
                      m_physicalNames[i].name << std::endl;
        }
        output << "$EndPhysicalNames" << std::endl;
    }

    if (m_nodeDataSets.size() > 0) {
        for (int i = 0; i < m_nodeDataSets.size(); ++i) {

            NodeDataSet& nodeDataSet = *m_nodeDataSets[i];
            output << "$NodeData" << std::endl;
            output << nodeDataSet.stringTags.size() << std::endl;
            for (int i = 0; i< nodeDataSet.stringTags.size(); ++i) {
                output << '\"'+nodeDataSet.stringTags[i]+'\"' << std::endl;
            }
            output << nodeDataSet.realTags.size() << std::endl;
            for (int i = 0; i < nodeDataSet.realTags.size(); ++i) {
                output << boost::lexical_cast<std::string>(nodeDataSet.realTags[i]) << std::endl;
            }
            output << 4 << std::endl; // Number of integer tags
            output << nodeDataSet.timeStep << std::endl;
            output << nodeDataSet.numberOfFieldComponents << std::endl;
            output << nodeDataSet.values.size() << std::endl;
            output << nodeDataSet.partition << std::endl;
            for (int i = 0; i < nodeDataSet.values.size(); ++i ) {
                output << nodeDataSet.nodeIndices[i];
                for (int j = 0; j < nodeDataSet.values[i].size(); ++j) {
                    output << " " << boost::lexical_cast<std::string>(nodeDataSet.values[i][j]);
                }
                output << std::endl;
            }
            output << "$EndNodeData" << std::endl;
        }
    }

    if (m_elementDataSets.size() > 0) {
        for (int i = 0; i < m_elementDataSets.size(); ++i) {

            ElementDataSet& elementDataSet = *m_elementDataSets[i];
            output << "$ElementData" << std::endl;
            output << elementDataSet.stringTags.size() << std::endl;
            for (int i = 0; i< elementDataSet.stringTags.size(); ++i) {
                output << '\"'+elementDataSet.stringTags[i]+'\"' << std::endl;
            }
            output << elementDataSet.realTags.size() << std::endl;
            for (int i = 0; i < elementDataSet.realTags.size(); ++i) {
                output << boost::lexical_cast<std::string>(elementDataSet.realTags[i]) << std::endl;
            }
            output << 4 << std::endl; // Number of integer tags
            output << elementDataSet.timeStep << std::endl;
            output << elementDataSet.numberOfFieldComponents << std::endl;
            output << elementDataSet.values.size() << std::endl;
            output << elementDataSet.partition << std::endl;
            for (int i = 0; i < elementDataSet.values.size(); ++i ) {
                output << elementDataSet.elementIndices[i];
                for (int j = 0; j < elementDataSet.values[i].size(); ++j) {
                    output << " " << boost::lexical_cast<std::string>(elementDataSet.values[i][j]);
                }
                output << std::endl;
            }
            output << "$EndElementData" << std::endl;
        }
    }

    if (m_elementNodeDataSets.size() > 0) {
        for (int i = 0; i < m_elementNodeDataSets.size(); ++i) {

            ElementNodeDataSet& elementNodeDataSet = *m_elementNodeDataSets[i];
            output << "$ElementNodeData" << std::endl;
            output << elementNodeDataSet.stringTags.size() << std::endl;
            for (int i = 0; i< elementNodeDataSet.stringTags.size(); ++i) {
                output << '\"'+elementNodeDataSet.stringTags[i]+'\"' << std::endl;
            }
            output << elementNodeDataSet.realTags.size() << std::endl;
            for (int i = 0; i < elementNodeDataSet.realTags.size(); ++i) {
                output << boost::lexical_cast<std::string>(elementNodeDataSet.realTags[i]) << std::endl;
            }
            output << 4 << std::endl; // Number of integer tags
            output << elementNodeDataSet.timeStep << std::endl;
            output << elementNodeDataSet.numberOfFieldComponents << std::endl;
            output << elementNodeDataSet.values.size() << std::endl;
            output << elementNodeDataSet.partition << std::endl;
            for (int i = 0; i < elementNodeDataSet.values.size(); ++i ) {
                output << elementNodeDataSet.elementIndices[i] << " " <<
                          elementNodeDataSet.values[i].size();
                for (int j = 0; j < elementNodeDataSet.values[i].size(); ++j) {
                    for (int k = 0; k < elementNodeDataSet.values[i][j].size(); ++k)
                        output << " " << boost::lexical_cast<std::string>(elementNodeDataSet.values[i][j][k]);
                }
                output << std::endl;
            }
            output << "$EndElementNodeData" << std::endl;
        }
    }


    for (int i = 0; i < m_interpolationSchemeSets.size(); ++i) {

        InterpolationSchemeSet& scheme = *m_interpolationSchemeSets[i];
        output << "$InterpolationScheme" << std::endl;
        output << '\"'+scheme.name+'\"' << std::endl;
        output << 1 << std::endl; // Only one element topology supported right now.
        output << scheme.topology << std::endl;
        output << scheme.values.size() << std::endl;

        for (int n = 0; n < scheme.values.size(); ++n) {
            output << scheme.nrows[n] << std::endl;
            output << scheme.ncols[n] << std::endl;

            for (int j = 0; j < scheme.nrows[n]; ++j) {
                output << scheme.values[n][j*scheme.ncols[n]];
                for (int k = 1; k < scheme.ncols[n]; ++k)
                    output << " " << scheme.values[n][j*scheme.ncols[n]+k];
                output << std::endl;
            }
        }
        output << "$EndInterpolationScheme" << std::endl;
    }





}
void GmshData::write(const std::string& fileName) const {

    std::ofstream out;
    out.open(fileName.c_str(), std::ios::trunc);
    write(out);
    out.close();

}


GmshData GmshData::read(std::istream& input, int elementType, int physicalEntity) {

    bool haveMeshFormat = false;
    bool haveNodes = false;
    bool haveElements = false;
    bool havePeriodic = false;
    bool havePhysicalNames = false;

    GmshData gmshData;

    std::string line;
    while (std::getline(input,line)) {

        if (line == "$MeshFormat") {
            if (haveMeshFormat) throw std::runtime_error(
                        "GmshData::read(): MeshFormat Section appears more than once.");
            std::cout << "Reading MeshFormat..." << std::endl;
            std::getline(input,line);
            StringVector tokens = stringTokens(line);
            if (tokens.size() != 3) throw std::runtime_error(
                        "GmshData::read(): Wrong format of MeshFormat");
            if (tokens[0] != "2" & tokens[0] != "2.2") throw std::runtime_error(
                        "GmshData::read(): Version of MSH file not supported.");
            if (boost::lexical_cast<int>(tokens[1]) != 0) throw std::runtime_error(
                        "GmshData::read(): File Type not supported.");
            int dataSize = boost::lexical_cast<int>(tokens[2]);
            if (dataSize != sizeof(double)) throw std::runtime_error(
                        "MeshFormat::read(): Data size not supported.");
            std::getline(input,line);
            if (line != "$EndMeshFormat") throw std::runtime_error(
                        "GmshData::read(): Error reading MeshFormat section.");
            haveMeshFormat = true;


        }
        else if (line == "$Nodes") {
            if (haveNodes) throw std::runtime_error(
                        "GmshData::read(): Nodes section appears more than once. ");
            std::cout << "Reading Nodes..." << std::endl;
            std::getline(input,line);
            int numberOfNodes = boost::lexical_cast<int>(line);
            gmshData.reserveNumberOfNodes(numberOfNodes);
            for (int i = 0; i < numberOfNodes; ++i){
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 4) throw std::runtime_error(
                            "GmshData::read(): Wrong format of node definition detected.");
                int index = boost::lexical_cast<int>(tokens[0]);
                double x = boost::lexical_cast<double>(tokens[1]);
                double y = boost::lexical_cast<double>(tokens[2]);
                double z = boost::lexical_cast<double>(tokens[3]);
                gmshData.addNode(index,x,y,z);
            }
            std::getline(input,line);
            if (line != "$EndNodes") throw std::runtime_error(
                        "GmshData::read(): Error reading Nodes section. ");
            haveNodes = true;
        }
        else if (line == "$Elements") {
            if (haveElements) throw std::runtime_error(
                        "GmshData::read(): Elements section appears more than once.");
            std::cout << "Reading Elements..." << std::endl;
            std::getline(input,line);
            int numberOfElements = boost::lexical_cast<int>(line);
            for (int i = 0; i < numberOfElements; ++i) {
                std::getline(input,line);
                    StringVector tokens = stringTokens(line);
                    int index = boost::lexical_cast<int>(tokens.at(0));
                    int currentElementType = boost::lexical_cast<int>(tokens.at(1));
                    int ntags = boost::lexical_cast<int>(tokens.at(2));
                    int currentPhysicalEntity = boost::lexical_cast<int>(tokens.at(3));
                    int elementaryEntity = boost::lexical_cast<int>(tokens.at(4));
                    int npartitions = 0;
                    if (ntags > 5) npartitions = boost::lexical_cast<int>(tokens.at(5));
                    std::vector<int> partitions;
                    for (int i = 0;i < npartitions; ++i)
                        partitions.push_back(boost::lexical_cast<int>(tokens.at(6+i)));
                    std::vector<int> nodes;
                    for (int i = 3 + ntags; i < tokens.size(); ++i)
                        nodes.push_back(boost::lexical_cast<int>(tokens.at(i)));
                    if ((elementType == -1 || currentElementType == elementType) &&
                            (physicalEntity == -1 || currentPhysicalEntity == physicalEntity)) {
                            gmshData.addElement(index,currentElementType,nodes,physicalEntity,elementaryEntity,partitions);
                            gmshData.m_elementIndices.insert(index);
                    }
                }
            std::getline(input,line);
            if (line != "$EndElements") throw std::runtime_error(
                        "GmshData::read(): Error reading Elements section. ");
            haveElements = true;

            }
        else if (line == "$Periodic") {
            if (havePeriodic) throw std::runtime_error(
                        "GmshData::read(): Periodic section appears more than once.");
            std::cout << "Reading Periodic..." << std::endl;
            std::getline(input,line);
            int numberOfPeriodicEntities = boost::lexical_cast<int>(line);
            for (int i = 0;i < numberOfPeriodicEntities; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size()!= 3) throw std::runtime_error(
                            "GmshData::read(): Wrong format for periodic entities.");
                int dimension = boost::lexical_cast<int>(tokens[0]);
                int slaveTag = boost::lexical_cast<int>(tokens[1]);
                int masterTag = boost::lexical_cast<int>(tokens[2]);
                gmshData.addPeriodicEntity(dimension,slaveTag,masterTag);
            }

            std::getline(input,line);
            int numberOfPeriodicNodes = boost::lexical_cast<int>(line);
            for (int i = 0; i < numberOfPeriodicNodes; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 2) throw std::runtime_error(
                            "GmshData::read(): Wrong format for periodic nodes.");
                int slaveNode = boost::lexical_cast<int>(tokens[0]);
                int masterNode = boost::lexical_cast<int>(tokens[1]);
                gmshData.addPeriodicNode(slaveNode,masterNode);
            }
            std::getline(input,line);
            if (line != "$EndPeriodic") throw std::runtime_error(
                        "GmshData::read(): Error reading Periodic section.");
            havePeriodic = true;
        }
        else if (line == "$PhysicalNames") {
            if (havePhysicalNames) throw std::runtime_error(
                        "GmshData::read(): PhysicalNames section appears more than once.");
            std::cout << "Reading PhysicalNames..." << std::endl;
            std::getline(input,line);
            int numberOfPhysicalNames = boost::lexical_cast<int>(line);
            for (int i = 0; i< numberOfPhysicalNames; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 3) throw std::runtime_error(
                            "PhysicalNamesSet::read(): Wrong format for physical names.");
                int dimension = boost::lexical_cast<int>(tokens[0]);
                int number = boost::lexical_cast<int>(tokens[1]);
                std::string name = tokens[2];
                gmshData.addPhysicalName(dimension,number,name);
            }
            std::getline(input,line);
            if (line != "$EndPhysicalNames") throw std::runtime_error(
                        "GmshData::read(): Error reading PhysicalNames section.");
            havePhysicalNames = true;


        }
        else if (line == "$NodeData") {

            std::cout << "Reading NodeData..." << std::endl;
            std::getline(input,line);
            int numberOfStringTags = boost::lexical_cast<int>(line);
            std::vector<std::string> stringTags;
            for (int i = 0; i < numberOfStringTags; ++i) {
                std::getline(input,line);
                line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
                stringTags.push_back(line);
            }

            // Real tags
            std::getline(input,line);
            int numberOfRealTags = boost::lexical_cast<int>(line);
            std::vector<double> realTags;
            for (int i = 0; i < numberOfRealTags; ++i ) {
                std::getline(input,line);
                realTags.push_back(boost::lexical_cast<double>(line));
            }

            // Integer tags
            std::getline(input,line);
            int numberOfIntegerTags = boost::lexical_cast<int>(line);
            if (numberOfIntegerTags < 3) throw std::runtime_error(
                        "GmshData::read(): At least 3 integer tags required.");
            std::vector<int> integerTags;
            for (int i = 0; i < numberOfIntegerTags; ++i) {
                std::getline(input,line);
                integerTags.push_back(boost::lexical_cast<int>(line));
            }
            int timeStep = integerTags[0];
            int numberOfFieldComponents = integerTags[1];
            int numberOfNodes = integerTags[2];
            int partition = 0;
            if (integerTags.size() > 3) partition = integerTags[3];
            int dataSetIndex = gmshData.numberOfNodeDataSets();
            gmshData.addNodeDataSet(stringTags,realTags,numberOfFieldComponents,numberOfNodes,timeStep,partition);

            for (int i = 0; i < numberOfNodes; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 1+numberOfFieldComponents) throw std::runtime_error(
                            "GmshData::read(): Data has wrong format.");
                int index = boost::lexical_cast<int>(tokens[0]);
                std::vector<double> values;
                for (int j = 1; j < tokens.size(); ++j) {
                    values.push_back(boost::lexical_cast<double>(tokens[j]));
                }
                gmshData.addNodeData(dataSetIndex,index,values);
            }

            std::getline(input,line);
            if (line != "$EndNodeData") throw std::runtime_error(
                        "GmshData::read(): Error reading NodeData section.");


        }
        else if (line == "$ElementData") {

            std::cout << "Reading ElementData..." << std::endl;
            std::getline(input,line);
            int numberOfStringTags = boost::lexical_cast<int>(line);
            std::vector<std::string> stringTags;
            for (int i = 0; i < numberOfStringTags; ++i) {
                std::getline(input,line);
                line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
                stringTags.push_back(line);
            }

            // Real tags
            std::getline(input,line);
            int numberOfRealTags = boost::lexical_cast<int>(line);
            std::vector<double> realTags;
            for (int i = 0; i < numberOfRealTags; ++i ) {
                std::getline(input,line);
                realTags.push_back(boost::lexical_cast<double>(line));
            }

            // Integer tags
            std::getline(input,line);
            int numberOfIntegerTags = boost::lexical_cast<int>(line);
            if (numberOfIntegerTags < 3) throw std::runtime_error(
                        "GmshData::read(): At least 3 integer tags required.");
            std::vector<int> integerTags;
            for (int i = 0; i < numberOfIntegerTags; ++i) {
                std::getline(input,line);
                integerTags.push_back(boost::lexical_cast<int>(line));
            }
            int timeStep = integerTags[0];
            int numberOfFieldComponents = integerTags[1];
            int numberOfElements = integerTags[2];
            int partition = 0;
            if (integerTags.size() > 3) partition = integerTags[3];
            int dataSetIndex = gmshData.numberOfElementDataSets();
            gmshData.addElementDataSet(stringTags,realTags,numberOfFieldComponents,numberOfElements,timeStep,partition);

            for (int i = 0; i < numberOfElements; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 1+numberOfFieldComponents) throw std::runtime_error(
                            "GmshData::read(): Data has wrong format.");
                int index = boost::lexical_cast<int>(tokens[0]);
                std::vector<double> values;
                for (int j = 1; j < tokens.size(); ++j) {
                    values.push_back(boost::lexical_cast<double>(tokens[j]));
                }
                if (gmshData.m_elementIndices.find(index)!=gmshData.m_elementIndices.end())
                    gmshData.addElementData(dataSetIndex,index,values);
            }

            std::getline(input,line);
            if (line != "$EndElementData") throw std::runtime_error(
                        "GmshData::read(): Error reading ElementData section.");


        }
        else if (line == "$ElementNodeData") {

            std::cout << "Reading ElementNodeData..." << std::endl;
            std::getline(input,line);
            int numberOfStringTags = boost::lexical_cast<int>(line);
            std::vector<std::string> stringTags;
            for (int i = 0; i < numberOfStringTags; ++i) {
                std::getline(input,line);
                line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
                stringTags.push_back(line);
            }

            // Real tags
            std::getline(input,line);
            int numberOfRealTags = boost::lexical_cast<int>(line);
            std::vector<double> realTags;
            for (int i = 0; i < numberOfRealTags; ++i ) {
                std::getline(input,line);
                realTags.push_back(boost::lexical_cast<double>(line));
            }

            // Integer tags
            std::getline(input,line);
            int numberOfIntegerTags = boost::lexical_cast<int>(line);
            if (numberOfIntegerTags < 3) throw std::runtime_error(
                        "GmshData::read(): At least 3 integer tags required.");
            std::vector<int> integerTags;
            for (int i = 0; i < numberOfIntegerTags; ++i) {
                std::getline(input,line);
                integerTags.push_back(boost::lexical_cast<int>(line));
            }
            int timeStep = integerTags[0];
            int numberOfFieldComponents = integerTags[1];
            int numberOfElements = integerTags[2];
            int partition = 0;
            if (integerTags.size() > 3) partition = integerTags[3];
            int dataSetIndex = gmshData.numberOfElementDataSets();
            gmshData.addElementNodeDataSet(stringTags,realTags,numberOfFieldComponents,numberOfElements,timeStep,partition);

            for (int i = 0; i < numberOfElements; ++i) {
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() < 2) throw std::runtime_error(
                            "GmshData::read(): Data has wrong format.");
                int numberOfNodes = boost::lexical_cast<int>(tokens[1]);
                if (tokens.size() != 2+numberOfFieldComponents*numberOfNodes) throw std::runtime_error(
                            "GmshData::read(): Data has wrong format.");
                int index = boost::lexical_cast<int>(tokens[0]);
                std::vector<std::vector<double> > values;
                for (int j = 0; j < numberOfNodes; ++j) {
                    values.push_back(std::vector<double>());
                    for (int k = 0; k < numberOfFieldComponents; ++k)
                        values[j].push_back(boost::lexical_cast<double>(tokens[2+j*numberOfFieldComponents+k]));
                }
                if (gmshData.m_elementIndices.find(index)!=gmshData.m_elementIndices.end())
                    gmshData.addElementNodeData(dataSetIndex,index,values);
            }

            std::getline(input,line);
            if (line != "$EndElementNodeData") throw std::runtime_error(
                        "GmshData::read(): Error reading ElementNodeData section.");


        }
        else if (line == "$InterpolationSchemeSet") {

            std::cout << "Reading InterpolationSchemSet..." << std::endl;
            std::getline(input,line);
            line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
            std::string name = line;
            std::getline(input,line);
            if (boost::lexical_cast<int>(line) != 1) throw std::runtime_error(
                        "GmshData::read(): Only one topology is currently supported.");
            std::getline(input,line);
            int topology = boost::lexical_cast<int>(line);
            int dataSetIndex = gmshData.numberOfInterpolationSchemeSets();
            gmshData.addInterpolationSchemeSet(name,topology);
            std::getline(input,line);
            int numberOfInterpolationMatrices = boost::lexical_cast<int>(line);
            for (int i = 0; i < numberOfInterpolationMatrices; ++i) {
                std::vector<double> matrix;
                std::getline(input,line);
                StringVector tokens = stringTokens(line);
                if (tokens.size() != 2) throw std::runtime_error(
                            "GmshData::read(): Wrong format of interpolation matrices.");
                int nrows = boost::lexical_cast<int>(tokens[0]);
                int ncols = boost::lexical_cast<int>(tokens[1]);
                matrix.reserve(nrows*ncols);
                for (int j = 0; j < nrows; ++j) {
                    std::getline(input,line);
                    StringVector tokens = stringTokens(line);
                    if (tokens.size() != ncols) throw std::runtime_error(
                                "GmshData::read(): Wrong format of interpolation matrices.");
                    for (int k = 0;k < ncols; ++k)
                        matrix.push_back(boost::lexical_cast<double>(tokens[k]));
                }
                gmshData.addInterpolationMatrix(dataSetIndex,nrows,ncols,matrix);
            }
            std::getline(input,line);
            if (line != "$EndInterpolationSchemeSet") throw std::runtime_error(
                        "GmshData::read(): Error reading InterpolationSchemeSet section.");
        }

    }
    return gmshData;

}

GmshData GmshData::read(const std::string& fileName, int elementType, int physicalEntity) {

    std::ifstream(input);
    input.open(fileName.c_str());
    GmshData gmshData = read(input,elementType,physicalEntity);
    input.close();
    return gmshData;

}

void GmshData::resetNodeDataSets() {
    m_nodeDataSets.clear();
}

void GmshData::resetElementDataSets() {
    m_elementDataSets.clear();
}

void GmshData::resetElementNodeDataSets(){
    m_elementNodeDataSets.clear();
}

void GmshData::resetDataSets() {

    resetNodeDataSets();
    resetElementDataSets();
    resetElementNodeDataSets();

}


GmshIo::GmshIo(const shared_ptr<const Grid>& grid) : m_grid(grid) {

    if (!m_grid) throw std::runtime_error(
                "GmshIo(): Grid is not initialized.");

    std::auto_ptr<GridView> viewPtr = m_grid->leafView();
    const GridView& view = *viewPtr;

    m_nodePermutation.resize(view.entityCount(2));
    m_inverseNodePermutation.resize(view.entityCount(2)+1,-1);
    m_elementPermutation.resize(view.entityCount(0));
    m_inverseElementPermutation.resize(view.entityCount(0)+1,-1);

    for (int i = 0; i < m_nodePermutation.size(); ++i)
        m_nodePermutation[i] = i+1;
    for (int i = 0; i < m_elementPermutation.size(); ++i)
        m_elementPermutation[i] = i+1;
    for (int i = 1; i < m_inverseElementPermutation.size(); ++i)
        m_inverseElementPermutation[i] = i-1;
    for (int i = 1; i < m_inverseNodePermutation.size(); ++i)
        m_inverseNodePermutation[i] = i-1;

    // Now write out the grid to the gmshData object.

    const IndexSet& indexSet = view.indexSet();

    std::auto_ptr<EntityIterator<2> > nodeIterator = view.entityIterator<2>();
    while (!nodeIterator->finished()) {
        const Entity<2>& node = nodeIterator->entity();
        int index = m_nodePermutation[indexSet.entityIndex(node)];
        const Geometry& geom = node.geometry();
        arma::Col<double> coords(3,1);
        geom.getCenter(coords);
        m_gmshData.addNode(index,coords(0),coords(1),coords(2));
        nodeIterator->next();
    }

    std::auto_ptr<EntityIterator<0> > elementIterator = view.entityIterator<0>();

    while (!elementIterator->finished()) {

        const Entity<0>& element = elementIterator->entity();
        int index = m_elementPermutation[indexSet.entityIndex(element)];
        std::vector<int> nodeIndices;
        int numberOfNodes = element.subEntityCount<2>();
        nodeIndices.reserve(numberOfNodes);
        for (int i = 0; i < numberOfNodes; ++i)
            nodeIndices.push_back(m_nodePermutation[indexSet.subEntityIndex(element,i,2)]);
        int physicalIndex = element.domain();
        m_gmshData.addElement(index,2,nodeIndices,physicalIndex,0);
        elementIterator->next();
    }


}

GmshIo::GmshIo(GmshData gmshData) : m_gmshData(gmshData){}

GmshIo::GmshIo(std::string fileName, int physicalEntity) : m_gmshData(GmshData::read(fileName,2,physicalEntity)){}


shared_ptr<const Grid> GmshIo::grid() const {

    if (m_grid) return m_grid;

    // Sanity check on the Gmsh Data

    std::vector<int> nodeIndices;
    m_gmshData.getNodeIndices(nodeIndices);
    std::vector<int> elementIndices;
    m_gmshData.getElementIndices(elementIndices);

    int numberOfGmshNodes = m_gmshData.numberOfNodes();
    int numberOfGmshElements = m_gmshData.numberOfElements();

    int maxNodeIndex = *(std::max_element(nodeIndices.begin(),nodeIndices.end()));
    int maxElementIndex = *(std::max_element(elementIndices.begin(),elementIndices.end()));

    std::vector<std::vector<int> > elements;

    m_inverseNodePermutation.resize(maxNodeIndex+1,-1);
    m_inverseElementPermutation.resize(maxElementIndex+1,-1);
    std::vector<int> domainIndices;

    for (int i = 0; i < elementIndices.size(); ++i) {
        int physicalEntity;
        int elementaryEntity;
        int elementType;
        int index = elementIndices[i];
        std::vector<int> nodes;

        m_gmshData.getElement(index,elementType,nodes,physicalEntity,elementaryEntity);
        if (elementType == 2) {
            int elementIndex = m_elementPermutation.size();
            m_elementPermutation.push_back(index);
            m_inverseElementPermutation[index] = elementIndex;
            elements.push_back(std::vector<int>());
            domainIndices.push_back(physicalEntity);
            std::vector<int>& currentNodes = elements.back();
            for (int j = 0; j < nodes.size(); ++j) {
                int nodeIndex;
                if (m_inverseNodePermutation[nodes[j]] == -1) { // Node not yet assigned
                    nodeIndex = m_nodePermutation.size();
                    m_inverseNodePermutation[nodes[j]] = nodeIndex;
                    m_nodePermutation.push_back(nodes[j]);
                }
                else {
                    nodeIndex = m_inverseNodePermutation[nodes[j]];
                }
                currentNodes.push_back(nodeIndex);
            }
        }
    }
    arma::Mat<double> armaNodes(3,m_nodePermutation.size());
    arma::Mat<int> armaElements(3,m_elementPermutation.size());
    for (int i = 0; i < m_nodePermutation.size();++i){
        double x,y,z;
        m_gmshData.getNode(m_nodePermutation[i],x,y,z);
        armaNodes(0,i) = x;
        armaNodes(1,i) = y;
        armaNodes(2,i) = z;
    }
    for (int i = 0; i < m_elementPermutation.size(); ++i) {
        armaElements(0,i) = elements[i][0];
        armaElements(1,i) = elements[i][1];
        armaElements(2,i) = elements[i][2];
    }
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    m_grid =  GridFactory::createGridFromConnectivityArrays
                (params,armaNodes,armaElements,domainIndices);
    return m_grid;

}

const std::vector<int>& GmshIo::nodePermutation() const {

    return m_nodePermutation;

}

const std::vector<int>& GmshIo::elementPermutation() const {

    return m_elementPermutation;

}

const std::vector<int>& GmshIo::inverseNodePermutation() const {

    return m_inverseNodePermutation;

}

const std::vector<int>& GmshIo::inverseElementPermutation() const {

    return m_inverseElementPermutation;

}

const GmshData& GmshIo::gmshData() const {

    return m_gmshData;

}

GmshData& GmshIo::gmshData() {

    return m_gmshData;

}

void GmshIo::write(std::string fileName) const {

    m_gmshData.write(fileName);

}

void GmshIo::resetNodeDataSets() {
    m_gmshData.resetNodeDataSets();
}

void GmshIo::resetElementDataSets() {
    m_gmshData.resetElementDataSets();
}

void GmshIo::resetElementNodeDataSets() {
    m_gmshData.resetElementNodeDataSets();
}

void GmshIo::resetDataSets() {
    m_gmshData.resetDataSets();
}


template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType,ResultType> gridFunctionFromGmsh(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            const GmshIo& gmshIo, shared_ptr<const Grid> grid,
            GmshPostData::Type gmshPostDataType,
            int index = 0) {

    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    typedef typename ScalarTraits<ResultType>::RealType ResultRealType;

    int maxIndex;
    int numberOfFieldComponents;
    std::vector<std::string> stringTags;
    std::vector<double> realTags;

    shared_ptr<const Grid> currentGrid;
    if (grid)
        currentGrid = grid;
    else
        currentGrid = gmshIo.grid();

    if (gmshPostDataType == GmshPostData::ELEMENT) {

        maxIndex = gmshIo.gmshData().numberOfElementDataSets();
        if (index >= maxIndex) throw std::runtime_error(
                    "gridFunctionFromGmsh(): There is no dataset of the given type with this index.");
        std::vector<int> elementIndices;
        std::vector<std::vector<double> > values;
        gmshIo.gmshData().getElementDataSet(index,stringTags,realTags,numberOfFieldComponents,elementIndices,values);
        for (int i = 0; i < values.size(); ++i)
            if (values[i].size() != 1) throw std::runtime_error(
                        "gridFunctionFromGmsh(): Currently, only single component data are supported.");
        const std::vector<int>& inverseElementPermutation = gmshIo.inverseElementPermutation();
        int numberOfElements = currentGrid->leafView()->entityCount(0);
        if (elementIndices.size() != numberOfElements) throw
            std::runtime_error("gridFunctionFromGmsh(): Number of element indices does not agree with the number of grid elements.");
        arma::Col<ResultType> coefficients(numberOfElements);
        for (int i = 0; i < elementIndices.size(); ++i)
            coefficients(acc(inverseElementPermutation,elementIndices[i])) =
                    (ResultRealType) values[i][0];

        const shared_ptr<const Space<BasisFunctionType> > space
                (new PiecewiseConstantScalarSpace<BasisFunctionType>(currentGrid));
        return GridFunction<BasisFunctionType,ResultType>(context,space,coefficients);

    }
    else if (gmshPostDataType == GmshPostData::ELEMENT_NODE) {

        maxIndex = gmshIo.gmshData().numberOfElementNodeDataSets();
        if (index >= maxIndex) throw std::runtime_error(
                    "gridFunctionFromGmsh(): There is no dataset of the given type with this index.");
        std::vector<int> elementIndices;
        std::vector<std::vector<std::vector<double> > > values;
        gmshIo.gmshData().getElementNodeDataSet(index,stringTags,realTags,numberOfFieldComponents,elementIndices,values);
        for (int i = 0; i < values.size(); ++i) {
            if (values[i].size() != 3) throw std::runtime_error(
                        "gridFunctionFromGmsh(): Only 3-noded triangles are supported.");
            for (int j = 0; j < values[i].size(); ++j)
                if (values[i][j].size() != 1) throw std::runtime_error(
                            "gridFunctionFromGmsh(): Currently, only single component data are supported.");
        }
        const std::vector<int>& inverseElementPermutation = gmshIo.inverseElementPermutation();
        int numberOfElements = currentGrid->leafView()->entityCount(0);
        std::cout << elementIndices.size() << " " << numberOfElements << std::endl;
        if (elementIndices.size() != numberOfElements) throw
            std::runtime_error("gridFunctionFromGmsh(): Number of element indices does not agree with the number of grid elements.");
        arma::Col<ResultType> coefficients(3*numberOfElements);
        for (int i = 0; i < elementIndices.size(); ++i)
            for (int j = 0; j < 3; ++j) {
                coefficients(3*acc(inverseElementPermutation,elementIndices[i])+j) =
                        (ResultRealType) values[i][j][0];
            }
        const shared_ptr<const Space<BasisFunctionType> > space
                (new PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>(currentGrid));
        return GridFunction<BasisFunctionType,ResultType>(context,space,coefficients);


    }
    else if (gmshPostDataType == GmshPostData::NODE) {

        maxIndex = gmshIo.gmshData().numberOfNodes();
        if (index >= maxIndex) throw std::runtime_error(
                    "gridFunctionFromGmsh(): There is no dataset of the given type with this index.");
        std::vector<int> nodeIndices;
        std::vector<std::vector<double> > values;
        gmshIo.gmshData().getNodeDataSet(index,stringTags,realTags,numberOfFieldComponents,nodeIndices,values);
        for (int i = 0; i < values.size(); ++i)
            if (values[i].size() != 1) throw std::runtime_error(
                        "gridFunctionFromGmsh(): Currently, only single component data are supported.");
        const std::vector<int>& inverseNodePermutation = gmshIo.inverseNodePermutation();
        int numberOfNodes = currentGrid->leafView()->entityCount(2);
        if (nodeIndices.size() != numberOfNodes) throw
            std::runtime_error("gridFunctionFromGmsh(): Number of node indices does not agree with the number of grid elements.");
        arma::Col<ResultType> coefficients(numberOfNodes);
        for (int i = 0; i < nodeIndices.size(); ++i)
            coefficients(acc(inverseNodePermutation,nodeIndices[i])) =
                    (ResultRealType) values[i][0];

        const shared_ptr<const Space<BasisFunctionType> > space
                (new PiecewiseLinearContinuousScalarSpace<BasisFunctionType>(currentGrid));
        return GridFunction<BasisFunctionType,ResultType>(context,space,coefficients);


    }




}

template <typename BasisFunctionType, typename ResultType>
void exportToGmsh(GridFunction<BasisFunctionType,ResultType> gridFunction,
                  const char* dataLabel, GmshIo& gmshIo,
                  GmshPostData::Type gmshPostDataType,
                  std::string complexMode){


    if (complexMode != "real" && complexMode != "imag" &&
            complexMode != "abs" && complexMode != "all")
        throw std::runtime_error("exportToGmsh(): complexMode must be one of 'real', 'imag', or 'abs'");

    shared_ptr<const Space<BasisFunctionType> > space = gridFunction.space();
    if (!space)
        throw std::runtime_error("exportToGmsh() must not be "
                                 "called on an uninitialized GridFunction object");
    const size_t dimWorld = 3;
    if (space->grid()->dimWorld() != 3 || space->grid()->dim() != 2)
        throw std::runtime_error("exportToGmsh(): currently only "
                                 "2D grids in 3D spaces are supported");

    if (!gmshIo.grid())
        gmshIo = GmshIo(gridFunction.grid());
    if (gmshIo.grid()!= gridFunction.grid()) throw
        std::runtime_error("exportToGmsh(): gmshIo.grid()!= gridFunction.grid().");


    if (complexMode =="all") {

        std::string realLabel = std::string(dataLabel)+"_real";
        std::string imagLabel = std::string(dataLabel)+"_imag";
        std::string absLabel = std::string(dataLabel)+"_abs";

        exportToGmsh(gridFunction,realLabel.c_str(),gmshIo,gmshPostDataType,"real");
        exportToGmsh(gridFunction,imagLabel.c_str(),gmshIo,gmshPostDataType,"imag");
        exportToGmsh(gridFunction,absLabel.c_str(),gmshIo,gmshPostDataType,"abs");

        return;

    }

    GmshData& gmshData = gmshIo.gmshData();

    std::vector<std::string> stringTags;
    stringTags.push_back(dataLabel);
    std::vector<double> realTags;

    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    boost::array<arma::Mat<CoordinateType>, 5> localCoordsOnTriangles;

    localCoordsOnTriangles[0].set_size(2, 3);
    localCoordsOnTriangles[0].fill(0.);
    localCoordsOnTriangles[0](0, 1) = 1.;
    localCoordsOnTriangles[0](1, 2) = 1.;

    localCoordsOnTriangles[1] = localCoordsOnTriangles[0];

    localCoordsOnTriangles[2].set_size(2, 6);
    localCoordsOnTriangles[2].fill(0.);
    localCoordsOnTriangles[2](0, 1) = 1.;
    localCoordsOnTriangles[2](1, 2) = 1.;
    localCoordsOnTriangles[2](0, 3) = 0.5;
    localCoordsOnTriangles[2](0, 4) = 0.5;
    localCoordsOnTriangles[2](1, 4) = 0.5;
    localCoordsOnTriangles[2](1, 5) = 0.5;

    localCoordsOnTriangles[3].set_size(2, 10);
    localCoordsOnTriangles[3].fill(0.);
    localCoordsOnTriangles[3](0, 1) = 1.;
    localCoordsOnTriangles[3](1, 2) = 1.;
    localCoordsOnTriangles[3](0, 3) = 1./3.;
    localCoordsOnTriangles[3](0, 4) = 2./3.;
    localCoordsOnTriangles[3](0, 5) = 2./3.;
    localCoordsOnTriangles[3](1, 5) = 1./3.;
    localCoordsOnTriangles[3](0, 6) = 1./3.;
    localCoordsOnTriangles[3](1, 6) = 2./3.;
    localCoordsOnTriangles[3](1, 7) = 2./3.;
    localCoordsOnTriangles[3](1, 8) = 1./3.;
    localCoordsOnTriangles[3](0, 9) = 1./3.;
    localCoordsOnTriangles[3](1, 9) = 1./3.;

    localCoordsOnTriangles[4].set_size(2, 15);
    localCoordsOnTriangles[4].fill(0.);
    localCoordsOnTriangles[4](0, 1) = 1.;
    localCoordsOnTriangles[4](1, 2) = 1.;
    localCoordsOnTriangles[4](0, 3) = 0.25;
    localCoordsOnTriangles[4](0, 4) = 0.5;
    localCoordsOnTriangles[4](0, 5) = 0.75;
    localCoordsOnTriangles[4](0, 6) = 0.75;
    localCoordsOnTriangles[4](1, 6) = 0.25;
    localCoordsOnTriangles[4](0, 7) = 0.5;
    localCoordsOnTriangles[4](1, 7) = 0.5;
    localCoordsOnTriangles[4](0, 8) = 0.25;
    localCoordsOnTriangles[4](1, 8) = 0.75;
    localCoordsOnTriangles[4](1, 9) = 0.75;
    localCoordsOnTriangles[4](1, 10) = 0.5;
    localCoordsOnTriangles[4](1, 11) = 0.25;
    localCoordsOnTriangles[4](0, 12) = 0.25;
    localCoordsOnTriangles[4](1, 12) = 0.25;
    localCoordsOnTriangles[4](0, 13) = 0.5;
    localCoordsOnTriangles[4](1, 13) = 0.25;
    localCoordsOnTriangles[4](0, 14) = 0.25;
    localCoordsOnTriangles[4](1, 14) = 0.5;

    boost::array<int, 5> elementTypes;
    elementTypes[0] = 2;
    elementTypes[1] = 2;
    elementTypes[2] = 9;
    elementTypes[3] = 21;
    elementTypes[4] = 25;

    arma::Mat<ResultType> values;
    arma::Mat<CoordinateType> globalCoords;
    const GridView& view = space->gridView();
    int numberOfNodes = view.entityCount(2);
    int numberOfElements = view.entityCount(0);
    const IndexSet& indexSet = view.indexSet();

    if (gmshPostDataType == GmshPostData::NODE) {

        const std::vector<int>& nodePermutation = gmshIo.nodePermutation();
        arma::Mat<ResultType> values;
        gridFunction.evaluateAtSpecialPoints(VtkWriter::VERTEX_DATA,values);
        arma::Mat<CoordinateType> modifiedValues;
        if (complexMode == "real")
            modifiedValues = arma::real(values);
        else if (complexMode == "imag")
            modifiedValues = arma::imag(values);
        else if (complexMode == "abs")
            modifiedValues = arma::abs(values);
        int dataSetIndex = gmshData.numberOfNodeDataSets();
        gmshData.addNodeDataSet(stringTags,realTags,values.n_rows,numberOfNodes);
        for (int i = 0; i < numberOfNodes; ++i ) {
            std::vector<double> vals(values.n_rows);
            for (int j = 0; j < values.n_rows; ++j)
                vals[j] = modifiedValues(j,i);
            gmshData.addNodeData(dataSetIndex,nodePermutation[i],vals);
        }

    }
    else if (gmshPostDataType == GmshPostData::ELEMENT) {

        const std::vector<int>& elementPermutation = gmshIo.elementPermutation();
        arma::Mat<ResultType> values;
        gridFunction.evaluateAtSpecialPoints(VtkWriter::CELL_DATA,values);
        arma::Mat<CoordinateType> modifiedValues;
        if (complexMode == "real")
            modifiedValues = arma::real(values);
        else if (complexMode == "imag")
            modifiedValues = arma::imag(values);
        else if (complexMode == "abs")
            modifiedValues = arma::abs(values);
        int dataSetIndex = gmshData.numberOfElementDataSets();
        gmshData.addElementDataSet(stringTags,realTags,values.n_rows,numberOfElements);
        for (int i = 0; i < numberOfElements; ++i ) {
            std::vector<double> vals(values.n_rows);
            for (int j = 0; j < values.n_rows; ++j)
                vals[j] = modifiedValues(j,i);
            gmshData.addElementData(dataSetIndex,elementPermutation[i],vals);
        }

    }
    else if (gmshPostDataType == GmshPostData::ELEMENT_NODE) {

        const std::vector<int>& elementPermutation = gmshIo.elementPermutation();
        std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();

        int dataSetIndex = gmshData.numberOfElementNodeDataSets();
        gmshData.addElementNodeDataSet(stringTags,realTags,gridFunction.componentCount(),numberOfElements);
        while (!it->finished()) {

            const Entity<0>& element = it->entity();
            int elementIndex = indexSet.entityIndex(element);
            const Geometry& geo = element.geometry();
            arma::Mat<ResultType> values;
            arma::Mat<CoordinateType> modifiedValues;
            gridFunction.evaluate(element, localCoordsOnTriangles[0], values);
            if (complexMode == "real")
                modifiedValues = arma::real(values);
            else if (complexMode == "imag")
                modifiedValues = arma::imag(values);
            else if (complexMode == "abs")
                modifiedValues = arma::abs(values);

            std::vector<std::vector<double> > vals(values.n_cols);
            for (int j = 0; j < values.n_cols; ++j) {
                vals[j].reserve(values.n_rows);
                for (int i = 0; i < values.n_rows; ++i)
                    vals[j].push_back(modifiedValues(i,j));
            }
            gmshData.addElementNodeData(dataSetIndex,elementPermutation[elementIndex],vals);
            it->next();
        }
    }


}


template <typename BasisFunctionType, typename ResultType>
void exportToGmsh(GridFunction<BasisFunctionType,ResultType> gridFunction,
                  const char* dataLabel,
                  const char* fileName,
                  GmshPostData::Type gmshPostDataType,
                  std::string complexMode) {

    GmshIo gmshIo(gridFunction.grid());
    exportToGmsh(gridFunction,dataLabel,gmshIo,gmshPostDataType,complexMode);
    gmshIo.write(fileName);

}


//template <typename BasisFunctionType, typename ResultType>
//void exportToGmsh(GridFunction<BasisFunctionType,ResultType> gridFunction,
//                  const char* dataLabel, const char* fileName)
//{
//    shared_ptr<const Space<BasisFunctionType> > space = gridFunction.space();
//    if (!space)
//        throw std::runtime_error("exportToGmsh() must not be "
//                                 "called on an uninitialized GridFunction object");
//    const size_t dimWorld = 3;
//    if (space->grid()->dimWorld() != 3 || space->grid()->dim() != 2)
//        throw std::runtime_error("exportToGmsh(): currently only "
//                                 "2D grids in 3D spaces are supported");

//    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
//    boost::array<arma::Mat<CoordinateType>, 5> localCoordsOnTriangles;

//    localCoordsOnTriangles[0].set_size(2, 3);
//    localCoordsOnTriangles[0].fill(0.);
//    localCoordsOnTriangles[0](0, 1) = 1.;
//    localCoordsOnTriangles[0](1, 2) = 1.;

//    localCoordsOnTriangles[1] = localCoordsOnTriangles[0];

//    localCoordsOnTriangles[2].set_size(2, 6);
//    localCoordsOnTriangles[2].fill(0.);
//    localCoordsOnTriangles[2](0, 1) = 1.;
//    localCoordsOnTriangles[2](1, 2) = 1.;
//    localCoordsOnTriangles[2](0, 3) = 0.5;
//    localCoordsOnTriangles[2](0, 4) = 0.5;
//    localCoordsOnTriangles[2](1, 4) = 0.5;
//    localCoordsOnTriangles[2](1, 5) = 0.5;

//    localCoordsOnTriangles[3].set_size(2, 10);
//    localCoordsOnTriangles[3].fill(0.);
//    localCoordsOnTriangles[3](0, 1) = 1.;
//    localCoordsOnTriangles[3](1, 2) = 1.;
//    localCoordsOnTriangles[3](0, 3) = 1./3.;
//    localCoordsOnTriangles[3](0, 4) = 2./3.;
//    localCoordsOnTriangles[3](0, 5) = 2./3.;
//    localCoordsOnTriangles[3](1, 5) = 1./3.;
//    localCoordsOnTriangles[3](0, 6) = 1./3.;
//    localCoordsOnTriangles[3](1, 6) = 2./3.;
//    localCoordsOnTriangles[3](1, 7) = 2./3.;
//    localCoordsOnTriangles[3](1, 8) = 1./3.;
//    localCoordsOnTriangles[3](0, 9) = 1./3.;
//    localCoordsOnTriangles[3](1, 9) = 1./3.;

//    localCoordsOnTriangles[4].set_size(2, 15);
//    localCoordsOnTriangles[4].fill(0.);
//    localCoordsOnTriangles[4](0, 1) = 1.;
//    localCoordsOnTriangles[4](1, 2) = 1.;
//    localCoordsOnTriangles[4](0, 3) = 0.25;
//    localCoordsOnTriangles[4](0, 4) = 0.5;
//    localCoordsOnTriangles[4](0, 5) = 0.75;
//    localCoordsOnTriangles[4](0, 6) = 0.75;
//    localCoordsOnTriangles[4](1, 6) = 0.25;
//    localCoordsOnTriangles[4](0, 7) = 0.5;
//    localCoordsOnTriangles[4](1, 7) = 0.5;
//    localCoordsOnTriangles[4](0, 8) = 0.25;
//    localCoordsOnTriangles[4](1, 8) = 0.75;
//    localCoordsOnTriangles[4](1, 9) = 0.75;
//    localCoordsOnTriangles[4](1, 10) = 0.5;
//    localCoordsOnTriangles[4](1, 11) = 0.25;
//    localCoordsOnTriangles[4](0, 12) = 0.25;
//    localCoordsOnTriangles[4](1, 12) = 0.25;
//    localCoordsOnTriangles[4](0, 13) = 0.5;
//    localCoordsOnTriangles[4](1, 13) = 0.25;
//    localCoordsOnTriangles[4](0, 14) = 0.25;
//    localCoordsOnTriangles[4](1, 14) = 0.5;

//    boost::array<int, 5> elementTypes;
//    elementTypes[0] = 2;
//    elementTypes[1] = 2;
//    elementTypes[2] = 9;
//    elementTypes[3] = 21;
//    elementTypes[4] = 25;

//    arma::Mat<ResultType> values;
//    arma::Mat<CoordinateType> globalCoords;
//    const GridView& view = space->gridView();

//    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();

//    size_t nodeCount = 0;
//    size_t elementCount = 0;
//    std::stringstream nodes, elements, data;
//    while (!it->finished()) {
//        const Entity<0>& element = it->entity();
//        const Geometry& geo = element.geometry();
//        if (!element.type().isTriangle())
//            throw std::runtime_error(
//                "GridFunction::exportToGmsh(): "
//                "at present only triangular elements are supported");
//        int order = space->shapeset(element).order();
//        if (order >= localCoordsOnTriangles.size())
//            throw std::runtime_error(
//                "GridFunction::exportToGmsh(): "
//                "Gmsh does not support triangular elements of order larger than 5");
//        gridFunction.evaluate(element, localCoordsOnTriangles[order], values);
//        geo.local2global(localCoordsOnTriangles[order], globalCoords);

//        const size_t pointCount = localCoordsOnTriangles[order].n_cols;
//        for (size_t p = 0; p < pointCount; ++p) {
//            nodes << nodeCount + 1 + p;
//            for (size_t d = 0; d < dimWorld; ++d)
//                nodes  << ' ' << globalCoords(d, p);
//            nodes << '\n';
//        }

//        elements << ++elementCount << ' ' << elementTypes[order] << " 2 1 1";
//        for (size_t p = 0; p < pointCount; ++p)
//            elements << ' ' << nodeCount + 1 + p;
//        elements << '\n';

//        for (size_t p = 0; p < pointCount; ++p) {
//            data << nodeCount + 1 + p;
//            for (size_t d = 0; d < values.n_rows; ++d)
//                data << ' ' << realPart(values(d, p)); // TODO: export imag. part too
//            data << '\n';
//        }
//        nodeCount += pointCount;
//        it->next();
//    }

//    std::ofstream fout(fileName);
//    fout << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
//    fout << nodeCount << '\n' << nodes.str();
//    fout << "$EndNodes\n$Elements\n";
//    fout << elementCount << '\n' << elements.str();
//    fout << "$EndElements\n$NodeData\n";
//    fout << "1\n\"" << dataLabel << "\"\n";
//    fout << "1\n0.\n";
//    fout << "3\n0\n" << gridFunction.componentCount() << '\n';
//    fout << nodeCount << '\n' << data.str();
//    fout << "$EndNodeData\n";
//}


#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
template GridFunction<BASIS, RESULT> gridFunctionFromGmsh( \
                const shared_ptr<const Context<BASIS, RESULT> >& context, \
                const GmshIo& gmshIo, shared_ptr<const Grid> grid, \
                GmshPostData::Type gmshPostData, \
                int index = 0); \
template void exportToGmsh(GridFunction<BASIS, RESULT> gridFunction, \
                      const char* dataLabel, const char* fileName, \
                      GmshPostData::Type gmshPostDataType, std::string complexMode); \
template void exportToGmsh(GridFunction<BASIS, RESULT> gridFunction, \
                      const char* dataLabel, \
                      GmshIo& gmshIo, GmshPostData::Type gmshPostDataType, \
                      std::string complexMode)


FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FREE_FUNCTIONS)


}



