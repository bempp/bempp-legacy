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
                    int capacity, int numberOfFieldComponents, int timeStep, int partition) {

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
                    std::vector<int>& nodeIndices, std::vector<std::vector<double> >& values) {

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


    output << "$MeshFormat";
    output << "2.2" << " " <<
              0 << " " <<
              sizeof(double) << std::endl;
    output << "$EndMeshFormat";

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
void GmshData::write(const std::string& fname) const {

    std::ofstream out;
    out.open(fname.c_str(), std::ios::trunc);
    write(out);
    out.close();

}


GmshData GmshData::read(std::istream& input) {

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
            if (tokens[0] != "2.2") throw std::runtime_error(
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
                    int elementType = boost::lexical_cast<int>(tokens.at(1));
                    int ntags = boost::lexical_cast<int>(tokens.at(2));
                    int physicalEntity = boost::lexical_cast<int>(tokens.at(3));
                    int elementaryEntity = boost::lexical_cast<int>(tokens.at(4));
                    int npartitions = 0;
                    if (ntags > 5) npartitions = boost::lexical_cast<int>(tokens.at(5));
                    std::vector<int> partitions;
                    for (int i = 0;i < npartitions; ++i)
                        partitions.push_back(boost::lexical_cast<int>(tokens.at(6+i)));
                    std::vector<int> nodes;
                    for (int i = 3 + ntags; i < tokens.size(); ++i)
                        nodes.push_back(boost::lexical_cast<int>(tokens.at(i)));
                    gmshData.addElement(index,elementType,nodes,physicalEntity,elementaryEntity,partitions);
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
                gmshData.addElementData(dataSetIndex,index,values);
            }

            std::getline(input,line);
            if (line != "$EndElementData") throw std::runtime_error(
                        "GmshData::read(): Error reading ElementData section.");


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

GmshData GmshData::read(const std::string& fname) {

    std::ifstream(input);
    input.open(fname.c_str());
    GmshData gmshData = read(input);
    input.close();
    return gmshData;

}

} // namespace

//MeshFormat::MeshFormat() : versionNumber("2.2"),
//    fileType(0),
//    dataSize(sizeof(double)) {}

//void MeshFormat::write(std::ostream& output) const {

//    output << "$MeshFormat";
//    output << versionNumber << " " <<
//              fileType << " " <<
//              dataSize << std::endl;
//    output << "$EndMeshFormat";
//}


//MeshFormat MeshFormat::read(std::istream& input) {

//    MeshFormat meshFormat;
//    std::string line;
//    std::getline(input,line);
//    StringVector tokens = stringTokens(line);
//    if (tokens.size() != 3) throw std::runtime_error(
//                "MeshFormat::read(): Wrong format of MeshFormat");

//    meshFormat.versionNumber = tokens[0];
//    if (meshFormat.versionNumber != "2.2") throw std::runtime_error(
//                "MeshFormat::read(): Version of MSH file not supported.");
//    meshFormat.fileType = boost::lexical_cast<int>(tokens[1]);
//    if (meshFormat.fileType != 0) throw std::runtime_error(
//                "MeshFormat::read(): File Type not supported.");
//    meshFormat.dataSize = boost::lexical_cast<int>(tokens[2]);
//    if (meshFormat.dataSize != sizeof(double)) throw std::runtime_error(
//                "MeshFormat::read(): Data size not supported.");
//    return meshFormat;

//}

//void NodeSet::write(std::ostream& output) const{

//    int nnodes = nodes.size();
//    if (nnodes >0){
//        output << "$Nodes" << std::endl;
//        output << nnodes << std::endl;
//        for (int i=0; i < nnodes;i++){
//            output << nodes[i].index << " " <<
//                      boost::lexical_cast<std::string>(nodes[i].coordinates[0]) << " " <<
//                      boost::lexical_cast<std::string>(nodes[i].coordinates[1]) << " " <<
//                      boost::lexical_cast<std::string>(nodes[i].coordinates[2]) << " " << std::endl;
//        }
//        output << "$EndNodes" << std::endl;
//    }
//}

//NodeSet NodeSet::read(std::istream& input) {

//    NodeSet nodeSet;
//    std::string line;
//    std::getline(input,line);
//    int numberOfNodes = boost::lexical_cast<int>(line);
//    nodeSet.nodes.resize(numberOfNodes);

//    int count = 0;
//    if (numberOfNodes>0) {
//        while (std::getline(input,line)){
//            if (count>=numberOfNodes) break;
//            NodeSet::Node& node = nodeSet.nodes[count];
//            StringVector tokens = stringTokens(line);
//            if (tokens.size() != 4) throw std::runtime_error(
//                        "NodeSet::read(): Wrong format of node definition detected.");
//            node.index = boost::lexical_cast<int>(tokens[0]);
//            node.coordinates[0] = boost::lexical_cast<double>(tokens[1]);
//            node.coordinates[1] = boost::lexical_cast<double>(tokens[2]);
//            node.coordinates[2] = boost::lexical_cast<double>(tokens[3]);
//            ++count;
//        }
//    }
//    if (count < numberOfNodes)
//        throw std::runtime_error("NodeSet::read(): Wrong number of nodes.");
//    return nodeSet;
//}

//void ElementSet::write(std::ostream &output) const {

//    int nelements = elements.size();
//    if (nelements>0) {
//        output << "$Elements" << std::endl;
//        output << nelements << std::endl;
//        for (int i = 0; i<nelements; i++) {
//            output << elements[i].index << " " <<
//                      elements[i].elementType << " " <<
//                      elements[i].tags.size();
//            for (int j = 0; j < elements[i].tags.size(); j++)
//                output << " " <<  elements[i].tags[j];
//            for (int j = 0; j < elements[i].nodes.size(); j++)
//                output << " " << elements[i].nodes[j];
//            output << std::endl;
//        }
//        output << "$EndElements" << std::endl;
//    }
//}

//ElementSet ElementSet::read(std::istream &input) {

//    ElementSet elementSet;
//    std::string line;
//    std::getline(input,line);
//    int numberOfElements = boost::lexical_cast<int>(line);

//    elementSet.elements.resize(numberOfElements);
//    if (numberOfElements>0){
//        int count = 0;
//        while (std::getline(input,line)) {
//            if (count >= numberOfElements) break;
//            ElementSet::Element& element = elementSet.elements[count];
//            StringVector tokens = stringTokens(line);
//            element.index = boost::lexical_cast<int>(tokens.at(0));
//            element.elementType = boost::lexical_cast<int>(tokens.at(1));
//            int ntags = boost::lexical_cast<int>(tokens.at(2));
//            element.tags.resize(ntags);
//            for (int i = 2;i < 2 + ntags; ++i)
//                element.tags[i] = boost::lexical_cast<int>(tokens.at(i));
//            for (int i = 2 + ntags; i < tokens.size(); ++i)
//                element.nodes.push_back(boost::lexical_cast<int>(tokens.at(i)));
//            ++count;
//        }
//        if (count < numberOfElements)
//            throw std::runtime_error("ElementSet::read(): Wrong number of elements.");
//        return elementSet;

//    }

//}

//void PeriodicSet::write(std::ostream &output) const {

//    if ((periodicEntities.size() == 0) && (periodicNodes.size() == 0) )
//        return;

//    output << "$Periodic" << std::endl;
//    output << periodicEntities.size() << std::endl;
//    for (int i = 0; i < periodicEntities.size(); ++i){
//        output << periodicEntities[i].dimension << " " <<
//                  periodicEntities[i].slaveEntityTag << " " <<
//                  periodicEntities[i].masterEntityTag << std::endl;
//    }
//    output << periodicNodes.size() << std::endl;
//    for (int i = 0; i < periodicNodes.size(); ++i){
//        output << periodicNodes[i].slaveNode << " " <<
//                  periodicNodes[i].masterNode << std::endl;
//    }
//    output << "$EndPeriodic" << std::endl;

//}

//PeriodicSet PeriodicSet::read(std::istream &input) {

//    PeriodicSet periodicSet;
//    std::string line;
//    std::getline(input,line);
//    int numberOfPeriodicEntities = boost::lexical_cast<int>(line);
//    periodicSet.periodicEntities.resize(numberOfPeriodicEntities);
//    if (numberOfPeriodicEntities > 0) {
//        int count = 0;
//        while (std::getline(input,line)) {
//            if (count >= numberOfPeriodicEntities) break;
//            PeriodicSet::PeriodicEntity& periodicEntity = periodicSet.periodicEntities[count];
//            StringVector tokens = stringTokens(line);
//            if (tokens.size()!= 3) throw std::runtime_error(
//                        "PeriodicSet::read(): Wrong format for periodic entities.");
//            periodicEntity.dimension = boost::lexical_cast<int>(tokens[0]);
//            periodicEntity.slaveEntityTag = boost::lexical_cast<int>(tokens[1]);
//            periodicEntity.masterEntityTag = boost::lexical_cast<int>(tokens[2]);
//            ++count;
//        }
//        if (count<numberOfPeriodicEntities) throw std::runtime_error(
//                    "PeriodicSet::read(): Wrong number of periodic entities.");
//    }
//    std::getline(input,line);
//    int numberOfPeriodicNodes = boost::lexical_cast<int>(line);
//    periodicSet.periodicNodes.resize(numberOfPeriodicNodes);
//    if (numberOfPeriodicNodes > 0) {
//        int count = 0;
//        while (std::getline(input,line)) {
//            if (count >= numberOfPeriodicNodes) break;
//            PeriodicSet::PeriodicNode& periodicNode = periodicSet.periodicNodes[count];
//            StringVector tokens = stringTokens(line);
//            if (tokens.size() != 2) throw std::runtime_error(
//                        "PeriodicSet::read(): Wrong format for periodic nodes.");
//            periodicNode.slaveNode = boost::lexical_cast<int>(tokens[0]);
//            periodicNode.masterNode = boost::lexical_cast<int>(tokens[1]);
//            ++count;
//        }
//        if (count<numberOfPeriodicEntities) throw std::runtime_error(
//                    "PeriodicSet::read(): Wrong number of periodic nodes.");

//    }
//    return periodicSet;
//}

//void PhysicalNamesSet::write(std::ostream &output) const {

//    if (physicalNames.size() == 0) return;
//    output << "$PhysicalNames" << std::endl;
//    for (int i = 0;i < physicalNames.size(); ++i) {
//        output << physicalNames[i].dimension << " " <<
//                  physicalNames[i].number << " " <<
//                  physicalNames[i].name << std::endl;
//    }
//    output << "$EndPhysicalNames" << std::endl;
//}

//PhysicalNamesSet PhysicalNamesSet::read(std::istream &input) {

//    PhysicalNamesSet physicalNamesSet;
//    std::string line;
//    std::getline(input,line);
//    int numberOfPhysicalNames = boost::lexical_cast<int>(line);
//    physicalNamesSet.physicalNames.resize(numberOfPhysicalNames);
//    int count = 0;
//    while(std::getline(input,line)) {
//        if (count >= numberOfPhysicalNames) break;
//        StringVector tokens = stringTokens(line);
//        if (tokens.size() != 3) throw std::runtime_error(
//                    "PhysicalNamesSet::read(): Wrong format for physical names.");
//        PhysicalNamesSet::PhysicalName& physicalName = physicalNamesSet.physicalNames[count];
//        physicalName.dimension = boost::lexical_cast<int>(tokens[0]);
//        physicalName.number = boost::lexical_cast<int>(tokens[1]);
//        physicalName.name = tokens[2];
//        ++count;
//    }
//    if (count<numberOfPhysicalNames) throw std::runtime_error(
//                "PhysicalNamesSet::read(): Wrong number of physical names.");
//    return physicalNamesSet;
//}

//void NodeDataSet::write(std::ostream &output) const {

//    if (nodeValues.size() == 0) return;
//    output << "$NodeData" << std::endl;
//    output << stringTags.size() << std::endl;
//    for (int i = 0; i< stringTags.size(); ++i) {
//        output << '\"'+stringTags[i]+'\"' << std::endl;
//    }
//    output << realTags.size() << std::endl;
//    for (int i = 0; i < realTags.size(); ++i) {
//        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
//    }
//    output << integerTags.size() << std::endl;
//    for (int i = 0; i < integerTags.size(); ++i) {
//        output << integerTags[i] << std::endl;
//    }
//    for (int i = 0; i < nodeValues.size(); ++i ) {
//        output << nodeValues[i].index;
//        for (int j = 0; j < nodeValues[i].values.size(); ++j) {
//            output << " " << boost::lexical_cast<std::string>(nodeValues[i].values[j]);
//        }
//        output << std::endl;
//    }
//    output << "$EndNodeData" << std::endl;
//}

//NodeDataSet NodeDataSet::read(std::istream &input) {

//    NodeDataSet nodeDataSet;
//    std::string line;

//    // String tags
//    std::getline(input,line);
//    int numberOfStringTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfStringTags; ++i) {
//        std::getline(input,line);
//        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
//        nodeDataSet.stringTags.push_back(line);
//    }

//    // Real tags
//    std::getline(input,line);
//    int numberOfRealTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfRealTags; ++i ) {
//        std::getline(input,line);
//        nodeDataSet.realTags.push_back(boost::lexical_cast<double>(line));
//    }

//    // Integer tags
//    std::getline(input,line);
//    int numberOfIntegerTags = boost::lexical_cast<int>(line);
//    if (numberOfIntegerTags < 3) throw std::runtime_error(
//                "NodeDataSet::read(): At least 3 integer tags required.");
//    for (int i = 0; i < numberOfIntegerTags; ++i) {
//        std::getline(input,line);
//        nodeDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
//    }
//    int numberOfNodes = nodeDataSet.integerTags[2];
//    int valuesPerNode = nodeDataSet.integerTags[1];
//    nodeDataSet.nodeValues.resize(numberOfNodes);
//    for (int i = 0; i < numberOfNodes; ++i) {
//        std::getline(input,line);
//        NodeDataSet::NodeValue& nodeValue = nodeDataSet.nodeValues[i];
//        StringVector tokens = stringTokens(line);
//        if (tokens.size() != 1+valuesPerNode) throw std::runtime_error(
//                    "NodesDataSet::read(): Data has wrong format.");
//        nodeValue.index = boost::lexical_cast<int>(tokens[0]);
//        for (int j = 1; j < 1+valuesPerNode; ++j) {
//            nodeValue.values.push_back(boost::lexical_cast<double>(tokens[j]));
//        }
//    }
//    return nodeDataSet;
//}

//void ElementDataSet::write(std::ostream &output) const {

//    if (elementValues.size() == 0) return;
//    output << "$ElementData" << std::endl;
//    output << stringTags.size() << std::endl;
//    for (int i = 0; i< stringTags.size(); ++i) {
//        output << '\"'+stringTags[i]+'\"' << std::endl;
//    }
//    output << realTags.size() << std::endl;
//    for (int i = 0; i < realTags.size(); ++i) {
//        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
//    }
//    output << integerTags.size() << std::endl;
//    for (int i = 0; i < integerTags.size(); ++i) {
//        output << integerTags[i] << std::endl;
//    }
//    for (int i = 0; i < elementValues.size(); ++i ) {
//        output << elementValues[i].index;
//        for (int j = 0; j < elementValues[i].values.size(); ++j) {
//            output << " " << boost::lexical_cast<std::string>(elementValues[i].values[j]);
//        }
//        output << std::endl;
//    }
//    output << "$EndElementData" << std::endl;
//}

//ElementDataSet ElementDataSet::read(std::istream &input) {

//    ElementDataSet elementDataSet;
//    std::string line;

//    // String tags
//    std::getline(input,line);
//    int numberOfStringTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfStringTags; ++i) {
//        std::getline(input,line);
//        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
//        elementDataSet.stringTags.push_back(line);
//    }

//    // Real tags
//    std::getline(input,line);
//    int numberOfRealTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfRealTags; ++i ) {
//        std::getline(input,line);
//        elementDataSet.realTags.push_back(boost::lexical_cast<double>(line));
//    }

//    // Integer tags
//    std::getline(input,line);
//    int numberOfIntegerTags = boost::lexical_cast<int>(line);
//    if (numberOfIntegerTags < 3) throw std::runtime_error(
//                "ElementDataSet::read(): At least 3 integer tags required.");
//    for (int i = 0; i < numberOfIntegerTags; ++i) {
//        std::getline(input,line);
//        elementDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
//    }
//    int numberOfElements = elementDataSet.integerTags[2];
//    int valuesPerElement = elementDataSet.integerTags[1];
//    elementDataSet.elementValues.resize(numberOfElements);
//    for (int i = 0; i < numberOfElements; ++i) {
//        std::getline(input,line);
//        ElementDataSet::ElementValue& elementValue = elementDataSet.elementValues[i];
//        StringVector tokens = stringTokens(line);
//        if (tokens.size() != 1+valuesPerElement) throw std::runtime_error(
//                    "ElementDataSet::read(): Data has wrong format.");
//        elementValue.index = boost::lexical_cast<int>(tokens[0]);
//        for (int j = 1; j < 1+valuesPerElement; ++j) {
//            elementValue.values.push_back(boost::lexical_cast<double>(tokens[j]));
//        }
//    }
//    return elementDataSet;
//}

//void ElementNodeDataSet::write(std::ostream &output) const {

//    if (elementNodeValues.size() == 0) return;
//    output << "$ElementNodeData" << std::endl;
//    output << stringTags.size() << std::endl;
//    for (int i = 0; i< stringTags.size(); ++i) {
//        output << '\"'+stringTags[i]+'\"' << std::endl;
//    }
//    output << realTags.size() << std::endl;
//    for (int i = 0; i < realTags.size(); ++i) {
//        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
//    }
//    output << integerTags.size() << std::endl;
//    for (int i = 0; i < integerTags.size(); ++i) {
//        output << integerTags[i] << std::endl;
//    }
//    for (int i = 0; i < elementNodeValues.size(); ++i ) {
//        output << elementNodeValues[i].index << " " <<
//                  elementNodeValues[i].values.size();
//        for (int j = 0; j < elementNodeValues[i].values.size(); ++j) {
//            for (int k = 0; k < elementNodeValues[i].values[j].size(); ++k)
//                output << " " << boost::lexical_cast<std::string>(elementNodeValues[i].values[j][k]);
//        }
//        output << std::endl;
//    }
//    output << "$EndElementNodeData" << std::endl;
//}

//ElementNodeDataSet ElementNodeDataSet::read(std::istream &input) {

//    ElementNodeDataSet elementNodeDataSet;
//    std::string line;

//    // String tags
//    std::getline(input,line);
//    int numberOfStringTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfStringTags; ++i) {
//        std::getline(input,line);
//        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
//        elementNodeDataSet.stringTags.push_back(line);
//    }

//    // Real tags
//    std::getline(input,line);
//    int numberOfRealTags = boost::lexical_cast<int>(line);
//    for (int i = 0; i < numberOfRealTags; ++i ) {
//        std::getline(input,line);
//        elementNodeDataSet.realTags.push_back(boost::lexical_cast<double>(line));
//    }

//    // Integer tags
//    std::getline(input,line);
//    int numberOfIntegerTags = boost::lexical_cast<int>(line);
//    if (numberOfIntegerTags < 3) throw std::runtime_error(
//                "ElementNodeDataSet::read(): At least 3 integer tags required.");
//    for (int i = 0; i < numberOfIntegerTags; ++i) {
//        std::getline(input,line);
//        elementNodeDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
//    }
//    int numberOfElements = elementNodeDataSet.integerTags[2];
//    int valuesPerNode = elementNodeDataSet.integerTags[1];
//    elementNodeDataSet.elementNodeValues.resize(numberOfElements);
//    for (int i = 0; i < numberOfElements; ++i) {
//        std::getline(input,line);
//        ElementNodeDataSet::ElementNodeValue& elementNodeValue = elementNodeDataSet.elementNodeValues[i];
//        StringVector tokens = stringTokens(line);
//        if (tokens.size() < 2) throw std::runtime_error(
//                    "ElementNodeDataSet::read(): Data has wrong format.");
//        elementNodeValue.index = boost::lexical_cast<int>(tokens[0]);
//        int numberOfNodes = boost::lexical_cast<int>(tokens[1]);
//        if (tokens.size() != 2+numberOfNodes*valuesPerNode) throw std::runtime_error(
//                    "ElementNodeDataSet::read(): Data has wrong format.");
//        elementNodeValue.values.resize(numberOfNodes);
//        for (int j = 0; j < numberOfNodes; ++j) {
//            elementNodeValue.values[j].resize(valuesPerNode);
//            for (int k = 0; k < valuesPerNode; ++k)
//                elementNodeValue.values[j][k] = boost::lexical_cast<double>(tokens[j]);
//        }
//    }
//    return elementNodeDataSet;
//}

//void InterpolationSchemeSet::write(std::ostream &output) const {

//    if (interpolationMatrices.size() == 0) return;
//    output << "$InterpolationScheme" << std::endl;
//    output << '\"'+name+'\"' << std::endl;
//    output << 1 << std::endl; // Only one element topology supported right now.
//    output << topology << std::endl;
//    output << interpolationMatrices.size() << std::endl;
//    for (int i = 0; i < interpolationMatrices.size(); ++i) {
//        output << interpolationMatrices[i].nrows << std::endl;
//        output << interpolationMatrices[i].ncols << std::endl;
//        for (int j = 0; j < interpolationMatrices[i].nrows; ++j) {
//            output << interpolationMatrices[i].values[j*interpolationMatrices[i].nrows];
//            for (int k = 1; k < interpolationMatrices[i].ncols; ++k)
//                output << " " << interpolationMatrices[i].values[j*interpolationMatrices[i].nrows+k];
//            output << std::endl;
//        }
//    }
//    output << "$EndInterpolationScheme" << std::endl;
//}

//InterpolationSchemeSet InterpolationSchemeSet::read(std::istream &input) {

//    InterpolationSchemeSet interpolationSchemeSet;
//    std::string line;
//    std::getline(input,line);
//    line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
//    interpolationSchemeSet.name = line;
//    std::getline(input,line);
//    if (boost::lexical_cast<int>(line) != 1) throw std::runtime_error(
//                "InterpolationSchemSet::read(): Only one topology is currently supported.");
//    std::getline(input,line);
//    interpolationSchemeSet.topology = boost::lexical_cast<int>(line);
//    std::getline(input,line);
//    interpolationSchemeSet.interpolationMatrices.resize(boost::lexical_cast<int>(line));
//    for (int i = 0; i < interpolationSchemeSet.interpolationMatrices.size(); ++i) {
//        InterpolationSchemeSet::InterpolationMatrix& interpolationMatrix = interpolationSchemeSet.interpolationMatrices[i];
//        std::getline(input,line);
//        StringVector tokens = stringTokens(line);
//        if (tokens.size() != 2) throw std::runtime_error(
//                    "InterpolationSchemeSet::read(): Wrong format of interpolation matrices.");
//        interpolationMatrix.nrows = boost::lexical_cast<int>(tokens[0]);
//        interpolationMatrix.ncols = boost::lexical_cast<int>(tokens[1]);
//        interpolationMatrix.values.reserve(interpolationMatrix.nrows*interpolationMatrix.ncols);
//        for (int j = 0; j < interpolationMatrix.nrows; ++j) {
//            std::getline(input,line);
//            StringVector tokens = stringTokens(line);
//            if (tokens.size() != interpolationMatrix.ncols) throw std::runtime_error(
//                        "InterpolationSchemeSet::read(): Wrong format of interpolation matrices.");
//            for (int k = 0;k < interpolationMatrix.ncols; ++k)
//                interpolationMatrix.values.push_back(boost::lexical_cast<double>(tokens[k]));
//        }
//    }
//    return interpolationSchemeSet;
//}

//void GmshFile::write(std::ostream &output) const {

//    meshFormat.write(output);
//    nodeSet.write(output);
//    elementSet.write(output);
//    periodicSet.write(output);
//    physicalNamesSet.write(output);
//    for (int i = 0; i < nodeDataSets.size(); ++i)
//        nodeDataSets[i].write(output);
//    for (int i = 0; i < elementDataSets.size(); ++i)
//        elementDataSets[i].write(output);
//    for (int i = 0; i < elementNodeDataSets.size(); ++i)
//        elementNodeDataSets[i].write(output);
//    for (int i = 0; i < interpolationSchemeSets.size(); ++i)
//        interpolationSchemeSets[i].write(output);
//}

//void GmshFile::write(std::string fname) const {

//    std::ofstream out;
//    out.open(fname.c_str(), std::ios::trunc);
//    write(out);
//    out.close();

//}

//GmshFile GmshFile::read(std::istream& input) {

//    GmshFile gmshFile;
//    std::string line;

//    bool haveMeshFormat = false;
//    bool haveNodes = false;
//    bool haveElements = false;
//    bool havePeriodic = false;
//    bool havePhysicalNames = false;

//    while (std::getline(input,line)) {
//        if (line == "$MeshFormat") {
//            if (haveMeshFormat) throw std::runtime_error(
//                        "GmshFile::read(): '$MeshFormat' section can only appear once.");
//            gmshFile.meshFormat = MeshFormat::read(input);
//            std::getline(input,line);
//            if (line != "$EndMeshFormat") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$MeshFormat' section");
//            haveMeshFormat = true;
//        }
//        else if (line == "$Nodes") {
//            if (haveNodes) throw std::runtime_error(
//                        "GmshFile::read(): '$Nodes' section can only appear once.");
//            gmshFile.nodeSet = NodeSet::read(input);
//            std::getline(input,line);
//            if (line != "$EndNodes") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$Nodes' section");
//            haveNodes = true;
//        }
//        else if (line == "$Elements") {
//            if (haveElements) throw std::runtime_error(
//                        "GmshFile::read(): '$Elements' section can only appear once.");
//            gmshFile.elementSet = ElementSet::read(input);
//            std::getline(input,line);
//            if (line != "$EndElements") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$Elements' section.");
//            haveElements = true;
//        }
//        else if (line == "$Periodic") {
//            if (havePeriodic) throw std::runtime_error(
//                        "GmshFile::read(): '$Periodic' section can only appear once.");
//            gmshFile.periodicSet = PeriodicSet::read(input);
//            std::getline(input,line);
//            if (line != "$EndPeriodic") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$Periodic' section.");
//            havePeriodic = true;
//        }
//        else if (line == "$PhysicalNames") {
//            if (havePhysicalNames) throw std::runtime_error(
//                        "GmshFile::read(): '$PhysicalNames' section can only appear once.");
//            gmshFile.physicalNamesSet = PhysicalNamesSet::read(input);
//            std::getline(input,line);
//            if (line != "$EndPhysicalNames") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$PhysicalNames' section.");
//            havePhysicalNames = true;
//        }
//        else if (line == "$NodeData") {
//            gmshFile.nodeDataSets.push_back(NodeDataSet::read(input));
//            std::getline(input,line);
//            if (line != "$EndNodeData") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$NodeData' section.");
//        }
//        else if (line == "$ElementData") {
//            gmshFile.elementDataSets.push_back(ElementDataSet::read(input));
//            std::getline(input,line);
//            if (line != "$EndElementData") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$ElementData' section.");
//        }
//        else if (line == "$ElementNodeData") {
//            gmshFile.elementNodeDataSets.push_back(ElementNodeDataSet::read(input));
//            std::getline(input,line);
//            if (line != "$EndElementNodeData") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$ElementNodeData' section. ");
//        }
//        else if (line == "$InterpolationScheme") {
//            gmshFile.interpolationSchemeSets.push_back(InterpolationSchemeSet::read(input));
//            if (line != "$EndInterpolationScheme") throw std::runtime_error(
//                        "GmshFile::read(): Error reading '$InterpolationScheme' section.");
//        }
//    }
//    return gmshFile;
//}

//GmshFile GmshFile::read(std::string fname) {

//    std::ifstream in;
//    in.open(fname.c_str());
//    GmshFile gmshFile = read(in);
//    in.close();
//    return gmshFile;

//}


