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
    boost::tokenizer<> tok(input);
    boost::tokenizer<>::iterator it;
    for (it = tok.begin(); it != tok.end(); ++it)
        result.push_back(*it);
    return result;
}

} // namespace

namespace bempp {

void NodeSet::write(std::ostream& output) const{

    int nnodes = nodes.size();
    if (nnodes >0){
        output << "$Nodes" << std::endl;
        output << nnodes << std::endl;
        for (int i=0; i < nnodes;i++){
            output << nodes[i].index << " " <<
                      boost::lexical_cast<std::string>(nodes[i].coordinates[0]) << " " <<
                      boost::lexical_cast<std::string>(nodes[i].coordinates[1]) << " " <<
                      boost::lexical_cast<std::string>(nodes[i].coordinates[2]) << " " << std::endl;
        }
        output << "$EndNodes" << std::endl;
    }
}

NodeSet NodeSet::read(std::istream& input) {

    NodeSet nodeSet;
    std::string line;
    std::getline(input,line);
    int numberOfNodes = boost::lexical_cast<int>(line);
    nodeSet.nodes.resize(numberOfNodes);

    int count = 0;
    if (numberOfNodes>0) {
        while (std::getline(input,line)){
            if (count>=numberOfNodes) break;
            NodeSet::Node& node = nodeSet.nodes[count];
            StringVector tokens = stringTokens(line);
            if (tokens.size() != 4) throw std::runtime_error(
                        "NodeSet::read(): Wrong format of node definition detected.");
            node.index = boost::lexical_cast<int>(tokens[0]);
            node.coordinates[0] = boost::lexical_cast<double>(tokens[1]);
            node.coordinates[1] = boost::lexical_cast<double>(tokens[2]);
            node.coordinates[2] = boost::lexical_cast<double>(tokens[3]);
            ++count;
        }
    }
    if (count < numberOfNodes)
        throw std::runtime_error("NodeSet::read(): Wrong number of nodes.");
    return nodeSet;
}

void ElementSet::write(std::ostream &output) const {

    int nelements = elements.size();
    if (nelements>0) {
        output << "$Elements" << std::endl;
        output << nelements << std::endl;
        for (int i = 0; i<nelements; i++) {
            output << elements[i].index << " " <<
                      elements[i].elementType << " " <<
                      elements[i].tags.size();
            for (int j = 0; j < elements[i].tags.size(); j++)
                output << " " <<  elements[i].tags[j];
            for (int j = 0; j < elements[i].nodes.size(); j++)
                output << " " << elements[i].nodes[j];
            output << std::endl;
        }
        output << "$EndElements" << std::endl;
    }
}

ElementSet ElementSet::read(std::istream &input) {

    ElementSet elementSet;
    std::string line;
    std::getline(input,line);
    int numberOfElements = boost::lexical_cast<int>(line);

    elementSet.elements.resize(numberOfElements);
    if (numberOfElements>0){
        int count = 0;
        while (std::getline(input,line)) {
            if (count >= numberOfElements) break;
            ElementSet::Element& element = elementSet.elements[count];
            StringVector tokens = stringTokens(line);
            element.index = boost::lexical_cast<int>(tokens.at(0));
            element.elementType = boost::lexical_cast<int>(tokens.at(1));
            int ntags = boost::lexical_cast<int>(tokens.at(2));
            element.tags.resize(ntags);
            for (int i = 2;i < 2 + ntags; ++i)
                element.tags[i] = boost::lexical_cast<int>(tokens.at(i));
            for (int i = 2 + ntags; i < tokens.size(); ++i)
                element.nodes.push_back(boost::lexical_cast<int>(tokens.at(i)));
            ++count;
        }
        if (count < numberOfElements)
            throw std::runtime_error("ElementSet::read(): Wrong number of elements.");
        return elementSet;

    }

}

void PeriodicSet::write(std::ostream &output) const {

    if ((periodicEntities.size() == 0) && (periodicNodes.size() == 0) )
        return;

    output << "$Periodic" << std::endl;
    output << periodicEntities.size() << std::endl;
    for (int i = 0; i < periodicEntities.size(); ++i){
        output << periodicEntities[i].dimension << " " <<
                  periodicEntities[i].slaveEntityTag << " " <<
                  periodicEntities[i].masterEntityTag << std::endl;
    }
    output << periodicNodes.size() << std::endl;
    for (int i = 0; i < periodicNodes.size(); ++i){
        output << periodicNodes[i].slaveNode << " " <<
                  periodicNodes[i].masterNode << std::endl;
    }
    output << "$EndPeriodic" << std::endl;

}

PeriodicSet PeriodicSet::read(std::istream &input) {

    PeriodicSet periodicSet;
    std::string line;
    std::getline(input,line);
    int numberOfPeriodicEntities = boost::lexical_cast<int>(line);
    periodicSet.periodicEntities.resize(numberOfPeriodicEntities);
    if (numberOfPeriodicEntities > 0) {
        int count = 0;
        while (std::getline(input,line)) {
            if (count >= numberOfPeriodicEntities) break;
            PeriodicSet::PeriodicEntity& periodicEntity = periodicSet.periodicEntities[count];
            StringVector tokens = stringTokens(line);
            if (tokens.size()!= 3) throw std::runtime_error(
                        "PeriodicSet::read(): Wrong format for periodic entities.");
            periodicEntity.dimension = boost::lexical_cast<int>(tokens[0]);
            periodicEntity.slaveEntityTag = boost::lexical_cast<int>(tokens[1]);
            periodicEntity.masterEntityTag = boost::lexical_cast<int>(tokens[2]);
            ++count;
        }
        if (count<numberOfPeriodicEntities) throw std::runtime_error(
                    "PeriodicSet::read(): Wrong number of periodic entities.");
    }
    std::getline(input,line);
    int numberOfPeriodicNodes = boost::lexical_cast<int>(line);
    periodicSet.periodicNodes.resize(numberOfPeriodicNodes);
    if (numberOfPeriodicNodes > 0) {
        int count = 0;
        while (std::getline(input,line)) {
            if (count >= numberOfPeriodicNodes) break;
            PeriodicSet::PeriodicNode& periodicNode = periodicSet.periodicNodes[count];
            StringVector tokens = stringTokens(line);
            if (tokens.size() != 2) throw std::runtime_error(
                        "PeriodicSet::read(): Wrong format for periodic nodes.");
            periodicNode.slaveNode = boost::lexical_cast<int>(tokens[0]);
            periodicNode.masterNode = boost::lexical_cast<int>(tokens[1]);
            ++count;
        }
        if (count<numberOfPeriodicEntities) throw std::runtime_error(
                    "PeriodicSet::read(): Wrong number of periodic nodes.");

    }
    return periodicSet;
}

void PhysicalNamesSet::write(std::ostream &output) const {

    if (physicalNames.size() == 0) return;
    output << "$PhysicalNames" << std::endl;
    for (int i = 0;i < physicalNames.size(); ++i) {
        output << physicalNames[i].dimension << " " <<
                  physicalNames[i].number << " " <<
                  physicalNames[i].name << std::endl;
    }
    output << "$EndPhysicalNames" << std::endl;
}

PhysicalNamesSet PhysicalNamesSet::read(std::istream &input) {

    PhysicalNamesSet physicalNamesSet;
    std::string line;
    std::getline(input,line);
    int numberOfPhysicalNames = boost::lexical_cast<int>(line);
    physicalNamesSet.physicalNames.resize(numberOfPhysicalNames);
    int count = 0;
    while(std::getline(input,line)) {
        if (count >= numberOfPhysicalNames) break;
        StringVector tokens = stringTokens(line);
        if (tokens.size() != 3) throw std::runtime_error(
                    "PhysicalNamesSet::read(): Wrong format for physical names.");
        PhysicalNamesSet::PhysicalName& physicalName = physicalNamesSet.physicalNames[count];
        physicalName.dimension = boost::lexical_cast<int>(tokens[0]);
        physicalName.number = boost::lexical_cast<int>(tokens[1]);
        physicalName.name = tokens[2];
        ++count;
    }
    if (count<numberOfPhysicalNames) throw std::runtime_error(
                "PhysicalNamesSet::read(): Wrong number of physical names.");
    return physicalNamesSet;
}

void NodeDataSet::write(std::ostream &output) const {

    if (nodeValues.size() == 0) return;
    output << "$NodeData" << std::endl;
    output << stringTags.size() << std::endl;
    for (int i = 0; i< stringTags.size(); ++i) {
        output << '\"'+stringTags[i]+'\"' << std::endl;
    }
    output << realTags.size() << std::endl;
    for (int i = 0; i < realTags.size(); ++i) {
        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
    }
    output << integerTags.size() << std::endl;
    for (int i = 0; i < integerTags.size(); ++i) {
        output << integerTags[i] << std::endl;
    }
    for (int i = 0; i < nodeValues.size(); ++i ) {
        output << nodeValues[i].index;
        for (int j = 0; j < nodeValues[i].values.size(); ++j) {
            output << " " << boost::lexical_cast<std::string>(nodeValues[i].values[j]);
        }
        output << std::endl;
    }
    output << "$EndNodeData" << std::endl;
}

NodeDataSet NodeDataSet::read(std::istream &input) {

    NodeDataSet nodeDataSet;
    std::string line;

    // String tags
    std::getline(input,line);
    int numberOfStringTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfStringTags; ++i) {
        std::getline(input,line);
        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
        nodeDataSet.stringTags.push_back(line);
    }

    // Real tags
    std::getline(input,line);
    int numberOfRealTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfRealTags; ++i ) {
        std::getline(input,line);
        nodeDataSet.realTags.push_back(boost::lexical_cast<double>(line));
    }

    // Integer tags
    std::getline(input,line);
    int numberOfIntegerTags = boost::lexical_cast<int>(line);
    if (numberOfIntegerTags < 3) throw std::runtime_error(
                "NodeDataSet::read(): At least 3 integer tags required.");
    for (int i = 0; i < numberOfIntegerTags; ++i) {
        std::getline(input,line);
        nodeDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
    }
    int numberOfNodes = nodeDataSet.integerTags[2];
    int valuesPerNode = nodeDataSet.integerTags[1];
    nodeDataSet.nodeValues.resize(numberOfNodes);
    for (int i = 0; i < numberOfNodes; ++i) {
        std::getline(input,line);
        NodeDataSet::NodeValue& nodeValue = nodeDataSet.nodeValues[i];
        StringVector tokens = stringTokens(line);
        if (tokens.size() != 1+valuesPerNode) throw std::runtime_error(
                    "NodesDataSet::read(): Data has wrong format.");
        nodeValue.index = boost::lexical_cast<int>(tokens[0]);
        for (int j = 1; j < 1+valuesPerNode; ++j) {
            nodeValue.values.push_back(boost::lexical_cast<double>(tokens[j]));
        }
    }
    return nodeDataSet;
}

void ElementDataSet::write(std::ostream &output) const {

    if (elementValues.size() == 0) return;
    output << "$ElementData" << std::endl;
    output << stringTags.size() << std::endl;
    for (int i = 0; i< stringTags.size(); ++i) {
        output << '\"'+stringTags[i]+'\"' << std::endl;
    }
    output << realTags.size() << std::endl;
    for (int i = 0; i < realTags.size(); ++i) {
        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
    }
    output << integerTags.size() << std::endl;
    for (int i = 0; i < integerTags.size(); ++i) {
        output << integerTags[i] << std::endl;
    }
    for (int i = 0; i < elementValues.size(); ++i ) {
        output << elementValues[i].index;
        for (int j = 0; j < elementValues[i].values.size(); ++j) {
            output << " " << boost::lexical_cast<std::string>(elementValues[i].values[j]);
        }
        output << std::endl;
    }
    output << "$EndElementData" << std::endl;
}

ElementDataSet ElementDataSet::read(std::istream &input) {

    ElementDataSet elementDataSet;
    std::string line;

    // String tags
    std::getline(input,line);
    int numberOfStringTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfStringTags; ++i) {
        std::getline(input,line);
        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
        elementDataSet.stringTags.push_back(line);
    }

    // Real tags
    std::getline(input,line);
    int numberOfRealTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfRealTags; ++i ) {
        std::getline(input,line);
        elementDataSet.realTags.push_back(boost::lexical_cast<double>(line));
    }

    // Integer tags
    std::getline(input,line);
    int numberOfIntegerTags = boost::lexical_cast<int>(line);
    if (numberOfIntegerTags < 3) throw std::runtime_error(
                "ElementDataSet::read(): At least 3 integer tags required.");
    for (int i = 0; i < numberOfIntegerTags; ++i) {
        std::getline(input,line);
        elementDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
    }
    int numberOfElements = elementDataSet.integerTags[2];
    int valuesPerElement = elementDataSet.integerTags[1];
    elementDataSet.elementValues.resize(numberOfElements);
    for (int i = 0; i < numberOfElements; ++i) {
        std::getline(input,line);
        ElementDataSet::ElementValue& elementValue = elementDataSet.elementValues[i];
        StringVector tokens = stringTokens(line);
        if (tokens.size() != 1+valuesPerElement) throw std::runtime_error(
                    "ElementDataSet::read(): Data has wrong format.");
        elementValue.index = boost::lexical_cast<int>(tokens[0]);
        for (int j = 1; j < 1+valuesPerElement; ++j) {
            elementValue.values.push_back(boost::lexical_cast<double>(tokens[j]));
        }
    }
    return elementDataSet;
}

void ElementNodeDataSet::write(std::ostream &output) const {

    if (elementNodeValues.size() == 0) return;
    output << "$ElementNodeData" << std::endl;
    output << stringTags.size() << std::endl;
    for (int i = 0; i< stringTags.size(); ++i) {
        output << '\"'+stringTags[i]+'\"' << std::endl;
    }
    output << realTags.size() << std::endl;
    for (int i = 0; i < realTags.size(); ++i) {
        output << boost::lexical_cast<std::string>(realTags[i]) << std::endl;
    }
    output << integerTags.size() << std::endl;
    for (int i = 0; i < integerTags.size(); ++i) {
        output << integerTags[i] << std::endl;
    }
    for (int i = 0; i < elementNodeValues.size(); ++i ) {
        output << elementNodeValues[i].index << " " <<
                  elementNodeValues[i].values.size();
        for (int j = 0; j < elementNodeValues[i].values.size(); ++j) {
            for (int k = 0; k < elementNodeValues[i].values[j].size(); ++k)
                output << " " << boost::lexical_cast<std::string>(elementNodeValues[i].values[j][k]);
        }
        output << std::endl;
    }
    output << "$EndElementNodeData" << std::endl;
}

ElementNodeDataSet ElementNodeDataSet::read(std::istream &input) {

    ElementNodeDataSet elementNodeDataSet;
    std::string line;

    // String tags
    std::getline(input,line);
    int numberOfStringTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfStringTags; ++i) {
        std::getline(input,line);
        line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
        elementNodeDataSet.stringTags.push_back(line);
    }

    // Real tags
    std::getline(input,line);
    int numberOfRealTags = boost::lexical_cast<int>(line);
    for (int i = 0; i < numberOfRealTags; ++i ) {
        std::getline(input,line);
        elementNodeDataSet.realTags.push_back(boost::lexical_cast<double>(line));
    }

    // Integer tags
    std::getline(input,line);
    int numberOfIntegerTags = boost::lexical_cast<int>(line);
    if (numberOfIntegerTags < 3) throw std::runtime_error(
                "ElementNodeDataSet::read(): At least 3 integer tags required.");
    for (int i = 0; i < numberOfIntegerTags; ++i) {
        std::getline(input,line);
        elementNodeDataSet.integerTags.push_back(boost::lexical_cast<int>(line));
    }
    int numberOfElements = elementNodeDataSet.integerTags[2];
    int valuesPerNode = elementNodeDataSet.integerTags[1];
    elementNodeDataSet.elementNodeValues.resize(numberOfElements);
    for (int i = 0; i < numberOfElements; ++i) {
        std::getline(input,line);
        ElementNodeDataSet::ElementNodeValue& elementNodeValue = elementNodeDataSet.elementNodeValues[i];
        StringVector tokens = stringTokens(line);
        if (tokens.size() < 2) throw std::runtime_error(
                    "ElementNodeDataSet::read(): Data has wrong format.");
        elementNodeValue.index = boost::lexical_cast<int>(tokens[0]);
        int numberOfNodes = boost::lexical_cast<int>(tokens[1]);
        if (tokens.size() != 2+numberOfNodes*valuesPerNode) throw std::runtime_error(
                    "ElementNodeDataSet::read(): Data has wrong format.");
        elementNodeValue.values.resize(numberOfNodes);
        for (int j = 0; j < numberOfNodes; ++j) {
            elementNodeValue.values[j].resize(valuesPerNode);
            for (int k = 0; k < valuesPerNode; ++k)
                elementNodeValue.values[j][k] = boost::lexical_cast<double>(tokens[j]);
        }
    }
    return elementNodeDataSet;
}

void InterpolationSchemeSet::write(std::ostream &output) const {

    if (interpolationMatrices.size() == 0) return;
    output << "$InterpolationScheme" << std::endl;
    output << '\"'+name+'\"' << std::endl;
    output << 1 << std::endl; // Only one element topology supported right now.
    output << topology << std::endl;
    output << interpolationMatrices.size() << std::endl;
    for (int i = 0; i < interpolationMatrices.size(); ++i) {
        output << interpolationMatrices[i].nrows << std::endl;
        output << interpolationMatrices[i].ncols << std::endl;
        for (int j = 0; j < interpolationMatrices[i].nrows; ++j) {
            output << interpolationMatrices[i].values[j*interpolationMatrices[i].nrows];
            for (int k = 1; k < interpolationMatrices[i].ncols; ++k)
                output << " " << interpolationMatrices[i].values[j*interpolationMatrices[i].nrows+k];
            output << std::endl;
        }
    }
    output << "$EndInterpolationScheme" << std::endl;
}

InterpolationSchemeSet InterpolationSchemeSet::read(std::istream &input) {

    InterpolationSchemeSet interpolationSchemeSet;
    std::string line;
    std::getline(input,line);
    line.erase(std::remove(line.begin(),line.end(),'\"'),line.end());
    interpolationSchemeSet.name = line;
    std::getline(input,line);
    if (boost::lexical_cast<int>(line) != 1) throw std::runtime_error(
                "InterpolationSchemSet::read(): Only one topology is currently supported.");
    std::getline(input,line);
    interpolationSchemeSet.topology = boost::lexical_cast<int>(line);
    std::getline(input,line);
    interpolationSchemeSet.interpolationMatrices.resize(boost::lexical_cast<int>(line));
    for (int i = 0; i < interpolationSchemeSet.interpolationMatrices.size(); ++i) {
        InterpolationSchemeSet::InterpolationMatrix& interpolationMatrix = interpolationSchemeSet.interpolationMatrices[i];
        std::getline(input,line);
        StringVector tokens = stringTokens(line);
        if (tokens.size() != 2) throw std::runtime_error(
                    "InterpolationSchemeSet::read(): Wrong format of interpolation matrices.");
        interpolationMatrix.nrows = boost::lexical_cast<int>(tokens[0]);
        interpolationMatrix.ncols = boost::lexical_cast<int>(tokens[1]);
        interpolationMatrix.values.reserve(interpolationMatrix.nrows*interpolationMatrix.ncols);
        for (int j = 0; j < interpolationMatrix.nrows; ++j) {
            std::getline(input,line);
            StringVector tokens = stringTokens(line);
            if (tokens.size() != interpolationMatrix.ncols) throw std::runtime_error(
                        "InterpolationSchemeSet::read(): Wrong format of interpolation matrices.");
            for (int k = 0;k < interpolationMatrix.ncols; ++k)
                interpolationMatrix.values.push_back(boost::lexical_cast<double>(tokens[k]));
        }
    }
    return interpolationSchemeSet;
}

void GmshFile::write(std::ostream &output) const {

    nodeSet.write(output);
    elementSet.write(output);
    periodicSet.write(output);
    physicalNamesSet.write(output);
    for (int i = 0; i < nodeDataSets.size(); ++i)
        nodeDataSets[i].write(output);
    for (int i = 0; i < elementDataSets.size(); ++i)
        elementDataSets[i].write(output);
    for (int i = 0; i < elementNodeDataSets.size(); ++i)
        elementNodeDataSets[i].write(output);
    for (int i = 0; i < interpolationSchemeSets.size(); ++i)
        interpolationSchemeSets[i].write(output);
}

void GmshFile::write(std::string fname) const {

    std::ofstream out;
    out.open(fname.c_str(), std::ios::trunc);
    write(out);
    out.close();

}



} // namespace bempp
