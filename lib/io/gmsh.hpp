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

#ifndef gmsh_hpp
#define gmsh_hpp

#include <cstddef>
#include <vector>
#include <iostream>
#include <string>
#include <memory>

namespace bempp {

struct NodeSet {

    struct Node {
        int index;
        double coordinates[3];
    };

    std::vector<Node> nodes;

    void write(std::ostream& output) const;
    static NodeSet read(std::istream& input);

};

struct ElementSet {

    struct Element {
        int index;
        int elementType;
        std::vector<int> tags;
        std::vector<int> nodes;
    };

    std::vector<Element> elements;

    void write(std::ostream& output) const;
    static ElementSet read(std::istream& input);

};

struct PeriodicSet {

    struct PeriodicEntity {
        int dimension;
        int slaveEntityTag;
        int masterEntityTag;
    };

    struct PeriodicNode {
        int slaveNode;
        int masterNode;
    };

    std::vector<PeriodicEntity> periodicEntities;
    std::vector<PeriodicNode> periodicNodes;

    void write(std::ostream& output) const;
    static PeriodicSet read(std::istream& input);

};

struct PhysicalNamesSet {

    struct PhysicalName {
        int dimension;
        int number;
        std::string name;
    };

    std::vector<PhysicalName> physicalNames;

    void write(std::ostream& output) const;
    static PhysicalNamesSet read(std::istream& input);

};

struct NodeDataSet {

    struct NodeValue {
        int index;
        std::vector<double> values;
    };

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    std::vector<NodeValue> nodeValues;

    void write(std::ostream& output) const;
    static NodeDataSet read(std::istream& input);

};

struct ElementDataSet {

    struct ElementValue {
        int index;
        std::vector<double> values;
    };

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    std::vector<ElementValue> elementValues;

    void write(std::ostream& output) const;
    static ElementDataSet read(std::istream& input);

};

struct ElementNodeDataSet {

    struct ElementNodeValue {
        int index;
        std::vector<std::vector<double> > values;
    };

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    std::vector<ElementNodeValue> elementNodeValues;

    void write(std::ostream& output) const;
    static ElementNodeDataSet read(std::istream& input);

};

struct InterpolationSchemeSet {

    struct InterpolationMatrix {
        int nrows;
        int ncols;
        std::vector<double> values; // Stored in row-major order
    };

    std::string name;
    int topology;
    std::vector<InterpolationMatrix> interpolationMatrices;

    void write(std::ostream& output) const;
    static InterpolationSchemeSet read(std::istream& input);

};

struct GmshFile {

    NodeSet nodeSet;
    ElementSet elementSet;
    PeriodicSet periodicSet;
    PhysicalNamesSet physicalNamesSet;
    std::vector<NodeDataSet> nodeDataSets;
    std::vector<ElementDataSet> elementDataSets;
    std::vector<ElementNodeDataSet> elementNodeDataSets;
    std::vector<InterpolationSchemeSet> interpolationSchemeSets;

    void write(std::ostream& output) const;
    void write(std::string fname) const;
    static GmshFile read(std::istream& input);
    static GmshFile read(std::string fname);

};



} // namespace
#endif
