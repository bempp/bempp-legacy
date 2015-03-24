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
#include <set>
#include "../common/shared_ptr.hpp"
#include "../assembly/grid_function.hpp"
#include "../common/eigen_support.hpp"

namespace Bempp {

/** \cond FORWARD_DECL */
class Grid;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
/** \endcond */

struct GmshPostData {

  enum Type {
    NODE,
    ELEMENT,
    ELEMENT_NODE
  };
};

class GmshData {

public:
  enum GmshPostDataType {
    NODE,
    ELEMENT,
    ELEMENT_NODE
  };

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
  void addElement(int index, int elementType, const std::vector<int> &nodes,
                  int physicalEntity, int elementaryEntity,
                  const std::vector<int> &partitions = std::vector<int>());

  void addPeriodicEntity(int dimension, int slaveEntityTag,
                         int masterEntityTag);
  void addPeriodicNode(int slaveNode, int masterNode);

  void addPhysicalName(int dimension, int number, std::string name);

  void addNodeDataSet(const std::vector<std::string> &stringTags,
                      const std::vector<double> &realTags,
                      int numberOfFieldComponents, int capacity = 0,
                      int timeStep = 0, int partition = 0);
  void addElementDataSet(const std::vector<std::string> &stringTags,
                         const std::vector<double> &realTags,
                         int numberOfFieldComponents, int capacity = 0,
                         int timeStep = 0, int partition = 0);
  void addElementNodeDataSet(const std::vector<std::string> &stringTags,
                             const std::vector<double> &realTags,
                             int numberOfFieldComponents, int capacity = 0,
                             int timeStep = 0, int partition = 0);

  void addNodeData(int dataSetIndex, int node,
                   const std::vector<double> &values);
  void addElementData(int dataSetIndex, int element,
                      const std::vector<double> &values);
  void addElementNodeData(int dataSetIndex, int element,
                          const std::vector<std::vector<double>> &values);

  void addInterpolationSchemeSet(std::string name, int topology);
  void addInterpolationMatrix(int dataSetIndex, int nrows, int ncols,
                              const std::vector<double> &values);

  void getNodeIndices(std::vector<int> &indices) const;
  void getElementIndices(std::vector<int> &indices) const;

  void getNode(int index, double &x, double &y, double &z) const;
  void getElement(int index, int &elementType, std::vector<int> &nodes,
                  int &physicalEntity, int &elementaryEntity,
                  std::vector<int> &partitions) const;
  void getElement(int index, int &elementType, std::vector<int> &nodes,
                  int &physicalEntity, int &elementaryEntity) const;
  void getPeriodicEntity(int index, int &dimension, int &slaveEntityTag,
                         int &masterEntityTag) const;
  void getPeriodicNode(int index, int &slaveNode, int &masterNode) const;
  void getPhysicalName(int index, int &dimension, int &number,
                       std::string &name) const;
  void getNodeDataSet(int index, std::vector<std::string> &stringTags,
                      std::vector<double> &realTags,
                      int &numberOfFieldComponents,
                      std::vector<int> &nodeIndices,
                      std::vector<std::vector<double>> &values, int &timeStep,
                      int &partition) const;
  void getNodeDataSet(int index, std::vector<std::string> &stringTags,
                      std::vector<double> &realTags,
                      int &numberOfFieldComponents,
                      std::vector<int> &nodeIndices,
                      std::vector<std::vector<double>> &values) const;
  void getElementDataSet(int index, std::vector<std::string> &stringTags,
                         std::vector<double> &realTags,
                         int &numberOfFieldComponents,
                         std::vector<int> &elementIndices,
                         std::vector<std::vector<double>> &values,
                         int &timeStep, int &partition) const;
  void getElementDataSet(int index, std::vector<std::string> &stringTags,
                         std::vector<double> &realTags,
                         int &numberOfFieldComponents,
                         std::vector<int> &elementIndices,
                         std::vector<std::vector<double>> &values) const;
  void
  getElementNodeDataSet(int index, std::vector<std::string> &stringTags,
                        std::vector<double> &realTags,
                        int &numberOfFieldComponents,
                        std::vector<int> &elementIndices,
                        std::vector<std::vector<std::vector<double>>> &values,
                        int &timeStep, int &partition) const;
  void getElementNodeDataSet(
      int index, std::vector<std::string> &stringTags,
      std::vector<double> &realTags, int &numberOfFieldComponents,
      std::vector<int> &elementIndices,
      std::vector<std::vector<std::vector<double>>> &values) const;
  void getInterpolationSchemeSet(int index, std::string &name, int &topology,
                                 std::vector<int> &nrows,
                                 std::vector<int> &ncols,
                                 std::vector<std::vector<double>> &values);

  void reserveNumberOfNodes(int n);
  void reserveNumberOfElements(int n);

  void resetNodeDataSets();
  void resetElementDataSets();
  void resetElementNodeDataSets();

  void resetDataSets();

  void write(std::ostream &output) const;
  void write(const std::string &fileName) const;
  static GmshData read(std::istream &input, int elementType = 2,
                       int physicalEntity = -1);
  static GmshData read(const std::string &fileName, int elementType = 2,
                       int physicalEntity = -1);

private:
  struct NodeDataSet {

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    int timeStep;
    int numberOfFieldComponents;
    int partition;

    std::vector<int> nodeIndices;
    std::vector<std::vector<double>> values;
  };

  struct ElementDataSet {

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    int timeStep;
    int numberOfFieldComponents;
    int partition;

    std::vector<int> elementIndices;
    std::vector<std::vector<double>> values;
  };

  struct ElementNodeDataSet {

    std::vector<std::string> stringTags;
    std::vector<double> realTags;
    std::vector<int> integerTags;
    int timeStep;
    int numberOfFieldComponents;
    int partition;

    std::vector<int> elementIndices;
    std::vector<std::vector<std::vector<double>>> values;
  };

  struct InterpolationSchemeSet {

    std::vector<int> nrows;
    std::vector<int> ncols;
    std::vector<std::vector<double>> values; // Each outer array element stores
                                             // one matrix in row-major order

    std::string name;
    int topology;
  };

  struct Node {
    double x;
    double y;
    double z;
  };

  struct Element {
    int type;
    int physicalEntity;
    int elementaryEntity;
    std::vector<int> nodes;
    std::vector<int> partitions;
  };

  struct PeriodicEntity {
    int dimension;
    int slaveTag;
    int masterTag;
  };

  struct PeriodicNode {
    int slaveNode;
    int masterNode;
  };

  struct PhysicalName {
    int dimension;
    int number;
    std::string name;
  };

  std::string m_versionNumber;
  int m_fileType;
  int m_dataSize;
  int m_numberOfNodes;
  int m_numberOfElements;

  std::vector<shared_ptr<Node>> m_nodes;
  std::vector<shared_ptr<Element>> m_elements;

  std::set<int> m_elementIndices;

  std::vector<PeriodicEntity> m_periodicEntities;
  std::vector<PeriodicNode> m_periodicNodes;

  std::vector<PhysicalName> m_physicalNames;

  std::vector<shared_ptr<NodeDataSet>> m_nodeDataSets;
  std::vector<shared_ptr<ElementDataSet>> m_elementDataSets;
  std::vector<shared_ptr<ElementNodeDataSet>> m_elementNodeDataSets;

  std::vector<shared_ptr<InterpolationSchemeSet>> m_interpolationSchemeSets;
};

class GmshIo {

public:
  GmshIo(const shared_ptr<const Grid> &grid);
  GmshIo(GmshData gmshData);
  GmshIo(std::string fileName, int physicalEntity = -1);

  shared_ptr<const Grid> grid() const;
  const std::vector<int> &nodePermutation() const;
  const std::vector<int> &elementPermutation() const;
  const std::vector<int> &inverseNodePermutation() const;
  const std::vector<int> &inverseElementPermutation() const;
  const GmshData &gmshData() const;
  void write(std::string fileName) const;
  GmshData &gmshData();

  void resetNodeDataSets();
  void resetElementDataSets();
  void resetElementNodeDataSets();

  void resetDataSets();

private:
  mutable std::vector<int> m_nodePermutation;
  mutable std::vector<int> m_inverseNodePermutation;
  mutable std::vector<int> m_elementPermutation;
  mutable std::vector<int> m_inverseElementPermutation;
  GmshData m_gmshData;
  mutable shared_ptr<const Grid> m_grid;
};

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> gridFunctionFromGmsh(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const GmshIo &gmshIo, shared_ptr<const Grid> grid = shared_ptr<Grid>(),
    GmshPostData::Type gmshPostDataType = GmshPostData::ELEMENT_NODE,
    int index = 0);

template <typename BasisFunctionType, typename ResultType>
void exportToGmsh(GridFunction<BasisFunctionType, ResultType> gridFunction,
                  const char *dataLabel, GmshIo &gmshIo,
                  GmshPostData::Type gmshPostDataType =
                      GmshPostData::ELEMENT_NODE,
                  std::string complexMode = "real");

template <typename BasisFunctionType, typename ResultType>
void exportToGmsh(GridFunction<BasisFunctionType, ResultType> gridFunction,
                  const char *dataLabel, const char *fileName,
                  GmshPostData::Type gmshPostDataType =
                      GmshPostData::ELEMENT_NODE,
                  std::string complexMode = "real");

} // namespace
#endif
