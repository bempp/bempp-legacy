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

#include "../common/common.hpp"

#include "bempp/common/config_opencl.hpp"

#ifndef fiber_opencl_handler_hpp
#define fiber_opencl_handler_hpp

#define __CL_ENABLE_EXCEPTIONS

#include <string>
#include <vector>
#include "../common/armadillo_fwd.hpp"
#include "opencl_options.hpp"

#ifdef WITH_OPENCL
#include "CL/cl.hpp"
#endif

// TODO: rewrite the constructor of OpenClHandler.
// It should take a bool useOpenCl and *in addition to that* openClOptions.
// The role of the latter should be to e.g. select the device to use
// and other configurable execution parameters.
// If there are no such parameters, OpenClOptions should just be removed.

// Justification: right now there can be a conflict: the user can invoke
// AssemblyOptions::switchToOpenCl() and pass to it an instance of OpenClOptions
// with useOpenCl set to false. This makes no sense.

namespace Fiber {

#ifdef WITH_OPENCL

class OpenClHandler {
public:
  /**
   * \brief Default constructor.
   * Initializes OpenCL context and automatically chooses platform and device
   */
  OpenClHandler(const OpenClOptions &options);

  /**
   * Default destructor.
   * Releases OpenCL objects and frees device memory
   */
  ~OpenClHandler();

  bool UseOpenCl() const { return useOpenCl; }

  /**
   * \brief Returns a string that defines the 'ValueType' type used by most
   * device
   *   functions, plus some commonly used definitions
   */
  const std::pair<const char *, int> initStr() const;

  /**
   * \brief Load an OpenCL program from a string.
   * Pass in the kernel source code as a string.
   * \param strSource OpenCL kernel source
   */
  void loadProgramFromString(std::string strSource);

  /**
   * \brief Load an OpenCL program as a concatenation of multiple strings
   * \param sources Vector of source strings
   */
  void loadProgramFromStringArray(cl::Program::Sources sources) const;

  cl::Kernel &setKernel(const char *kernelName) const;

  /**
   * \brief Enqueue the currently set kernel for execution
   */
  void enqueueKernel(const cl::NDRange &global) const;

  /**
   * \brief Push the mesh geometry to device memory
   * \param mesh mesh instance
   * \note This method populates the cl_meshvtx and cl_meshidx buffers
   */
  template <typename CoordinateType, typename IndexType>
  void pushGeometry(const arma::Mat<CoordinateType> &vtx,
                    const arma::Mat<IndexType> &idx) const;

  /**
   * \brief Pass mesh geometry arguments to a kernel.
   * \param kernel OpenCL kernel to receive the arguments
   * \param argid Starting argument index for the parameters
   * \return First index after the geometry parameters
   */
  int SetGeometryArgs(cl::Kernel &kernel, int argid) const;

  /**
   * \brief Creates an OpenCL buffer for holding 'size' ValueType elements
   * \param size number of elements
   * \param usage CL memory flags, including CL_MEM_READ_ONLY, CL_MEM_WRITE_ONLY
   */
  template <typename BufferType>
  cl::Buffer *createBuffer(int size, cl_mem_flags usage) const;

  /**
   * \brief Push an index vector into an OpenCL buffer and return a pointer to
   * the buffer
   */
  template <typename BufferType>
  cl::Buffer *pushVector(const std::vector<BufferType> &vec) const;

  /**
   * \brief Push a vector of index pairs into an OpenCL buffer and return a
   * pointer to the buffer
   */
  template <typename BufferType>
  cl::Buffer *pushBuffer(const BufferType *buf, int size) const;

  /**
   * \brief Push a value vector into an OpenCL buffer and return a pointer to
   * the buffer
   */
  // cl::Buffer *pushValueVector (const std::vector<ValueType> &vec) const;

  /**
   * \brief Push a row into an OpenCL buffer and return a pointer to the buffer
   */
  template <typename BufferType>
  cl::Buffer *pushRow(const arma::Row<BufferType> &row) const;

  /**
   * \brief Push an index array to the OpenCL buffer and return the buffer
   */
  // cl::Buffer *pushIndexList (const arma::Row<IndexType>& idx);

  /**
   * \brief Push a matrix into an OpenCL buffer and return a pointer to the
   * buffer
   */
  template <typename BufferType>
  cl::Buffer *pushMatrix(const arma::Mat<BufferType> &mat) const;

  /**
   * \brief Push a matrix into an OpenCL buffer and return a pointer to the
   * buffer
   */
  template <typename BufferType>
  cl::Buffer *pushCube(const arma::Cube<BufferType> &cube) const;

  template <typename BufferType>
  void pullVector(const cl::Buffer &clbuf, std::vector<BufferType> &vec,
                  int size) const;

  template <typename BufferType>
  void pullCube(const cl::Buffer &clbuf, arma::Cube<BufferType> &cube) const;

  struct MeshGeom {
    struct MeshDims {
      int dim;  ///< Mesh dimension
      int nvtx; ///< Number of vertices
      int nels; ///< Number of elements
      int nidx; ///< Max number of indices per element
    } size;
    cl::Buffer cl_vtxbuf; ///< Mesh geometry: vertex list
    cl::Buffer cl_elbuf;  ///< Mesh geometry: node index list
  };

  const MeshGeom &meshGeom() const { return meshgeom; }

private:
  const std::pair<const char *, int> typedefStr() const;

  bool useOpenCl;
  unsigned int deviceUsed;
  std::vector<cl::Device> devices;
  cl::CommandQueue queue;
  cl::Context context;
  mutable cl::Program program;
  mutable cl::Kernel kernel;
  mutable cl::Event event;

  mutable MeshGeom meshgeom;

  mutable const char **progBuf;
  mutable int nProgBuf;
};

#else

// Dummy implementation for the OpenCL handler

class OpenClHandler {
public:
  OpenClHandler(const OpenClOptions &options) {}

  bool UseOpenCl() const { return false; }

  template <typename CoordinateType, typename IndexType>
  void pushGeometry(const arma::Mat<CoordinateType> &vtx,
                    const arma::Mat<IndexType> &idx) const {}
};

#endif // WITH_OPENCL

} // namespace Fiber

#endif // fiber_opencl_handler_hpp
