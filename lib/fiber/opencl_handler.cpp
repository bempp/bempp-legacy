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

#include "bempp/common/config_data_types.hpp"
#include "bempp/common/config_opencl.hpp"
#include "opencl_handler.hpp"

#ifdef WITH_OPENCL

#include "cl_util.hpp"
#include "CL/commontypes.h.str"

//#define CL_DIAGNOSTICS

#ifdef CL_DIAGNOSTICS
#define CALLECHO()                                                             \
  { std::cout << "***** CL: called " << __FUNCTION__ << std::endl; }
#else
#define CALLECHO()
#endif

namespace Fiber {

OpenClHandler::OpenClHandler(const OpenClOptions &options) {
  cl_int err;
  useOpenCl = options.useOpenCl;
  nProgBuf = 0;

  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  if (platforms.size() == 0)
    throw std::runtime_error("Could not obtain OpenCL platforms");

  deviceUsed = 0;
  cl_context_properties properties[] = {
      CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};

  context = cl::Context(CL_DEVICE_TYPE_GPU, properties);
  // devices = context.getInfo<CL_CONTEXT_DEVICES>();
  devices.push_back(context.getInfo<CL_CONTEXT_DEVICES>()[deviceUsed]);
  // note: we collect only a single device, because OpenCL runs into trouble
  // if 'devices' contains multiple entries with different compute capabilities
  // and it tries to generate a program that runs on all of them

  try {
    queue = cl::CommandQueue(context, devices[0 /*deviceUsed*/], 0, &err);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%d)\n", er.what(), er.err());
  }
}

OpenClHandler::~OpenClHandler() {
  if (nProgBuf)
    delete[] progBuf;
}

// const std::pair<const char*,int> OpenClHandler::typedefStr () const
//{
//    static const char *str =
//        "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\ntypedef float
// ValueType;\n";
//    static int len = strlen(str);
//    return std::make_pair(str,len);
//}

const std::pair<const char *, int> OpenClHandler::typedefStr() const {
  static const char *str = "#pragma OPENCL EXTENSION cl_khr_fp64 : "
                           "enable\ntypedef double ValueType;\n";
  static int len = strlen(str);
  return std::make_pair(str, len);
}

const std::pair<const char *, int> OpenClHandler::initStr() const {
  CALLECHO();

  static std::pair<char *, int> str;
  static bool need_setup = true;

  if (need_setup) {
    std::pair<const char *, int> tdef = typedefStr();
    str.first = new char[tdef.second + commontypes_h_len + 1];
    strcpy(str.first, tdef.first);
    strncpy(str.first + tdef.second, commontypes_h, commontypes_h_len);
    str.first[tdef.second + commontypes_h_len] = '\0';
    str.second = tdef.second + commontypes_h_len;
    need_setup = false;
  }
  return str;
}

void OpenClHandler::loadProgramFromString(std::string strSource) {
  CALLECHO();

  int pl = strSource.size();
  try {
    cl::Program::Sources source(1, std::make_pair(strSource.c_str(), pl));
    program = cl::Program(context, source);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }

  try {
    program.build(devices);
  }
  catch (cl::Error er) {
    printf("program.build: %s\n", oclErrorString(er.err()));
  }

#ifdef CL_DIAGNOSTICS
  std::cout << "OpenCL: Build program" << std::endl;
  std::cout << "Build Status: "
            << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0])
            << std::endl;
  std::cout << "Build Options:\t"
            << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0])
            << std::endl;
  std::cout << "Build Log:\t "
            << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
            << std::endl;
#endif
  std::cout << "x11" << std::endl;
}

void OpenClHandler::loadProgramFromStringArray(cl::Program::Sources strSources)
    const {
  CALLECHO();

  int i, nStr = strSources.size();

  // Check whether program is already loaded
  if (nStr == nProgBuf) {
    for (i = 0; i < nProgBuf; i++)
      if (progBuf[i] != strSources[i].first)
        break;
    if (i == nProgBuf) {
      return; // already loaded
    }
  } else {
    if (nProgBuf < nStr) {
      if (nProgBuf)
        delete[] progBuf;
      progBuf = new const char *[nStr];
    }
  }
  for (i = 0; i < nStr; i++)
    progBuf[i] = strSources[i].first;
  nProgBuf = nStr;

  try {
    program = cl::Program(context, strSources);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }

  try {
    program.build(devices);
  }
  catch (cl::Error er) {
    printf("program.build: %s\n", oclErrorString(er.err()));
  }

#ifdef CL_DIAGNOSTICS
  std::cout << "OpenCL: Build program" << std::endl;
  std::cout << "Build Status: "
            << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0])
            << std::endl;
  std::cout << "Build Options:\t"
            << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0])
            << std::endl;
  std::cout << "Build Log:\t "
            << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
            << std::endl;
#endif
}

cl::Kernel &OpenClHandler::setKernel(const char *kernelName) const {
  CALLECHO();

  cl_int err;
  try {
    kernel = cl::Kernel(program, kernelName, &err);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%d)\n", er.what(), er.err());
  }
  return kernel;
}

void OpenClHandler::enqueueKernel(const cl::NDRange &global) const {
  CALLECHO();

  cl_int err;
  err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, cl::NullRange,
                                   NULL, &event);
}

int OpenClHandler::SetGeometryArgs(cl::Kernel &kernel, int argid) const {
  CALLECHO();

  kernel.setArg(argid++, meshgeom.size);
  kernel.setArg(argid++, meshgeom.cl_vtxbuf);
  kernel.setArg(argid++, meshgeom.cl_elbuf);
  return argid;
}

template <typename CoordinateType, typename IndexType>
void OpenClHandler::pushGeometry(const arma::Mat<CoordinateType> &vtx,
                                 const arma::Mat<IndexType> &idx) const {
  CALLECHO();
  cl_int err;

  // Allocate the buffers
  meshgeom.size.dim = vtx.n_rows;
  meshgeom.size.nvtx = vtx.n_cols;
  size_t vtxbuf_size =
      meshgeom.size.dim * meshgeom.size.nvtx * sizeof(CoordinateType);
  meshgeom.cl_vtxbuf =
      cl::Buffer(context, CL_MEM_READ_ONLY, vtxbuf_size, NULL, &err);

  meshgeom.size.nels = idx.n_cols;
  meshgeom.size.nidx = idx.n_rows;
  size_t idxbuf_size =
      meshgeom.size.nels * meshgeom.size.nidx * sizeof(IndexType);
  meshgeom.cl_elbuf =
      cl::Buffer(context, CL_MEM_READ_ONLY, idxbuf_size, NULL, &err);

  // Copy the geometry data to device memory
  try {
    queue.enqueueWriteBuffer(meshgeom.cl_vtxbuf, CL_TRUE, 0, vtxbuf_size,
                             vtx.memptr(), NULL, &event);
    queue.enqueueWriteBuffer(meshgeom.cl_elbuf, CL_TRUE, 0, idxbuf_size,
                             idx.memptr(), NULL, &event);
  }
  catch (cl::Error er) {
    printf("Push mesh geometry: %s\n", oclErrorString(er.err()));
  }
  queue.finish();
}

template <typename BufferType>
cl::Buffer *OpenClHandler::createBuffer(int size, cl_mem_flags usage) const {
  CALLECHO();

  cl_int err;
  size_t bufsize = size * sizeof(BufferType);
  return new cl::Buffer(context, usage, bufsize, NULL, &err);
}

template <typename BufferType>
cl::Buffer *
OpenClHandler::pushVector(const std::vector<BufferType> &vec) const {
  CALLECHO();

  cl_int err;
  size_t bufsize = vec.size() * sizeof(BufferType);
  const BufferType *vecptr = &vec[0];
  cl::Buffer *clbuf =
      new cl::Buffer(context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
  queue.enqueueWriteBuffer(*clbuf, CL_TRUE, 0, bufsize, vecptr, NULL, &event);
  return clbuf;
}

template <typename BufferType>
cl::Buffer *OpenClHandler::pushBuffer(const BufferType *buf, int size) const {
  CALLECHO();

  cl_int err;
  size_t bufsize = size * sizeof(BufferType);
  cl::Buffer *clbuf =
      new cl::Buffer(context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
  queue.enqueueWriteBuffer(*clbuf, CL_TRUE, 0, bufsize, buf, NULL, &event);
  return clbuf;
}

template <typename BufferType>
cl::Buffer *OpenClHandler::pushRow(const arma::Row<BufferType> &row) const {
  CALLECHO();

  cl_int err;
  size_t bufsize = row.n_rows * row.n_cols * sizeof(BufferType);
  cl::Buffer *clbuf =
      new cl::Buffer(context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
  queue.enqueueWriteBuffer(*clbuf, CL_TRUE, 0, bufsize, row.memptr(), NULL,
                           &event);
  return clbuf;
}

template <typename BufferType>
cl::Buffer *OpenClHandler::pushMatrix(const arma::Mat<BufferType> &mat) const {
  CALLECHO();

  cl_int err;
  size_t bufsize = mat.n_rows * mat.n_cols * sizeof(BufferType);
  cl::Buffer *clbuf =
      new cl::Buffer(context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
  queue.enqueueWriteBuffer(*clbuf, CL_TRUE, 0, bufsize, mat.memptr(), NULL,
                           &event);
  return clbuf;
}

template <typename BufferType>
cl::Buffer *OpenClHandler::pushCube(const arma::Cube<BufferType> &cube) const {
  CALLECHO();

  cl_int err;
  size_t bufsize =
      cube.n_rows * cube.n_cols * cube.n_slices * sizeof(BufferType);
  cl::Buffer *clbuf =
      new cl::Buffer(context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
  queue.enqueueWriteBuffer(*clbuf, CL_TRUE, 0, bufsize, cube.memptr(), NULL,
                           &event);
  return clbuf;
}

template <typename BufferType>
void OpenClHandler::pullVector(const cl::Buffer &clbuf,
                               std::vector<BufferType> &vec, int size) const {
  CALLECHO();

  BufferType *vecptr = &vec[0];
  size_t bufsize = size * sizeof(BufferType);
  queue.enqueueReadBuffer(clbuf, CL_TRUE, 0, bufsize, vecptr, NULL, &event);
}

template <typename BufferType>
void OpenClHandler::pullCube(const cl::Buffer &clbuf,
                             arma::Cube<BufferType> &cube) const {
  CALLECHO();

  size_t bufsize =
      cube.n_rows * cube.n_cols * cube.n_slices * sizeof(BufferType);
  queue.enqueueReadBuffer(clbuf, CL_TRUE, 0, bufsize, cube.memptr(), NULL,
                          &event);
}

// Integer instantiations
template cl::Buffer *
OpenClHandler::pushVector<int>(const std::vector<int> &vec) const;

// Single-precision instantiations
template void
OpenClHandler::pushGeometry<float, int>(const arma::Mat<float> &vtx,
                                        const arma::Mat<int> &idx) const;
template cl::Buffer *
OpenClHandler::createBuffer<float>(int size, cl_mem_flags usage) const;
template cl::Buffer *
OpenClHandler::pushVector<float>(const std::vector<float> &vec) const;
template cl::Buffer *
OpenClHandler::pushMatrix<float>(const arma::Mat<float> &mat) const;
template void OpenClHandler::pullVector<float>(const cl::Buffer &clbuf,
                                               std::vector<float> &vec,
                                               int size) const;
template void OpenClHandler::pullCube<float>(const cl::Buffer &clbuf,
                                             arma::Cube<float> &cube) const;

// Double-precision instantiations
template void
OpenClHandler::pushGeometry<double, int>(const arma::Mat<double> &vtx,
                                         const arma::Mat<int> &idx) const;
template cl::Buffer *
OpenClHandler::createBuffer<double>(int size, cl_mem_flags usage) const;
template cl::Buffer *
OpenClHandler::pushVector<double>(const std::vector<double> &vec) const;
template cl::Buffer *
OpenClHandler::pushMatrix<double>(const arma::Mat<double> &mat) const;
template void OpenClHandler::pullVector<double>(const cl::Buffer &clbuf,
                                                std::vector<double> &vec,
                                                int size) const;
template void OpenClHandler::pullCube<double>(const cl::Buffer &clbuf,
                                              arma::Cube<double> &cube) const;

// Single-precision complex instantiations
template cl::Buffer *
OpenClHandler::createBuffer<std::complex<float>>(int size,
                                                 cl_mem_flags usage) const;
template cl::Buffer *OpenClHandler::pushVector<std::complex<float>>(
    const std::vector<std::complex<float>> &vec) const;
template void OpenClHandler::pullVector<std::complex<float>>(
    const cl::Buffer &clbuf, std::vector<std::complex<float>> &vec,
    int size) const;
template void OpenClHandler::pullCube<std::complex<float>>(
    const cl::Buffer &clbuf, arma::Cube<std::complex<float>> &cube) const;

// Double-precision complex instantiations
template cl::Buffer *
OpenClHandler::createBuffer<std::complex<double>>(int size,
                                                  cl_mem_flags usage) const;
template cl::Buffer *OpenClHandler::pushVector<std::complex<double>>(
    const std::vector<std::complex<double>> &vec) const;
template void OpenClHandler::pullVector<std::complex<double>>(
    const cl::Buffer &clbuf, std::vector<std::complex<double>> &vec,
    int size) const;
template void OpenClHandler::pullCube<std::complex<double>>(
    const cl::Buffer &clbuf, arma::Cube<std::complex<double>> &cube) const;

} // namespace Fiber

#endif // WITH_OPENCL
