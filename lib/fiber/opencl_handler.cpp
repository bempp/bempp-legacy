#include "opencl_handler.hpp"

#ifdef WITH_OPENCL

#include "cl_util.hpp"
#include "CL/commontypes.h.str"

#define CL_DIAGNOSTICS

#ifdef CL_DIAGNOSTICS
#define CALLECHO() { std::cout << "***** CL: called " << __FUNCTION__ << std::endl; }
#else
#define CALLECHO()
#endif

namespace Fiber {

template<typename ValueType, typename IndexType>
OpenClHandler<ValueType,IndexType>::OpenClHandler(const OpenClOptions& options)
{
    cl_int err;

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.size() == 0)
        throw std::runtime_error("Could not obtain OpenCL platforms");

    deviceUsed = 0;
    cl_context_properties properties[] =
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};

    context = cl::Context(CL_DEVICE_TYPE_GPU, properties);
    //devices = context.getInfo<CL_CONTEXT_DEVICES>();
    devices.push_back(context.getInfo<CL_CONTEXT_DEVICES>()[deviceUsed]);
    // note: we collect only a single device, because OpenCL runs into trouble
    // if 'devices' contains multiple entries with different compute capabilities
    // and it tries to generate a program that runs on all of them

    try{
        queue = cl::CommandQueue(context, devices[0/*deviceUsed*/], 0, &err);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%d)\n", er.what(), er.err());
    }
}

template<typename ValueType, typename IndexType>
OpenClHandler<ValueType,IndexType>::~OpenClHandler()
{
}

template<>
const std::string OpenClHandler<float,int>::typedefStr () const
{
    return std::string("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\ntypedef float ValueType;\n");
}

template<>
const std::string OpenClHandler<double,int>::typedefStr () const
{
    return std::string("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\ntypedef double ValueType;\n");
}

template<typename ValueType, typename IndexType>
const std::string OpenClHandler<ValueType,IndexType>::initStr () const
{
    CALLECHO();

    std::string str(typedefStr());
    return str.append (commontypes_h, commontypes_h_len);
}

template<typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::loadProgramFromString(
    std::string strSource)
{
    CALLECHO();

    int pl = strSource.size();
    try
    {
        cl::Program::Sources source(1,
            std::make_pair(strSource.c_str(), pl));
        program = cl::Program(context, source);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }

    try
    {
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

template<typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::loadProgramFromStringArray (
    std::vector<std::string> strSources) const
{
    CALLECHO();

    try {
        cl::Program::Sources source;
	for (int i = 0; i < strSources.size(); i++) {
 	    source.push_back (std::make_pair (strSources[i].c_str(),
					      strSources[i].size()));
	}
	program = cl::Program(context, source);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }

    try
    {
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

template <typename ValueType, typename IndexType>
cl::Kernel &OpenClHandler<ValueType,IndexType>::setKernel (const char *kernelName) const
{
    CALLECHO();

    cl_int err;
    try{
        kernel = cl::Kernel(program, kernelName, &err);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%d)\n", er.what(), er.err());
    }
    return kernel;
}

template <typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::enqueueKernel (const cl::NDRange &global) const
{
    CALLECHO();

    cl_int err;
    err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, cl::NullRange,
				     NULL, &event);
}

template <typename ValueType, typename IndexType>
int OpenClHandler<ValueType,IndexType>::SetGeometryArgs (cl::Kernel &kernel, int argid) const
{
    CALLECHO();

    kernel.setArg (argid++, meshgeom.size);
    kernel.setArg (argid++, meshgeom.cl_vtxbuf);
    kernel.setArg (argid++, meshgeom.cl_elbuf);
    return argid;
}

template <typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::pushGeometry (
    const arma::Mat<ValueType>& vtx,
    const arma::Mat<IndexType>& idx) const
{
    CALLECHO();
    cl_int err;

    // Allocate the buffers
    meshgeom.size.dim = vtx.n_rows;
    meshgeom.size.nvtx = vtx.n_cols;
    size_t vtxbuf_size = meshgeom.size.dim*meshgeom.size.nvtx*sizeof(ValueType);
    meshgeom.cl_vtxbuf = cl::Buffer (context, CL_MEM_READ_ONLY,
				     vtxbuf_size, NULL, &err);

    meshgeom.size.nels = idx.n_cols;
    meshgeom.size.nidx = idx.n_rows;
    size_t idxbuf_size = meshgeom.size.nels*meshgeom.size.nidx*sizeof(IndexType);
    meshgeom.cl_elbuf = cl::Buffer (context, CL_MEM_READ_ONLY, idxbuf_size,
				    NULL, &err);
			     
    // Copy the geometry data to device memory
    try
    {
	  queue.enqueueWriteBuffer (meshgeom.cl_vtxbuf, CL_TRUE, 0,
				    vtxbuf_size, vtx.memptr(), NULL, &event);
	  queue.enqueueWriteBuffer (meshgeom.cl_elbuf, CL_TRUE, 0,
				    idxbuf_size, idx.memptr(), NULL, &event);
    }
    catch (cl::Error er) {
        printf ("Push mesh geometry: %s\n", oclErrorString(er.err()));
    }
    queue.finish();
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::createValueBuffer (int size,
    cl_mem_flags usage) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = size * sizeof(ValueType);
    return new cl::Buffer (context, usage, bufsize, NULL, &err);
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::pushValueVector (
    const std::vector<ValueType>& vec) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = vec.size() * sizeof(ValueType);
    const ValueType *vecptr = &vec[0];
    cl::Buffer *clbuf = new cl::Buffer (context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
    queue.enqueueWriteBuffer (*clbuf, CL_TRUE, 0, bufsize, vecptr, NULL, &event);
    return clbuf;
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::pushIndexVector (
    const std::vector<IndexType>& vec) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = vec.size() * sizeof(IndexType);
    const IndexType *vecptr = &vec[0];
    cl::Buffer *clbuf = new cl::Buffer (context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
    queue.enqueueWriteBuffer (*clbuf, CL_TRUE, 0, bufsize, vecptr, NULL, &event);
    return clbuf;
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::pushValueRow (
    const arma::Row<ValueType>& row) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = row.n_rows * row.n_cols * sizeof(ValueType);
    cl::Buffer *clbuf = new cl::Buffer (context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
    queue.enqueueWriteBuffer (*clbuf, CL_TRUE, 0, bufsize, row.memptr(), NULL, &event);
    return clbuf;
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::pushValueMatrix (
    const arma::Mat<ValueType>& mat) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = mat.n_rows * mat.n_cols * sizeof(ValueType);
    cl::Buffer *clbuf = new cl::Buffer (context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
    queue.enqueueWriteBuffer (*clbuf, CL_TRUE, 0, bufsize, mat.memptr(), NULL, &event);
    return clbuf;
}

template <typename ValueType, typename IndexType>
cl::Buffer *OpenClHandler<ValueType,IndexType>::pushValueCube (
    const arma::Cube<ValueType>& cube) const
{
    CALLECHO();

    cl_int err;
    size_t bufsize = cube.n_rows * cube.n_cols * cube.n_slices * sizeof(ValueType);
    cl::Buffer *clbuf = new cl::Buffer (context, CL_MEM_READ_ONLY, bufsize, NULL, &err);
    queue.enqueueWriteBuffer (*clbuf, CL_TRUE, 0, bufsize, cube.memptr(), NULL, &event);
    return clbuf;
}

template <typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::pullValueVector (
    const cl::Buffer &clbuf,
    std::vector<ValueType>& vec,
    int size) const
{
    CALLECHO();

    ValueType *vecptr = &vec[0];
    size_t bufsize = size*sizeof(ValueType);
    queue.enqueueReadBuffer (clbuf, CL_TRUE, 0, bufsize, vecptr, NULL, &event);
}

template <typename ValueType, typename IndexType>
void OpenClHandler<ValueType,IndexType>::pullValueCube (
    const cl::Buffer &clbuf,
    arma::Cube<ValueType>& cube) const
{
    CALLECHO();

    size_t bufsize = cube.n_rows * cube.n_cols * cube.n_slices * sizeof(ValueType);
    queue.enqueueReadBuffer (clbuf, CL_TRUE, 0, bufsize, cube.memptr(), NULL, &event);
}


#ifdef COMPILE_FOR_FLOAT
template class OpenClHandler<float, int>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class OpenClHandler<double, int>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class OpenClHandler<std::complex<float>, int>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class OpenClHandler<std::complex<double>, int>;
#endif

} // namespace Fiber

#endif // WITH_OPENCL
