#ifndef __COMMONTYPES_H
#define __COMMONTYPES_H

typedef struct {
  int dim;
  int nvtx;
  int nels;
  int nidx;
} MeshParam;

typedef struct { // see opencl_framework.hpp for details
  int dim;
  int nvtx;
  int nels;
  int nidx;
  __global const ValueType *cl_vtxbuf;
  __global const int *cl_elbuf;
} MeshGeom;

#endif
