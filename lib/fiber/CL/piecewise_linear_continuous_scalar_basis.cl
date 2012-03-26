// -*-C++-*-

/**
 * \file basisf_lintri.cl
 * OpenCL implementation for shape function computation on a linear
 * triangle type.
 */

// Compute basis functions for a single element and a specified range of dofs
__kernel void clBasisfElement (
        MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	int elIdx,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	int nDof,
	int dofIdx,   // only used if nDof==1
	__global ValueType *result
)
{
    int pt = get_global_id(0);
    if (pt >= nPoints) return;

    if (nDof==1) {
        switch (dofIdx) {
	case 0:
	    result[pt] = 1.0-g_localPoints[pt*pointDim]-g_localPoints[pt*pointDim+1];
	    break;
	case 1:
	    result[pt] = g_localPoints[pt*pointDim];
	    break;
	case 2:
	    result[pt] = g_localPoints[pt*pointDim+1];
	    break;
	}
    } else {
        ValueType px = g_localPoints[pt*pointDim];
	ValueType py = g_localPoints[pt*pointDim+1];
	result[pt*nDof] = 1.0-px-py;
	result[pt*nDof+1] = px;
	result[pt*nDof+2] = py;
    }
    // note 1: basis functions are assumed independent of global geometry (for now)
    // note 3: component count currently assumed 1
}

// Compute basis functions for multiple elements
__kernel void clBasisfElements (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const int *g_elIdx,
	int nEls,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	int nDof,
	__global ValueType *result
)
{
    int el = get_global_id(0);
    if (el >= nEls) return;

    int pt = get_global_id(1);
    if (pt >= nPoints) return;

    int ofs = el*nPoints*nDof + pt*nDof;
    ValueType px = g_localPoints[pt*pointDim];
    ValueType py = g_localPoints[pt*pointDim+1];
    result[ofs++] = 1.0-px-py;
    result[ofs++] = px;
    result[ofs++] = py;
    // note 1: basis functions are assumed independent of global geometry (for now)
    // note 2: this computes basis functions for all dofs
    // note 3: component count currently assumed 1
}


#ifdef UNDEF
/**
 * \brief Evaluate linear basis function on a reference triangle
 * \param [in] point barycentric coordinates of point
 * \param [out] data basis function values for all DOFs
 * \note Definition of reference triangle: v = {{0,0,0},{1,0,0},{0,1,0}}
 * \note z-coordinates of all input points must be 0
 */
void devBasisEval (const ValueType *point,
		   ValueType *data)
{
    data[0] = 1 - point[0] - point[1];
    data[1] = point[0];
    data[2] = point[1];
}

/**
 * \brief Kernel function: evaluate products of basis functions on
 *   reference triangle
 * \param [in] g_p1 Array of reference points on first triangle [np1]
 * \param [in] g_p2 Array of reference points on second triangle [np2]
 * \param [in] np Number of points in g_p1 and g_p2
 * \param [out] bfprod Array of basis function products [np*3*3]
 * \note Elements in bfprod are ordered such that
 *   bfprod[i*9 + n*3+m]
 *   refers to point i and basis function m in first triangle, point i and
 *   basis function n in second triangle,
 */
__kernel void clBasisProd (__global const ValueType *g_p1,
			   __global const ValueType *g_p2,
			   int np,
			   __global ValueType *g_bfprod)
{
    int i = get_global_id(0);

    // range check
    if (i >= np) return;
    
    ValueType p1[3];
    ValueType p2[3];
    ValueType bf1[3];
    ValueType bf2[3];
    for (int j = 0; j < 3; j++) {
        p1[j] = g_p1[i*3+j];
	p2[j] = g_p2[i*3+j];
    }

    devBasisEval (p1, bf1);
    devBasisEval (p2, bf2);

    int m, n, ofs = i * 9;
    for (n = 0; n < 3; n++)
        for (m = 0; m < 3; m++)
	    g_bfprod[ofs + n*3 + m] = bf1[m] * bf2[n];
}
#endif
