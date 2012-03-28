// -*-C++-*-

/**
 * \file separable_numerical_double_integrator.cl
 * CL code for integrating DOFs in a row or column
 */

// This should be moved to a general mesh geometry utilities file
void devGetGeometry (
	__global const ValueType *g_meshVtx,
	int nMeshVtx,
	__global const int *g_meshIdx,
	int meshDim,
	int nElVertices,
	int elIdx,
	ValueType *vtx)
{
    int i, j, k, vtxIdx;
    for (i = k = 0; i < nElVertices; i++) {
        vtxIdx = g_meshIdx[elIdx*nElVertices + i];
        for (j = 0; j < meshDim; j++)
  	    vtx[k++] = g_meshVtx[vtxIdx*meshDim + j];
    }
}

//// This should be moved to the appropriate mapping class
//void devMapPoint (const ValueType *elgeom, const __global ValueType *g_localPoint, ValueType *globalPoint)
//{
//}


// Map a list of reference points to a single global element
// result (from slowest to fastest changing index): nPoints*meshDim
__kernel void clMapPointsToElement (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	int elIdx,
	__global ValueType *mappedPoints,
	__global ValueType *integrationElements
)
{
    int i, j;
    int pt = get_global_id(0);
    if (pt >= nPoints) return;

    int vtxIdx[4];
    for (i = 0; i < prm.nidx; i++)
        vtxIdx[i] = g_meshIdx[elIdx*prm.nidx+i]; // vertices associated with the element

    ValueType fun[3];
    ValueType px = g_localPoints[pt*pointDim];
    ValueType py = g_localPoints[pt*pointDim+1];
    fun[0] = 1.0-px-py;
    fun[1] = px;
    fun[2] = py;

    // map points
    int ofsPoint = pt*prm.dim;
    for (i = 0; i < prm.dim; i++) {
        ValueType res = 0.0;
        for (j = 0; j < prm.dim; j++) {
	    res += fun[j] * g_meshVtx[vtxIdx[j]*prm.dim+i];
	}
	mappedPoints[ofsPoint++] = res;
    }

    // normals and determinant of Jacobian. Note: for now, assumes 3 vertices and 3 coordinates
    ValueType ux = g_meshVtx[vtxIdx[1]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType uy = g_meshVtx[vtxIdx[1]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType uz = g_meshVtx[vtxIdx[1]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType vx = g_meshVtx[vtxIdx[2]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType vy = g_meshVtx[vtxIdx[2]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType vz = g_meshVtx[vtxIdx[2]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType nx = uy*vz - uz*vy;
    ValueType ny = uz*vx - ux*vz;
    ValueType nz = ux*vy - uy*vx;
    ValueType size = sqrt (nx*nx + ny*ny + nz*nz);
    integrationElements[pt] = size;
}


// Map a list of reference points to one or more global elements
// result (from slowest to fastest changing index): nEls*nPoints*meshDim
__kernel void clMapPointsToElements (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	__global const int *g_elIdx,
	int nEls,
	__global ValueType *mappedPoints,
	__global ValueType *integrationElements
)
{
    int i, j;
    int el = get_global_id(0);
    if (el >= nEls) return;

    int pt = get_global_id(1);
    if (pt >= nPoints) return;

    int elIdx = g_elIdx[el];  // the element to operate on

    int vtxIdx[4];
    for (i = 0; i < prm.nidx; i++)
        vtxIdx[i] = g_meshIdx[elIdx*prm.nidx+i]; // vertices associated with the element

    ValueType fun[3];
    ValueType px = g_localPoints[pt*pointDim];
    ValueType py = g_localPoints[pt*pointDim+1];
    fun[0] = 1.0-px-py;
    fun[1] = px;
    fun[2] = py;

    int ofsPoint = (el*nPoints + pt)*prm.dim;
    for (i = 0; i < prm.dim; i++) {
        ValueType res = 0.0;
        for (j = 0; j < prm.dim; j++)
	    res += fun[j] * g_meshVtx[vtxIdx[j]*prm.dim+i];
	mappedPoints[ofsPoint++] = res;
    }

    // determinant of Jacobian. Note: for now, assumes 3 vertices and 3 coordinates
    ValueType ux = g_meshVtx[vtxIdx[1]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType uy = g_meshVtx[vtxIdx[1]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType uz = g_meshVtx[vtxIdx[1]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType vx = g_meshVtx[vtxIdx[2]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType vy = g_meshVtx[vtxIdx[2]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType vz = g_meshVtx[vtxIdx[2]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType nx = uy*vz - uz*vy;
    ValueType ny = uz*vx - ux*vz;
    ValueType nz = ux*vy - uy*vx;
    ValueType size = sqrt (nx*nx + ny*ny + nz*nz);
    integrationElements[el*nPoints+pt] = size;
}


// Map a list of reference points to a single global element
// result (from slowest to fastest changing index): nPoints*meshDim
__kernel void clMapPointsAndNormalsToElement (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	int elIdx,
	__global ValueType *mappedPoints,
	__global ValueType *mappedNormals,
	__global ValueType *integrationElements
)
{
    int i, j;
    int pt = get_global_id(0);
    if (pt >= nPoints) return;

    int vtxIdx[4];
    for (i = 0; i < prm.nidx; i++)
        vtxIdx[i] = g_meshIdx[elIdx*prm.nidx+i]; // vertices associated with the element

    ValueType fun[3];
    ValueType px = g_localPoints[pt*pointDim];
    ValueType py = g_localPoints[pt*pointDim+1];
    fun[0] = 1.0-px-py;
    fun[1] = px;
    fun[2] = py;

    // map points
    int ofsPoint = pt*prm.dim;
    int ofsNormal = ofsPoint;
    for (i = 0; i < prm.dim; i++) {
        ValueType res = 0.0;
        for (j = 0; j < prm.dim; j++) {
	    res += fun[j] * g_meshVtx[vtxIdx[j]*prm.dim+i];
	}
	mappedPoints[ofsPoint++] = res;
    }

    // normals and determinant of Jacobian. Note: for now, assumes 3 vertices and 3 coordinates
    ValueType ux = g_meshVtx[vtxIdx[1]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType uy = g_meshVtx[vtxIdx[1]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType uz = g_meshVtx[vtxIdx[1]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType vx = g_meshVtx[vtxIdx[2]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType vy = g_meshVtx[vtxIdx[2]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType vz = g_meshVtx[vtxIdx[2]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType nx = uy*vz - uz*vy;
    ValueType ny = uz*vx - ux*vz;
    ValueType nz = ux*vy - uy*vx;
    ValueType size = sqrt (nx*nx + ny*ny + nz*nz);
    mappedNormals[ofsNormal++] = nx/size;
    mappedNormals[ofsNormal++] = ny/size;
    mappedNormals[ofsNormal++] = nz/size;
    integrationElements[pt] = size;
}


// Map a list of reference points to one or more global elements
// result (from slowest to fastest changing index): nEls*nPoints*meshDim
__kernel void clMapPointsAndNormalsToElements (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const ValueType *g_localPoints,
	int nPoints,
	int pointDim,
	__global const int *g_elIdx,
	int nEls,
	__global ValueType *mappedPoints,
	__global ValueType *mappedNormals,
	__global ValueType *integrationElements
)
{
    int i, j;
    int el = get_global_id(0);
    if (el >= nEls) return;

    int pt = get_global_id(1);
    if (pt >= nPoints) return;

    int elIdx = g_elIdx[el];  // the element to operate on

    int vtxIdx[4];
    for (i = 0; i < prm.nidx; i++)
        vtxIdx[i] = g_meshIdx[elIdx*prm.nidx+i]; // vertices associated with the element

    ValueType fun[3];
    ValueType px = g_localPoints[pt*pointDim];
    ValueType py = g_localPoints[pt*pointDim+1];
    fun[0] = 1.0-px-py;
    fun[1] = px;
    fun[2] = py;

    int ofsPoint = (el*nPoints + pt)*prm.dim;
    int ofsNormal = ofsPoint;
    for (i = 0; i < prm.dim; i++) {
        ValueType res = 0.0;
        for (j = 0; j < prm.dim; j++)
	    res += fun[j] * g_meshVtx[vtxIdx[j]*prm.dim+i];
	mappedPoints[ofsPoint++] = res;
    }

    // determinant of Jacobian. Note: for now, assumes 3 vertices and 3 coordinates
    ValueType ux = g_meshVtx[vtxIdx[1]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType uy = g_meshVtx[vtxIdx[1]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType uz = g_meshVtx[vtxIdx[1]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType vx = g_meshVtx[vtxIdx[2]*prm.dim+0] - g_meshVtx[vtxIdx[0]*prm.dim+0];
    ValueType vy = g_meshVtx[vtxIdx[2]*prm.dim+1] - g_meshVtx[vtxIdx[0]*prm.dim+1];
    ValueType vz = g_meshVtx[vtxIdx[2]*prm.dim+2] - g_meshVtx[vtxIdx[0]*prm.dim+2];
    ValueType nx = uy*vz - uz*vy;
    ValueType ny = uz*vx - ux*vz;
    ValueType nz = ux*vy - uy*vx;
    ValueType size = sqrt (nx*nx + ny*ny + nz*nz);
    mappedNormals[ofsNormal++] = nx/size;
    mappedNormals[ofsNormal++] = ny/size;
    mappedNormals[ofsNormal++] = nz/size;
    integrationElements[el*nPoints+pt] = size;
}


__kernel void clIntegrate (
	MeshParam prm,
	__global const ValueType *g_meshVtx,
	__global const int *g_meshIdx,
	__global const ValueType *g_globalTrialPoints,
	__global const ValueType *g_globalTestPoints,
	__global const ValueType *g_globalTrialNormals,
	__global const ValueType *g_trialIntegrationElements,
	__global const ValueType *g_testIntegrationElements,
	__global const ValueType *g_trialValues,
	__global const ValueType *g_testValues,
	__global const ValueType *g_trialWeights,
	__global const ValueType *g_testWeights,
	int trialPointCount,
	int testPointCount,
	int trialComponentCount,
	int testComponentCount,
	int trialDofCount,
	int testDofCount,
	int elementCount,
	int indexAIsTest,
	__global const int *g_elementIndicesA,
	int indexB,
	__global ValueType *g_result
)
{
    int id = get_global_id(0);
    if (id >= elementCount)
        return;

    int i;
    int indexA = g_elementIndicesA[id];
    int testPointIdx, trialPointIdx;

    ValueType geomA[9]; // 9: dim * vertices per element
    ValueType geomB[9]; // 9: dim * vertices per element
    ValueType globalTrialPoint[3];
    ValueType globalTestPoint[3];
    ValueType globalTrialNormal[3];
    ValueType *trialGeom;
    ValueType *testGeom;
    ValueType kval;

    devGetGeometry (g_meshVtx, prm.nvtx, g_meshIdx, prm.dim, prm.nidx, indexA, geomA);
    devGetGeometry (g_meshVtx, prm.nvtx, g_meshIdx, prm.dim, prm.nidx, indexB, geomB);

    __global const ValueType *trialIntegrationElements;
    __global const ValueType *testIntegrationElements;
    __global const ValueType *trialValues;
    __global const ValueType *testValues;

    if (indexAIsTest) {
        testPointIdx = id*testPointCount*prm.dim;
	trialPointIdx = 0;
        testGeom = geomA;
	trialGeom = geomB;
        trialIntegrationElements = g_trialIntegrationElements;
	testIntegrationElements = g_testIntegrationElements + id*testPointCount;
	trialValues = g_trialValues;
	testValues = g_testValues + id*testComponentCount*testDofCount*testPointCount;
    } else {
        testPointIdx = 0;
	trialPointIdx = id*trialPointCount*prm.dim;
        testGeom = geomB;
	trialGeom = geomA;
        trialIntegrationElements = g_trialIntegrationElements + id*trialPointCount;
	testIntegrationElements = g_testIntegrationElements;
	trialValues = g_trialValues + id*trialComponentCount*trialDofCount*trialPointCount;
	testValues = g_testValues;
    }

    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (int testDof = 0; testDof < testDofCount; ++testDof)
	{
	    ValueType sum = 0;
	    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
	        for (i = 0; i < prm.dim; i++) {
		    globalTrialPoint[i] = g_globalTrialPoints[trialPointIdx+trialPoint*prm.dim+i];
		    globalTrialNormal[i] = g_globalTrialNormals[trialPointIdx+trialPoint*prm.dim+i];
		}
	        for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
		    for (i = 0; i < prm.dim; i++)
		        globalTestPoint[i] = g_globalTestPoints[testPointIdx+testPoint*prm.dim+i];
		    kval = devKerneval (globalTestPoint,globalTrialPoint,globalTrialNormal,prm.dim);
		    // currently no dependency on trialDim,testDim
		    for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
		        for (int testDim = 0; testDim < testComponentCount; ++testDim) {
			    sum += g_testWeights[testPoint] *
			        testIntegrationElements[testPoint] *
			        testValues[testDim + testDof*testComponentCount + testPoint*testComponentCount*testDofCount] *
			        kval *
			        trialValues[trialDim + trialDof*trialComponentCount + trialPoint*trialComponentCount*trialDofCount] *
			        trialIntegrationElements[trialPoint] *
			        g_trialWeights[trialPoint];
			}
		}
	    }
	    g_result[testDof + trialDof*testDofCount + id*testDofCount*trialDofCount] = sum;
	}
}
