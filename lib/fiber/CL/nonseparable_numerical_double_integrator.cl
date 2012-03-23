// -*-C++-*-

/**
 * \file integrate_pairs.cl
 * CL code for integrating and array of element pairs
 */

/**
 * \brief Integration kernel
 *
 * This version performs integrals over a list of triangle pairs.
 * \param g_meshvtx Full mesh vertex list
 * \param nmeshvtx Size of vertex list
 * \param g_meshels Full mesh element index list
 * \param nmeshels Size of element index list
 * \param g_tri1 List of triangle indices for the first element of each pair
 * \param g_tri2 List of triangle indices for the second element of each pair
 * \param geometryPairCount Number of element pairs
 * \param g_refpt1 Reference quadrature points for first element
 * \param g_refpt2 Reference quadrature points for second element
 * \param g_weight %Quadrature weights
 * \param nquad Number of quadrature points
 * \param g_basisf Array of basis function values for reference element
 *    [nquad*nquad*3*3]
 * \param g_val Array of integration results [ntri*3*3]
 * \note The following device functions are required to be present in the
 *   CL program:
 *   - devDetJac(const Point3*)
 *   - devMapPoint(const Point3 *vtx,const Point3 *refpt,Point3 *pt)
 *   - devKerneval(const Point3 *pt1, const Point3 *pt2)
 */
__kernel void clIntegratePairs (
    __global const Point3 *g_meshvtx,
    int nmeshvtx,
    __global const Triangle *g_meshels,
    int nmeshels,
    __global const int *g_tri1,
    __global const int *g_tri2,
    int geometryPairCount,
    __global const Point3 *g_refpt1,
    __global const Point3 *g_refpt2,
    __global const FLOAT *g_weight,
    int nquad,
    __global const FLOAT *g_basisf,
    __global FLOAT *g_val
    )
{
    int i, j, qpt;
    int pairIndex = get_global_id(0);
    if (pairIndex >= geometryPairCount) return;

    Point3 vtx1[3];  // vertices of first triangle
    Point3 vtx2[3];  // vertices of second triangle
    Point3 refpt1, refpt2; // quadrature points on reference element
    Point3 pt1, pt2;       // mapped quadrature points

    int tri1 = g_tri1[pairIndex]; // index of first triangle
    int tri2 = g_tri2[pairIndex]; // index of second triangle

    // copy vertices from the global list
    for (i = 0; i < 3; i++) {
        vtx1[i] = g_meshvtx[g_meshels[tri1].idx[i]];
	vtx2[i] = g_meshvtx[g_meshels[tri2].idx[i]];
    }

    // basis function mapper
    FLOAT d1 = devDetJac (vtx1);
    FLOAT d2 = devDetJac (vtx2);
    FLOAT d12 = d1*d2;
    FLOAT val[3*3];
    FLOAT kval;

    for (i = 0; i < 9; i++)
        val[i] = (FLOAT)0;

    for (qpt = 0; qpt < nquad; qpt++) {
        refpt1 = g_refpt1[qpt];
	refpt2 = g_refpt2[qpt];

	// map points from reference to global coordinates
	devMapPoint (vtx1, &refpt1, &pt1);
	devMapPoint (vtx2, &refpt2, &pt2);

	// evaluate kernel
	kval = devKerneval (&pt1, &pt2) * d12 * g_weight[qpt];

	// fill the local matrix
	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
	        val[j + i*3] += kval * g_basisf[j + i*3 + qpt*9];
    }

    // copy results back to global matrix
    for (i = 0; i < 9; i++)
        g_val[pairIndex*9 + i] = val[i];
}
