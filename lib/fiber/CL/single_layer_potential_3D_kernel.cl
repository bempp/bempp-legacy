// -*-C++-*-

/**
 * \file single_layer_potential_3D_kernel.cl
 * OpenCL implementation for single layer potential kernel evaluation
 */

/**
 * \brief Single layer potential evaluation for a single point pair
 * \param testPoint test point coordinates
 * \param trialPoint trial point coordinates
 * \param coordCount number of coordinates for each point
 * \note testPoint and trialPoint must be of size coordCount
 */
ValueType devKerneval (const ValueType *testPoint,
		       const ValueType *trialPoint,
	               const ValueType *trialNormal,
		       int coordCount)
{
    int k;
    ValueType sum, diff;
    sum = 0;
    for (k = 0; k < coordCount; k++) {
        diff = testPoint[k] - trialPoint[k];
	sum += diff*diff;
    }
    return 1.0 / (4.0 * M_PI * sqrt(sum));
}


#ifdef UNDEF // For now

/**
 * \brief Single layer potential evaluation for two elements from a grid
 * \param testPoints array of test points
 * \param trialPoints array of trial points
 * \param result buffer for kernel evaluation results
 * \param testPointCount number of test points
 * \param trialPointCount number of trial points
 * \param coordCount number of coordinates for each point
 * \note testPoints must be of size testPointCount*coordCount
 * \note trialPoints must be of size trialPointCount*coordCount
 * \note result must be large enough to store testPointCount*trialPointCount
 *   values
 */
void devKernevalGrid (__global const ValueType *testPoints,
		      __global const ValueType *trialPoints,
		      __global ValueType *result,
		      int testPointCount,
		      int trialPointCount,
		      int coordCount)
{
    int i, j, k;
    ValueType testPoint[coordCount], trialPoint[coordCount];
    ValueType sum, diff;
    for (i = 0; i < trialPointCount; i++) {
        for (k = 0; k < coordCount; k++)
	    trialPoint[k] = trialPoints[k+i*coordCount];
        for (j = 0; j < testPointCount; j++) {
	    for (k = 0; k < coordCount; k++)
	        testPoint[k] = testPoints[k+j*coordCount];
	        // the test points could be buffered locally outside the loop
	    result[j+i*testPointCount] =
	        devKerneval (testPoint, trialPoint, coordCount);
	}
    }
}


/**
 * \brief Single layer potential evaluation for an element pair
 * \param testPoints array of test points
 * \param trialPoints array of trial points
 * \param result buffer for kernel evaluation results
 * \param pointCount number of test and trial points
 * \param coordCount number of coordinates for each point
 * \note testPoints and trialPoints must be of size pointCount*coordCount
 * \note result must be large enough to store pointCount values
 */
void devKernevalPair (__global const ValueType *testPoints,
		      __global const ValueType *trialPoints,
		      __global ValueType *result,
		      int pointCount,
		      int coordCount)
{
    int i, j;
    ValueType testPoint[coordCount], trialPoint[coordCount];
    ValueType sum, diff;
    for (i = 0; i < pointCount; i++) {
        for (k = 0; k < coordCount; k++) {
	    trialPoint[k] = trialPoints[k+i*coordCount];
	    testPoint[k] = testPoints[k+j*coordCount];
	}
	result[i] = devKerneval (testPoint, trialPoint, coordCount);
    }
}


#endif // UNDEF
