// -*-C++-*-

/**
 * \file laplace_3d_double_layer_potential_kernel.cl
 * OpenCL implementation for double layer potential kernel evaluation
 */

/**
 * \brief Double layer potential evaluation for a single point pair
 * \param testPoint test point coordinates
 * \param trialPoint trial point coordinates
 * \param trialNormal components of the vector normal to the surface at trial
 * point
 * \param coordCount number of coordinates for each point
 * \note testPoint, trialPoint and trialNormal must be of size coordCount
 */
ValueType devKerneval(const ValueType *testPoint, const ValueType *trialPoint,
                      const ValueType *trialNormal, int coordCount) {
  int k;
  ValueType diff, distance, denominatorSum, numeratorSum;
  denominatorSum = 0;
  numeratorSum = 0;
  for (k = 0; k < coordCount; k++) {
    diff = trialPoint[k] - testPoint[k];
    denominatorSum += diff * diff;
    numeratorSum += diff * trialNormal[k];
  }
  distance = sqrt(denominatorSum);
  return -numeratorSum / (4.0 * M_PI * distance * distance * distance);
}
