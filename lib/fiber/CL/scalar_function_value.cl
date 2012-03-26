/**
 * \brief Evaluate linear basis function on a reference triangle
 * \param [in] p pointer to reference point
 * \param [out] bf array of basis function values [3]
 * \note Definition of reference triangle: v = {{0,0,0},{1,0,0},{0,1,0}}
 * \note z-coordinates of all input points must be 0
 */
void devExpressionEvaluate (const ValueType *basisData,
	              ValueType *result,
		      int n)
{
    for (int i = 0; i < n; i++)
	result[i] = basisData[i];
}
