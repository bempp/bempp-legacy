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

#include "fmm_black_box.hpp"
#include "fmm_transform.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/math/constants/constants.hpp>
#include <complex>
#include <iostream>		// std::cout
#include <limits>

namespace Bempp
{

template <typename ValueType>
void FmmBlackBox<ValueType>::generateGaussPoints()
{
	// calculate the Chebyshev nodes, the N zeros of T_N(x)
	CoordinateType nodes[m_n];
	for (unsigned int m = 0; m < m_n; m++) {
		nodes[m] = cos( M_PI*CoordinateType(2*m+1)/CoordinateType(2*m_n) );
	}

	// form nodes in R^3 by tensor product
	unsigned int ind = 0;
	for (unsigned int i1 = 0; i1 < m_n; i1++) {
		for (unsigned int i2 = 0; i2 < m_n; i2++) {
			for (unsigned int i3 = 0; i3 < m_n; i3++) {
				this->m_quadraturePoints(0, ind) = nodes[i1]; // x
				this->m_quadraturePoints(1, ind) = nodes[i2]; // y
				this->m_quadraturePoints(2, ind) = nodes[i3]; // z
				this->m_quadratureWeights(ind) = 1.;
				ind++;
			} // for each node along z
		} // for each node along y
	} // for each node along x

	// compute T_k(x) for k between 0 and N-1 inclusive.
	for (unsigned int m = 0; m < m_n; m++) {
		m_Tk(m, 0) = 1;
		m_Tk(m, 1) = nodes[m];
		for (unsigned int k = 2; k < m_n; k++) {
			m_Tk(m, k) = 2*nodes[m]*m_Tk(m, k-1) - m_Tk(m, k-2);
		}
	}

} // FmmBlackBox<ValueType>::generateGaussPoints()

template <typename ValueType>
arma::Mat<ValueType> FmmBlackBox<ValueType>::M2M(
		const arma::Col<CoordinateType>& childPosition,
		const arma::Col<CoordinateType>& parentPosition) const
{
	// change to use interpolation at this stage in future
	arma::Mat<ValueType> T(this->quadraturePointCount(), 
		this->quadraturePointCount());

	arma::Col<CoordinateType> xvec = childPosition - parentPosition;

	// source point positions along x, y and z relative to parent in [-1,1]
	arma::Mat<CoordinateType> sourcePos(getN(), 3);
	const CoordinateType prec = std::numeric_limits<CoordinateType>::denorm_min();
	for (unsigned int dim = 0; dim < xvec.n_rows; dim++) {
		int childOffset = (xvec(dim) > prec) - (xvec(dim) < -prec); // sign
		sourcePos.col(dim) = (m_Tk.col(1) + childOffset)/2;  // T_1 are nodes
	}

	// calculate S(m,n) = \sum T_k(x_m) sources_k(x_n), independently for 
	// each dimension using Clenshaw algorithm to avoid forming T_k(x_s) explicitly
	CoordinateType S[getN()][sourcePos.n_rows][3];

	CoordinateType b[getN() + 2]; // reverse coefficients
	b[getN()] = b[getN() + 1] = 0;

	for (unsigned int dim = 0; dim < xvec.n_rows; dim++) {
		for (unsigned int m = 0; m < getN(); m++) {
			// might in general have more sources than m_n
			for (unsigned int n = 0; n < sourcePos.n_rows; n++) {
				// f(x) = a_0 T_0(x) + a_1 T_1(x) + a_2 T_2(x) + ...
				// N/2 f(x_s) = 1/2 T_0(x_m) T_0(x_s) + T_1(x_m) T_1(x_s) + ...
				// so set a_k = T_k(x_m) and halve first term in sum below
				CoordinateType x = sourcePos(n, dim);

				for (unsigned int k = getN()-1; k > 0; k--) { // sum over order k
					b[k] = m_Tk(m, k) + 2*x*b[k+1] - b[k+2];
				}
				// 0.5 due to (T_m, T_n) = \pi/2 (1+\delta_{0n})\delta_{mn}
				S[m][n][dim] = (0.5*m_Tk(m, 0) + x*b[1] - b[2])*2/getN();
			} // for each source point
		} // for each m
	} // for each dim

	// take tensor product along each dimension
	CoordinateType n3 = getN()*getN()*getN();
	arma::Mat<ValueType> S3(n3, n3);
	unsigned int m = 0;
	for (unsigned int m1 = 0; m1 < getN(); m1++) {
		for (unsigned int m2 = 0; m2 < getN(); m2++) {
			for (unsigned int m3 = 0; m3 < getN(); m3++) {
				unsigned int n = 0;
				for (unsigned int n1 = 0; n1 < getN(); n1++) {
					for (unsigned int n2 = 0; n2 < getN(); n2++) {
						for (unsigned int n3 = 0; n3 < getN(); n3++) {
							S3(m, n) = S[m1][n1][0]*S[m2][n2][1]*S[m3][n3][2];
							n++;
						} // for n3
					} // for n2
				} // for n1
				m++;
			} // for m3
		} // for m2
	} // for m1

	return S3;
}

template <typename ValueType>
arma::Mat<ValueType> FmmBlackBox<ValueType>::L2L(
		const arma::Col<CoordinateType>& parentPosition,
		const arma::Col<CoordinateType>& childPosition) const
{
	 // argument order swapped, L2L = M2M' exactly
	return M2M(childPosition, parentPosition).st();
}


template <typename ValueType>
ValueType 
evalKernel(
	const arma::Col<ValueType>& fieldPos, 
	const arma::Col<ValueType>& sourcePos)
{
	arma::Col<ValueType> R = fieldPos - sourcePos;
	return 1./(4*M_PI*norm(R, 2));
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmBlackBox<ValueType>::M2L(
		const arma::Col<CoordinateType>& sourceCentre, // loops over interaction list
		const arma::Col<CoordinateType>& fieldCentre,  // origin
		const arma::Col<CoordinateType>& boxSize) const
{
	arma::Mat<CoordinateType> fieldNodes (getN(), 3);
	arma::Mat<CoordinateType> sourceNodes(getN(), 3);

	for (unsigned int dim = 0; dim < boxSize.n_rows; dim++) {
		fieldNodes.col (dim) = (m_Tk.col(1)*boxSize(dim)/2 + fieldCentre (dim));
		sourceNodes.col(dim) = (m_Tk.col(1)*boxSize(dim)/2 + sourceCentre(dim));
	}

	CoordinateType n3 = getN()*getN()*getN();
	arma::Mat<ValueType> T3(n3, n3);

	arma::Col<CoordinateType> fieldPos (fieldCentre.n_rows);
	arma::Col<CoordinateType> sourcePos(sourceCentre.n_rows);
	unsigned int m = 0;
	for (unsigned int m1 = 0; m1 < getN(); m1++) {
		fieldPos(0) = fieldNodes(m1, 0);
		for (unsigned int m2 = 0; m2 < getN(); m2++) {
			fieldPos(1) = fieldNodes(m2, 1);
			for (unsigned int m3 = 0; m3 < getN(); m3++) {
				fieldPos(2) = fieldNodes(m3, 2);
				unsigned int n = 0;
				for (unsigned int n1 = 0; n1 < getN(); n1++) {
					sourcePos(0) = sourceNodes(n1, 0);
					for (unsigned int n2 = 0; n2 < getN(); n2++) {
						sourcePos(1) = sourceNodes(n2, 1);
						for (unsigned int n3 = 0; n3 < getN(); n3++) {
							sourcePos(2) = sourceNodes(n3, 2);
							T3(m, n) = evalKernel(fieldPos, sourcePos);
							n++;
						} // for n3
					} // for n2
				} // for n1
				m++;
			} // for m3
		} // for m2
	} // for m1

	return T3;
}


// FmmSingleLayerHighFreqFMM

template <typename ValueType>
void FmmSingleLayerBlackBox<ValueType>::evaluateTrial(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& multipole, // [-1,1]
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	arma::Col<CoordinateType> pointScaled = 2.*(point - nodeCentre) / nodeSize;

	// Gauss points can lie outside a given node, since triangles can intesect the faces.
	// Should ideally enlarge the box so that all mulipoles fully enclose Gauss points.
	// If this is not done, then end up with extrapolation, which increases the error

	// should optimise by only calculating matrix S along a single dimension of 
	// the mulipole coefficients and reusing the information, but where to calculate
	// (be careful since parallelised)
	CoordinateType S = 1.0;
	CoordinateType TkMultipole[this->getN()];
	CoordinateType TkPoint[this->getN()];
	for (unsigned int dim=0; dim<point.n_rows; dim++) {
		TkMultipole[0] = 1;
		TkMultipole[1] = multipole[dim];
		for (unsigned int j = 2; j < this->getN(); j++) {
			TkMultipole[j] = 2*multipole[dim]*TkMultipole[j-1] - TkMultipole[j-2];
		}
		TkPoint[0] = 1;
		TkPoint[1] = pointScaled[dim];
		for (unsigned int j = 2; j < this->getN(); j++) {
			TkPoint[j] = 2*pointScaled[dim]*TkPoint[j-1] - TkPoint[j-2];
		}
		CoordinateType localS = TkMultipole[0]*TkPoint[0];
		for (unsigned int j = 1; j < this->getN(); j++) {
			localS += 2*TkMultipole[j]*TkPoint[j];
		}
		S *= localS/(this->getN());
	}

	result(0) =  S;
}

// should really create a third function that both evaluateTest and evaluateTrial call
template <typename ValueType>
void FmmSingleLayerBlackBox<ValueType>::evaluateTest(
			const arma::Col<CoordinateType>& point,
			const arma::Col<CoordinateType>& normal,
			const arma::Col<CoordinateType>& multipole,
			const arma::Col<CoordinateType>& nodeCentre,
			const arma::Col<CoordinateType>& nodeSize,
			arma::Col<ValueType>& result) const
{
	evaluateTrial(point, normal, multipole, nodeCentre, nodeSize, result);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmBlackBox);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmSingleLayerBlackBox);

} // namespace Bempp
