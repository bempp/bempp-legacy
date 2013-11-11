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

#include "fmm_cache.hpp"
#include "fmm_transform.hpp"
#include "octree.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>		// std::cout
#include <complex>
#include <string>
#include <sstream>


namespace Bempp
{

template <typename ValueType>
FmmCache<ValueType>::FmmCache(
		const FmmTransform<ValueType>& fmmTransform,
		unsigned int levels)
	: m_fmmTransform(fmmTransform), m_levels(levels), m_topLevel(2)
{
}

template <typename ValueType>
//template <typename KernelType>
void 
FmmCache<ValueType>::initCache(
		const arma::Col<CoordinateType> &lowerBound,
		const arma::Col<CoordinateType> &upperBound)//,
//		const Fiber::CollectionOfKernels<KernelType>& kernels)
{
	arma::Col<CoordinateType> origin(3);
	origin.fill(0.);

	m_cacheM2L.resize(m_levels-m_topLevel+1);

	for (unsigned int level = m_topLevel; level<=m_levels; level++) {

		// there are 7^3-3^3=316 unique translation matrices for translation
		// invariant operators
		m_cacheM2L[level-m_topLevel].resize(316);

		unsigned int boxesPerSide = getNodesPerSide(level);
		arma::Col<CoordinateType> boxSize;
		boxSize = (upperBound - lowerBound)/boxesPerSide;

		arma::Col<CoordinateType> centre(3);

		unsigned int index = 0;
		for (int indx=-3; indx<=3; indx++) {
			centre[0] = indx*boxSize[0];
			for (int indy=-3; indy<=3; indy++) {
				centre[1] = indy*boxSize[1];
				for (int indz=-3; indz<=3; indz++) {
					centre[2] = indz*boxSize[2];

					if (abs(indx) > 1 || abs(indy) > 1 || abs(indz) > 1) {

						arma::Mat<ValueType> m2l = m_fmmTransform.M2L(centre, origin, boxSize, level);

						m_cacheM2L[level-m_topLevel][index++] = m2l;

					} // if not a nearest neighbour
				} // for each x offset
			} // for each x offset
		} // for each x offset
	} // for each level

	// M2M & L2L cache
	m_cacheM2M.resize(m_levels-m_topLevel);
	m_cacheL2L.resize(m_levels-m_topLevel);
	for (unsigned int level = m_topLevel; level<m_levels; level++) { //2->levels-1
		m_cacheM2M[level-m_topLevel].resize(8);
		m_cacheL2L[level-m_topLevel].resize(8);

		unsigned int boxesPerSide = getNodesPerSide(level);
		arma::Col<CoordinateType> boxSize;
		boxSize = (upperBound - lowerBound)/boxesPerSide;

		arma::Col<CoordinateType> Rchild(3); // child pos

		for (unsigned long child = 0;	child < 8; child++) {

			unsigned long ind[3];
			deMorton(&ind[0], &ind[1], &ind[2], child);

			Rchild(0) = ind[0]; Rchild(1) = ind[1]; Rchild(2) = ind[2];
			Rchild = (Rchild-0.5) % boxSize/2;

			arma::Mat<ValueType> m2m = m_fmmTransform.M2M(Rchild, origin, level);

			m_cacheM2M[level-m_topLevel][child] = m2m;

			arma::Mat<ValueType> l2l = m_fmmTransform.L2L(origin, Rchild, level);

			m_cacheL2L[level-m_topLevel][child] = l2l;
		} // for each child
	} // for each level

	compressM2L(true);
} // initCache


// For the bbFFM the M2L operator will be block Toeplitz (fully Toeplitz in 1D)
// when the kernel is symmetric. N.B that the order in which the interaction blocks
// are stored does not matter. However, if symmetry is to be exploited (so that one
// SVD suffices), the source interaction blocks must be arranged symmetrically about
// the field centre, e.g. -4 -3 x x f x x 3 4 not -4 -3 x x f x x 3 4 5.
// If symmetry is to be exploited note that Ufat_k'*Vthin_k = \pm 1. However, since
// we calculate multipoleCoefficients = Ufat*(Ufat'*K*Vthin)*Vthin'*multipoleCoefficients
// the \pm 1 does not matter. For symmetric kernels it is thus safe to assume
// Vthin = Ufat.
// If checking Ufat_k'*Vthin_k = \pm 1, note that if kernel is isotropic in x, y 
// and z, then triplet vectors result with the same sigma. Due to this degeneracy, the 
// vectors do not satisfy Ufat_k'*Vthin_k = \pm 1. If one scales x, y, and z slightly 
// differently, then the degeneracy is broken, and the equality is perfectly satisfied.
template <typename ValueType>
void 
FmmCache<ValueType>::compressM2L(bool isSymmetric)
{
	if ( !m_fmmTransform.isCompressedM2L() ) {
		return;
	}
	std::cout << "Compressing M2L matrices by SVD" << std::endl;

	arma::Mat<ValueType> kernelWeightMat;
	m_fmmTransform.getKernelWeight(kernelWeightMat, m_kernelWeightVec);
	
	int npt = m_fmmTransform.quadraturePointCount();

	m_Ufat.resize(m_levels-m_topLevel+1);
	m_Vthin.resize(m_levels-m_topLevel+1);
	arma::Mat<ValueType> Ufat(npt, npt);
	arma::Col<CoordinateType> sigma(npt);
	arma::Mat<ValueType> Vfat;
	arma::Mat<ValueType> kernelsFat(npt, 316*npt);

	for (unsigned int level = m_topLevel; level<=m_levels; level++) {
		// scale all kernel matrices by the weight and copy to flat structure
		for (unsigned int item = 0; item<316; item++) {
			m_cacheM2L[level-m_topLevel][item] %= kernelWeightMat;
			kernelsFat.cols(item*npt, (item+1)*npt-1)
			  = m_cacheM2L[level-m_topLevel][item];
		}
		// Compute the SVD of the scaled fat collection of Kernels
		if ( !arma::svd_econ(Ufat, sigma, Vfat, kernelsFat, 'l') ) {
			throw std::invalid_argument("FmmCache<ValueType>::compressM2L(): "
				"singular value decomposition failed");
		}

		// store the first few columns of U used for the reduced rank approximation
		// into m_Ured
		int cutoff = npt/2 - 1;
		m_Ufat[level-m_topLevel] = Ufat.cols(0, cutoff-1);

		if ( !isSymmetric ) { // if M2L or Kernel is asymmetric
			arma::Mat<ValueType> kernelsThin(316*npt, npt);
			arma::Mat<ValueType> Uthin;
			arma::Mat<ValueType> Vthin(npt, npt);
			for (unsigned int item = 0; item<316; item++) {
				kernelsThin.rows(item*npt, (item+1)*npt-1)
				  = m_cacheM2L[level-m_topLevel][item];
			}

			if ( !arma::svd_econ(Uthin, sigma, Vthin, kernelsThin, 'r') ) {
				throw std::invalid_argument("FmmCache<ValueType>::compressM2L(): "
					"singular value decomposition failed");			
			}
			m_Vthin[level-m_topLevel] = Vthin.cols(0, cutoff-1);
		} else {
			m_Vthin[level-m_topLevel] = m_Ufat[level-m_topLevel];
		}

		// Reduce down the M2L matrices from npt x npt to cutoff x cutoff
		for (unsigned int item = 0; item<316; item++) {
			m_cacheM2L[level-m_topLevel][item] = 
				 m_Ufat [level-m_topLevel].st()*m_cacheM2L[level-m_topLevel][item]
				*m_Vthin[level-m_topLevel];
		}

	} // for each level in octree
}

// call before M2L operation on all nodes on all levels
template <typename ValueType>
void 
FmmCache<ValueType>::compressMultipoleCoefficients(
	arma::Col<ValueType>& mcoefs,
	int level) const
{
	if ( m_fmmTransform.isCompressedM2L() ) {
		mcoefs = m_Vthin[level-m_topLevel].st()*(m_kernelWeightVec % mcoefs);
	}
}

// call after M2L operation on all nodes on all levels
template <typename ValueType>
void 
FmmCache<ValueType>::explodeLocalCoefficients(
	arma::Col<ValueType>& lcoefs,
	int level) const
{
	if ( m_fmmTransform.isCompressedM2L() ) {
		lcoefs = m_kernelWeightVec % (m_Ufat[level-m_topLevel]*lcoefs);
	}
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::M2M(unsigned int level, unsigned int item) const
{
	return m_cacheM2M[level-m_topLevel][item];
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::M2L(unsigned int level, unsigned int item) const
{
	return m_cacheM2L[level-m_topLevel][item];
	//unsigned int boxesPerSide = getNodesPerSide(level);
	//return m_cacheM2L[0][item]*(4./boxesPerSide); // 1/r singularity test
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmCache<ValueType>::L2L(unsigned int level, unsigned int item) const
{
	return m_cacheL2L[level-m_topLevel][item];
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmCache);

} // namespace Bempp
