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

#include "fmm_transform.hpp"
#include "lebedev.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <complex>
#include <iostream>		// std::cout
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>

#define USE_AMOS_SPECIAL_FUNCTIONS

#if defined USE_AMOS_SPECIAL_FUNCTIONS
#include "amos.hpp"
#endif

namespace Bempp
{

template <typename ValueType> // must be a real type
ValueType diff_legendre_p(int n, ValueType x);

template <typename ValueType>  // must be a real type
void legendre_roots(unsigned int N, ValueType *roots, ValueType *weights);

template <typename ValueType>
void FmmHighFreq<ValueType>::generateGaussPoints()
{
	CoordinateType pi = boost::math::constants::pi<CoordinateType>();

/*	std::vector<lebedev_point_t> lebedev = getLebedevSphere(this->quadraturePointCount());

	for (unsigned int p=0; p < this->quadraturePointCount(); p++) {
		this->m_s[3*p+0] = lebedev[p].x;
		this->m_s[3*p+1] = lebedev[p].y;
		this->m_s[3*p+2] = lebedev[p].z;
		this->m_w[p] = lebedev[p].w;
	}
	return;
*/
	// form Gauss-Legendre quadrature m_fmmTransformpoints along x= cos(theta)
	CoordinateType costheta[m_L]; // roots along x = cos(theta)
	CoordinateType sintheta[m_L];
	CoordinateType wtheta[m_L];
	legendre_roots(m_L, costheta, wtheta);
	for (unsigned int l=0; l<m_L; l++) {
		sintheta[l] = sqrt(1-costheta[l]*costheta[l]);
	}

	// form regular grid along phi (periodic function)
	CoordinateType cosphi[2*m_L+1];
	CoordinateType sinphi[2*m_L+1];
	CoordinateType wphi = 2*pi/(2*m_L+1);
	for (unsigned int l=0; l<=2*m_L; l++) {
		cosphi[l] = cos(2*pi*l/(2*m_L+1));
		sinphi[l] = sin(2*pi*l/(2*m_L+1));
	}
	// form pairs of Gauss points along (theta, phi)
	unsigned int ind = 0;
	for (unsigned int l=0; l<m_L; l++) {
		for (unsigned int m=0; m<=2*m_L; m++) {
			this->m_quadraturePoints(0, ind) = sintheta[l]*cosphi[m];	// x
			this->m_quadraturePoints(1, ind) = sintheta[l]*sinphi[m];	// y
			this->m_quadraturePoints(2, ind) = costheta[l];			// z
			this->m_quadratureWeights(ind)   = wtheta[l]*wphi;
			ind++;
		}
	}
}

template <typename ValueType>
arma::Mat<ValueType> FmmHighFreq<ValueType>::M2M(
		const arma::Col<CoordinateType>& xold, 
		const arma::Col<CoordinateType>& xnew) const
{
	arma::Col<ValueType> T(this->quadraturePointCount());

	arma::Col<CoordinateType> R = xnew - xold;
	arma::Col<CoordinateType> khat(3);
	for (unsigned int p=0; p < this->quadraturePointCount(); p++) {
		khat = this->getQuadraturePoint(p);
		CoordinateType rcostheta = dot(R, khat);
		T[p] = exp(-m_kappa*rcostheta);
	} // for each quadrature point
	return T;
}

template <typename ValueType>
arma::Mat<ValueType> FmmHighFreq<ValueType>::L2L(
		const arma::Col<CoordinateType>& xold, 
		const arma::Col<CoordinateType>& xnew) const
{
	return M2M(xold, xnew);
}

template <typename ValueType>
ValueType getI()
{
	throw std::invalid_argument("getI(): "
		"can only be called for complex result types");
}
template <>
std::complex<float> getI()
{
	return std::complex<float>(0.0, 1.0);
}
template <>
std::complex<double> getI()
{
	return std::complex<double>(0.0, 1.0);
}
double real(double x)
{
	return x;
}
double imag(double x)
{
	return 0.0;
}

template <typename ValueType>
arma::Mat<ValueType> 
FmmHighFreq<ValueType>::M2L(
		const arma::Col<CoordinateType>& x1, 
		const arma::Col<CoordinateType>& x2) const
{
	using namespace boost::math;
	CoordinateType pi = boost::math::constants::pi<CoordinateType>();

	arma::Col<ValueType> T(this->quadraturePointCount());
	T.fill(0.0);

	arma::Col<CoordinateType> xvec = x2 - x1;
	CoordinateType r = norm(xvec, 2);
	const arma::Col<CoordinateType>& Rhat = xvec/r;

	ValueType i = getI<ValueType>();

	for (unsigned int l=0; l<=m_L; l++) {
		ValueType hl;

#if defined USE_AMOS_SPECIAL_FUNCTIONS
		ValueType z = i*m_kappa*r;
		double zr = real(z);
		double zi = imag(z);
		double nu = l+0.5;
		int kode = 1;
		int N = 1;
		int kind = 1;

		double cyr,cyi; 	// Output values
		int nz,ierr;

		amos::zbesh(&zr,&zi,&nu,&kode,&kind,&N,&cyr,&cyi,&nz,&ierr);
		hl = sqrt(pi/(CoordinateType(2)*z))
			*(CoordinateType(cyr)+i*CoordinateType(cyi));

		std::string amosErrorMessages[6] = {
			"IERR=0, NORMAL RETURN - COMPUTATION COMPLETED",
			"IERR=1, INPUT ERROR   - NO COMPUTATION",
			"IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS"
			"        TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH",
			"IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE"
			"        BUT LOSSES OF SIGNIFCANCE BY ARGUMENT"
			"        REDUCTION PRODUCE LESS THAN HALF OF MACHINE"
			"        ACCURACY",
			"IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-"
			"        TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-"
			"        CANCE BY ARGUMENT REDUCTION",
			"IERR=5, ERROR              - NO COMPUTATION,"
			"        ALGORITHM TERMINATION CONDITION NOT MET"
		};
		if (ierr != 0) {
			throw std::invalid_argument(std::string("FmmHighFreq::M2L(x1, x2): "
				"AMOS: ")+amosErrorMessages[ierr] );
		}
#else // use boost special functions, only works for purely real or imag kappa
		if (real(m_kappa)==0) { // purely imaginary
			CoordinateType z = -imag(m_kappa)*r;
			hl = sph_bessel(l, z) + i*sph_neumann(l, z);
		} else if (imag(m_kappa)==0 && real(m_kappa)>0) { // purely real and decaying
			CoordinateType zi = real(m_kappa)*r;
			hl = -sqrt(ValueType(2)/(zi*pi))*pow(i,-l)
				*cyl_bessel_k(CoordinateType(l+0.5), zi);
		} else {
			throw std::invalid_argument("FmmHighFreq::M2L(x1, x2): "
     			"boost special functions only support purely real or imaginary args");
		}
#endif
		ValueType scaledhl = -m_kappa/(16*pi*pi)*pow(i,l)*CoordinateType(2*l+1)*hl;

		arma::Col<CoordinateType> khat(3); // slow! do not perform in loop!
		for (unsigned int p = 0; p < this->quadraturePointCount(); p++) {
			khat = this->getQuadraturePoint(p);
			CoordinateType costheta = dot(Rhat, khat);
			if(costheta> 1.0) costheta =  1.0;
			if(costheta<-1.0) costheta = -1.0;
			T(p) += scaledhl*legendre_p(l, costheta);
		}
	}
	return T;
}

// P'_n(x): Derivative of the n_{th} order Legendre polynomial w.r.t. x
template <typename ValueType> // must be a real type
ValueType diff_legendre_p(int n, ValueType x)
{
	using namespace boost::math;
	return n*( x*legendre_p(n, x) - legendre_p(n-1, x) ) / (x*x - 1);
}

// from http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
// find the Gauss Legendre quadrature points, P_n(xi) = 0
template <typename ValueType>  // must be a real type
void legendre_roots(unsigned int N, ValueType *roots, ValueType *weights)
{
	using namespace boost::math;
	ValueType pi = boost::math::constants::pi<ValueType>();

	roots[int(N/2)] = 0; // take advantage of symmetry
	for (unsigned int i = 1; i <= int(N/2); i++) {
		ValueType x, x1;
		x = cos(pi * (i - .25) / (N + .5)); // guess root position
		int iter=100;
		do { // apply Newton-Raphson method to find roots
			x1 = x;
			x -= legendre_p(N, x) / diff_legendre_p(N, x);
		} while (x != x1 && --iter); // well-behaved function, convergence guaranteed
		roots[i-1] =  x;
		roots[N-i] = -x;
 	}

	for (unsigned int i = 1; i <= int(N/2)+1; i++) {
		ValueType x = roots[i-1];
		ValueType diffPx = diff_legendre_p(N, x);
		weights[i-1] = 2 / ((1 - x*x) * diffPx*diffPx);
		weights[N-i] = weights[i-1];
	}

/*	for (unsigned int i = 1; i <= N; i++)
		std::cout << roots[i - 1] << ", ";
	std::cout <<std::endl;
 
	for (unsigned int i = 1; i <= N; i++)
		std::cout << weights[i - 1] << ", ";
	std::cout <<std::endl;
*/
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmTransform);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmHighFreq);

} // namespace Bempp
