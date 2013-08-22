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

//#define USE_AMOS_SPECIAL_FUNCTIONS

#if defined USE_AMOS_SPECIAL_FUNCTIONS
#include "amos.hpp"
#endif

namespace Bempp
{

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

template <typename ValueType>
ValueType pi()
{
	return boost::math::constants::pi<ValueType>();
}

template <typename ValueType> // must be a real type
ValueType diff_legendre_p(int n, ValueType x);

template <typename ValueType>  // must be a real type
void legendre_roots(unsigned int N, ValueType *roots, ValueType *weights);

template <typename ValueType>
void FmmHighFreq<ValueType>::generateGaussPoints()
{
	//m_k0 = ValueType(5.0)*getI<ValueType>();
/*	std::vector<lebedev_point_t> lebedev = getLebedevSphere(this->m_P);

	for (unsigned int p=0; p<this->m_P; p++) {
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
	CoordinateType wphi = 2*pi<CoordinateType>()/(2*m_L+1);
	for (unsigned int l=0; l<=2*m_L; l++) {
		cosphi[l] = cos(2*pi<CoordinateType>()*l/(2*m_L+1));
		sinphi[l] = sin(2*pi<CoordinateType>()*l/(2*m_L+1));
	}
	// form pairs of Gauss points along (theta, phi)
	unsigned int ind = 0;
	for (unsigned int l=0; l<m_L; l++) {
		for (unsigned int m=0; m<=2*m_L; m++) {
			this->m_s[3*ind+0] = sintheta[l]*cosphi[m];	// x
			this->m_s[3*ind+1] = sintheta[l]*sinphi[m];	// y
			this->m_s[3*ind+2] = costheta[l];			// z
			this->m_w[ind] = wtheta[l]*wphi;
			ind++;
		}
	}
}

template <typename ValueType>
arma::Col<ValueType> FmmHighFreq<ValueType>::M2M(
		CoordinateType xold[3], CoordinateType xnew[3]) const
{
	arma::Col<ValueType> T(this->m_P);//, std::complex<double>(0.0, 0.0));

	CoordinateType R[3] = {xnew[0]-xold[0], xnew[1]-xold[1], xnew[2]-xold[2]};
	for (unsigned int p=0; p<this->m_P; p++) { // each quadrature point
		CoordinateType rcostheta = R[0]*this->m_s[3*p+0]
			 + R[1]*this->m_s[3*p+1] + R[2]*this->m_s[3*p+2];
		T[p] = exp(-m_kappa*rcostheta);
	} // for each quadrature point
	return T;
}

template <typename ValueType>
arma::Col<ValueType> FmmHighFreq<ValueType>::L2L(
		CoordinateType xold[3], CoordinateType xnew[3]) const
{
	return M2M(xold, xnew);
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
arma::Col<ValueType> 
FmmHighFreq<ValueType>::M2L(
		CoordinateType x1[3], CoordinateType x2[3]) const
{
	using namespace boost::math;
	arma::Col<ValueType> T(this->m_P);
	T.fill(0.0);

	CoordinateType xvec[3] = {x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2]};
	CoordinateType r = sqrt(xvec[0]*xvec[0] + xvec[1]*xvec[1] + xvec[2]*xvec[2]);
	CoordinateType Rhat[3] = {xvec[0]/r, xvec[1]/r, xvec[2]/r};
//	std::cout << "r = " << r << std::endl;
	ValueType i = getI<ValueType>();
//	typename ScalarTraits<ValueType>::ComplexType ComplexType;
//	ComplexType i(0., 1.);

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
		hl = sqrt(pi<CoordinateType>()/(CoordinateType(2)*z))
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
#else // use boost special function, only works for purely real or imag kappa
		if (real(m_kappa)==0) { // purely real
			CoordinateType z = -imag(m_kappa)*r;
			hl = sph_bessel(l, z) + i*sph_neumann(l, z);
		} else if (imag(m_kappa)==0 && real(m_kappa)>0) { // purely imaginary and decaying
			CoordinateType zi = real(m_kappa)*r;
			hl = -sqrt(ValueType(2)/(zi*pi<CoordinateType>()))*pow(i,-l)
				*cyl_bessel_k(CoordinateType(l+0.5), zi);
		} else {
			throw std::invalid_argument("FmmHighFreq::M2L(x1, x2): "
     			"boost special functions only support purely real or imaginary args");
		}
#endif
		ValueType scaledhl = -m_kappa/(16*pi<CoordinateType>()
				*pi<CoordinateType>())*pow(i,l)*ValueType(2*l+1)*hl;

		for (unsigned int p=0; p<this->m_P; p++) { // each quadrature point
			CoordinateType costheta = Rhat[0]*this->m_s[3*p+0] 
				+ Rhat[1]*this->m_s[3*p+1] + Rhat[2]*this->m_s[3*p+2];
			if(costheta> 1.0) costheta =  1.0;
			if(costheta<-1.0) costheta = -1.0;
			T[p] += scaledhl*legendre_p(l, costheta);
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

	roots[int(N/2)] = 0; // take advantage of symmetry
	for (unsigned int i = 1; i <= int(N/2); i++) {
		ValueType x, x1;
		x = cos(pi<ValueType>() * (i - .25) / (N + .5)); // guess root position
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
