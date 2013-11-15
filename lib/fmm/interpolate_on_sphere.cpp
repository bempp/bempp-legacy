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

#include "interpolate_on_sphere.hpp"
#include "legendre_roots.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <complex>
#include <vector>		// std::vector
#include <iostream>		// std::cout

namespace Bempp
{

// should be called only for purely real data types
// assuming positive lmax, and |m|<=l
// returns Q_l^m (x) for l=|m|..lmax (stored column-wise), where
// Q_l^m (x) = P_l^m (x)/ sqrt(\int(P_l^m (u)^2 du)
// might adapt later to return Q_l^m (x) and Q_{l+1}^m (x) only
template <typename ValueType>
void normalisedAssociatedLegendre(
	unsigned int lmax, int m_signed,
	const arma::Row<ValueType> &x,
	arma::Mat<ValueType> &Q)
{
	unsigned int m = abs(m_signed);

	if (lmax<m) {
		Q.set_size(0, x.n_cols);
		return;
	}
	Q.set_size(lmax-m+1, x.n_cols);
	Q.fill(0.);

	// sectorial term, l=m
	Q.row(0).fill(1./sqrt(2));
	//Q.row(0) %= (-1)^m*(1-x.^2).^(m/2)/2^m*sqrt(factorial(2*m+1))/factorial(m);

	// avoid factorial in calculating sectorial term, since it is prone to 
	// overflows for large m. TODO: consider caching sectorial for each m
	arma::Row<ValueType> u = sqrt(1-x%x);
	for (unsigned k = 1; k <= m; k++) {
		Q.row(0) %= -sqrt((2.*k+1.)/(2.*k))*u;
	}

	// Q(m+1,m), formula given by setting l-1=m into the full recursive formula
	if (lmax > m) {
		Q.row(1) = sqrt(2*m+3)*(x % Q.row(0));
	}

	// non-sectorial terms, l>m
	for (unsigned int l = m+2; l <= lmax; l++) {
		ValueType alm = sqrt((2.*l-1)*(2*l+1)/((l-m)*(l+m)));
		ValueType blm = sqrt((2.*l+1)*(l+m-1)*(l-m-1)/((l-m)*(l+m)*(2*l-3)));
		Q.row(l-m) = alm*(x % Q.row(l-m-1)) - blm*Q.row(l-m-2);
	}

	if (m_signed < 0) {
		Q *= pow(-1, m);
	}
} // normalisedAssociatedLegendre

// same as above, but just returns Q_n^m, for the three degrees n={l-1,l,l+1}
// Need l-1 case for the evaluation of Q'_l^m(x)|_{x=0 }
template <typename ValueType>
void normalisedAssociatedLegendreThreeDegrees(
	unsigned int l, int m_signed,
	const arma::Row<ValueType> &x,
	arma::Mat<ValueType> &Q)
{
	unsigned int m = abs(m_signed);

	Q.set_size(3, x.n_cols);
	Q.fill(0.);

	if (l+1<m) {
		return;
	}

	// sectorial term, l=m
	Q.row(2).fill(1./sqrt(2));

	// avoid factorial in calculating sectorial term, since it is prone to 
	// overflows for large m. TODO: consider caching sectorial for each m
	arma::Row<ValueType> u = sqrt(1-x%x);
	for (unsigned k = 1; k <= m; k++) {
		Q.row(2) %= -sqrt((2.*k+1.)/(2.*k))*u;
	}

	// Q(m+1,m), formula given by setting l-1=m into the full recursive formula
	if (l+1 > m) {
		Q.row(1) = Q.row(2);
		Q.row(2) = sqrt(2*m+3)*(x % Q.row(1));
	}

	// non-sectorial terms, l>m
	for (unsigned int n = m+2; n <= l+1; n++) {
		ValueType alm = sqrt((2.*n-1)*(2*n+1)/((n-m)*(n+m)));
		ValueType blm = sqrt((2.*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
		Q = join_cols(Q.rows(1, 2), alm*(x % Q.row(2)) - blm*Q.row(1));
	}

	if (m_signed < 0) {
		Q *= pow(-1, m);
	}
} // normalisedAssociatedLegendre


// calculate the normalised associated Legendre polynomials and their derivatives
// check this is correct for -m using sage (gen_legendre_P only supports positive l&m?)
template <typename ValueType>
void normalisedAssociatedLegendre(
	unsigned int lmax, int m_signed,
	const arma::Row<ValueType> &x,
	arma::Mat<ValueType> &Q,
	arma::Mat<ValueType> &diffQ)
{
	normalisedAssociatedLegendre(lmax, m_signed, x, Q);

	unsigned int m = abs(m_signed);

	diffQ.set_size(Q.n_rows, Q.n_cols);
	arma::Row<ValueType> invu2 = 1./(1-x%x);

	// - 1/(1-x^2)*l*x*Q_l^m (x)
	for (int l = m; l <= lmax; l++) {
		diffQ.row(l-m) = -l* (invu2 % x % Q.row(l-m));
	}

	// + 1/(1-x^2)*a*Q_{l-1}^m (x)
	for (unsigned int l = m+1; l <= lmax; l++) {
		diffQ.row(l-m) += sqrt((2.*l+1.)*(l*l-m*m)/(2.*l-1.))*(invu2 % Q.row(l-m-1));
	}
}

	// test normalisedAssociatedLegendre and normalisedAssociatedLegendreRow
	/*arma::Row<CoordinateType> x(2);
	x(0) = 0.123; x(1) = 0.32;
	normalisedAssociatedLegendre(7, 3, x, Qold, diffQold);
	std::cout << Qold << std::endl;
	normalisedAssociatedLegendreRow(7-1, 3, x, Qold);
	std::cout << Qold << std::endl;*/
	//std::cout << diffQ << std::endl;


template <typename ResultType>
InterpolateOnSphere<ResultType>::InterpolateOnSphere(
	unsigned int Lold, unsigned int Lnew)
	: m_Lold(Lold), m_Lnew(Lnew), m_BT(std::min(Lold, Lnew)+1)
{
	if (m_Lold == m_Lnew || m_Lnew == 0) {
		return;
	}

	// can use precalculated Gauss-Legendre nodes, but just calculate for now
	CoordinateType xoldMem[Lold+1], woldMem[Lold+1];
	LegendreRoots<CoordinateType>(Lold+1, xoldMem, woldMem);
	const arma::Row<CoordinateType> xold(xoldMem, Lold+1);
	const arma::Row<CoordinateType> wold(woldMem, Lold+1);

	CoordinateType xnewMem[Lnew+1], wnewMem[Lnew+1];
	LegendreRoots<CoordinateType>(Lnew+1, xnewMem, wnewMem);
	const arma::Row<CoordinateType> xnew(xnewMem, Lnew+1);

	// As a starting point we have implemented the "semi-naive" implementation of 
	// Alpert for interpolation. Could be improved by 1D FFT (see Darve 2000).
	// N.B. weight function sometimes omitted from presentation, but it is essential
	// (omitted so that weight funtion can be used in up and down pass)
	// B. Alpert and R. Jakob-Chien, A fast spherical filter with uniform resolution, 
	// J. Comput. Phys. 136, 580 (1997).
	// E. Darve. The fast multipole method: Numerical implementation. J. Comp. Phys., 
	// 160:195–240, 2000.

	// instead of calculating the \sum_l=|m|^L^(l+1) Q_l^m(x)Q_l^m(x') directly,
	// use method of Christofell-Darboux (1.60) via Sylvand thesis

	// create B matrices that map from one set of Gauss-Legendre nodes to
	// another for a given order m. B^m = B^{-m} so only need to store the positive
	// orders. Potentially the symmetry of the Gauss nodes about x=0 can be expoited
	// to reduce further storage costs. If stored fully, it is a block diagonal matrix.
	// Store each block as an element of a vector
	//std::vector<arma::Mat<CoordinateType> > B(std::min(Lold, Lnew)+1);

	arma::Mat<CoordinateType> xdist = 
		 arma::ones<arma::Col<CoordinateType> >(Lold+1)*xnew
		-xold.st()*arma::ones<arma::Row<CoordinateType> >(Lnew+1);
	// When Lold and Lnew are both even, xold and xnew have coincident nodes at 0.
	// This must be treated as a special case (see below). N.B. that coincident
	// nodes only seem to appear at x=0 (although should probably check if they get 
	// close elsewhere)
	if ( (Lold&1)==0 && (Lnew&1)==0 ) { // avoid division by zero below
		xdist(Lold/2, Lnew/2) = 1;
	}

	arma::Mat<CoordinateType> xdistInvWold;
	xdistInvWold = (wold.st()*arma::ones<arma::Row<CoordinateType> >(Lnew+1))/xdist;

	for (unsigned int m = 0; m <= std::min(Lold, Lnew); m++) {
		arma::Mat<CoordinateType> Qold, Qnew; // take outside loop?
		normalisedAssociatedLegendreThreeDegrees(Lold, m, xold, Qold);
		normalisedAssociatedLegendreThreeDegrees(Lold, m, xnew, Qnew); // N.B. Lold

		CoordinateType epsilon = sqrt( ((Lold+1)*(Lold+1)-m*m)
			/ (4*(Lold+1)*(Lold+1)-1.) );

		m_BT[m] = epsilon*xdistInvWold % (Qold.row(1).st()*Qnew.row(2)
			 - Qold.row(2).st()*Qnew.row(1) );

		// Apply L'Hopital's rule to treat coincident nodes at xnew = xold = 0
		if ( (Lold&1)==0 && (Lnew&1)==0 ) {
			arma::Col<CoordinateType> diffQ0(2); //dQ_n^m(x)/dx|_{x=0} for n={l, l+1}
			// dQ_l^m(x)/dx|_{x=0} = sqrt((2*l+1)*(l^2-m^2)/(2*l-1))*Q_{l-1}^m(0)
			for (unsigned int l = Lold; l <= Lold+1; l++) {
				diffQ0(l-Lold) = sqrt((2*l+1.)*(l*l-m*m)/(2*l-1))*Qold(l-Lold, Lold/2);	
			}
			m_BT[m](Lold/2, Lnew/2) = epsilon*wold(Lold/2)*( 
				 Qold(1, Lold/2)*diffQ0(2-1) - Qold(2, Lold/2)*diffQ0(1-1));
			// alternatively set Qnew(1:2, Lnew/2) = diffQ0(0:1)
		}
	//	std::cout << m_B[m] << std::endl;
	} // for each m

}


template <typename ResultType>
void 
InterpolateOnSphere<ResultType>::interpolate(
	const arma::Col<ResultType>& coefsOld,
	arma::Col<ResultType>& coefsNew) const
{
	if (m_Lold == m_Lnew) {
		coefsNew = coefsOld;
		return;
	} else if (m_Lnew == 0) {
		coefsNew = coefsOld;
		return;
	}

	// test with an arbitrary function
/*	arma::Mat<ResultType> fold(2*m_Lold+1, m_Lold+1);
	CoordinateType xold[m_Lold+1], wold[m_Lold+1];
	legendre_roots(m_Lold+1, xold, wold);;
	for (unsigned int l=0; l<m_Lold+1; l++) {
		for (unsigned int m=0; m<=2*m_Lold; m++) {
			CoordinateType phi = 2*M_PI*m/(2*m_Lold+1);
			fold(m, l) = (cos(M_PI*xold[l])+4+xold[l])*(1+sin(2*phi)*0.1);
		}
	}
*/
	const arma::Mat<ResultType> fold(coefsOld.memptr(), 2*m_Lold+1, m_Lold+1);

	// take 1D Fourier transforms along phi of the multipole/local coefficients
	// N.B. the coefficients are ordered along phi, since Armadillo fft is columnwise
	// must take care, because potentially multipole coefficients can be real, for
	// instance in the black-box FMM case. However, fft result will always be complex
	arma::Mat<std::complex<CoordinateType> > foldtilde = 1/(2.*m_Lold+1)*arma::fft(fold);

	arma::Mat<std::complex<CoordinateType> > fnewtilde(2*m_Lnew+1, m_Lnew+1);
	fnewtilde.fill(0.);

	fnewtilde.row(0) = foldtilde.row(0)*m_BT[0];
	for (unsigned int m = 1; m <= std::min(m_Lold, m_Lnew); m++) {
		fnewtilde.row(m) = foldtilde.row(m)*m_BT[m];
		fnewtilde.row(2*m_Lnew+1-m) = foldtilde.row(2*m_Lold+1-m)*m_BT[m];
	}
	arma::Mat<std::complex<CoordinateType> > fnew = (2*m_Lnew+1)*arma::ifft(fnewtilde);

	// convert to ResultType and flatten. TODO: avoid copying, but how to avoid type 
	// conversion?
	coefsNew = vectorise( arma::conv_to<arma::Mat<ResultType> >::from(fnew) );

	//std::cout << fnew << std::endl;
	//exit(1);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(InterpolateOnSphere);

} // namespace Bempp
