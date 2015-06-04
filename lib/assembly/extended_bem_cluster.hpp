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

#ifndef bempp_extended_bem_cluster_hpp
#define bempp_extended_bem_cluster_hpp

#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "../common/common.hpp"
#include "ahmed_aux_fwd.hpp"

#include <bbxbemcluster.h>

namespace Bempp {

template <typename T> class ExtendedBemCluster : public bbxbemcluster<T> {
private:
  typedef bbxbemcluster<T> Base;
  enum { Dim = 3 }; // in future we might make this configurable

  // xminmax (inherited from cluster_bbx) is the bounding box of the
  // *reference points* of DOFs. In contrast, extMinmax is the bounding box
  // of the *bounding boxes* of DOFs. The former is used in cluster splitting;
  // the latter to determine the distance between clusters and check the
  // admissibility condition

protected:
  double extMinmax[2 * Dim];

public:
  // cluster with entries k <= i < l
  ExtendedBemCluster(
      T *dofs, unsigned *op_perm, unsigned k, unsigned l,
      unsigned int maximumBlockSize = std::numeric_limits<unsigned int>::max(),
      bool strongAdmissibility = false)
      : Base(dofs, k, l), m_maximumBlockSize(maximumBlockSize),
        m_strongAdmissibility(strongAdmissibility) {

    this->xminmax = new double[2 * Dim];

    for (int i = 0; i < Dim; ++i)
      this->xminmax[i] = std::numeric_limits<double>::max();
    for (int i = 0; i < Dim; ++i)
      this->xminmax[Dim + i] = -std::numeric_limits<double>::max();

    for (int i = 0; i < Dim; ++i)
      this->extMinmax[i] = std::numeric_limits<double>::max();
    for (int i = 0; i < Dim; ++i)
      this->extMinmax[Dim + i] = -std::numeric_limits<double>::max();

    // Compute the bounding box
    for (unsigned j = this->nbeg; j < this->nend; ++j) {
      const T *v = this->dofs + op_perm[j];
      for (int i = 0; i < Dim; ++i) {
        this->xminmax[i] = std::min<double>(this->xminmax[i], v->getcenter(i));
        this->xminmax[Dim + i] =
            std::max<double>(this->xminmax[Dim + i], v->getcenter(i));
        this->extMinmax[i] =
            std::min<double>(this->extMinmax[i], v->getlbound(i));
        this->extMinmax[Dim + i] =
            std::max<double>(this->extMinmax[Dim + i], v->getubound(i));
      }
    }

    // Calculate diam2 and main direction
    this->diam2 = 0.;
    this->maindir = 0;
    double maxExtent = 0.;
    for (int i = 0; i < Dim; ++i) {
      const double e = this->xminmax[Dim + i] - this->xminmax[i];
      this->diam2 += e * e;
      if (e > maxExtent) {
        maxExtent = e;
        this->maindir = i;
      }
    }

    // Calculate cntrdir
    this->cntrdir = 0.5 * (this->xminmax[this->maindir] +
                           this->xminmax[Dim + this->maindir]);

    // Compute the index of the dof clostest to centre of mass.
    // Code copied from (patched) bemcluster.h.

    this->seticom(
        op_perm[this->nbeg]); // original index. Will be replaced by permuted
                              // index in createClusterTree()

    double smin = 0.;
    for (int j = 0; j < Dim; ++j) {
      const double e = this->getcom(j) - dofs[op_perm[this->nbeg]].getcenter(j);
      smin += e * e;
    }

    for (int i = this->nbeg + 1; i < this->nend; ++i) {
      double s = 0.;
      for (int j = 0; j < Dim; ++j) {
        const double e = this->getcom(j) - dofs[op_perm[i]].getcenter(j);
        s += e * e;
      }
      if (s < smin) {
        this->seticom(op_perm[i]);
        smin = s;
      }
    }
  }

  virtual ExtendedBemCluster *clone(unsigned int *op_perm, unsigned int beg,
                                    unsigned int end) const {
    return new ExtendedBemCluster(this->dofs, op_perm, beg, end,
                                  m_maximumBlockSize, m_strongAdmissibility);
  }

  virtual unsigned dim() const { return Dim; }

  virtual bool isadm(double eta2, cluster *cl, bl_info &info) {
    if (this->size() > m_maximumBlockSize || cl->size() > m_maximumBlockSize)
      return (info.is_adm = false);
    else {
      ExtendedBemCluster<T> *p = dynamic_cast<ExtendedBemCluster<T> *>(cl);
      assert(p != NULL);
      info.is_sep = false;

      const double d2 = std::min(this->diam2, p->diam2);

      bool strongAdmissibility =
          m_strongAdmissibility || p->usingStrongAdmissibilityCondition();
      if (d2 < eta2 * centerToCenterDist2(p) &&
          (!strongAdmissibility || 0 < this->extDist2(p)))
        return (info.is_adm = true);
      else
        return (info.is_adm = false);
    }
  }

  // once all descendants are created (so that the array po_perm[nbeg...nend)
  // won't change any more), set icom to the permuted index of the centroid
  virtual void createClusterTree(const unsigned bmin, unsigned *op_perm,
                                 unsigned *po_perm) {
    Base::createClusterTree(bmin, op_perm, po_perm);
    this->seticom(po_perm[this->geticom()]);
  }

  double getcom(unsigned i) const {
    return (this->xminmax[i + 3] + this->xminmax[i]) / 2.;
  }

  //! returns square of the distance between the "exterior" bounding boxes of
  //! this cluster and \p cl (dependent on extMinmax rather than xminmax)
  double extDist2(const ExtendedBemCluster *cl) const {
    double d = 0.;
    const unsigned n = dim();

    for (unsigned i = 0; i < n; ++i) {
      const double e = cl->extMinmax[i] - extMinmax[i + n];
      if (e > 0.)
        d += e * e;
      else {
        const double e1 = extMinmax[i] - cl->extMinmax[i + n];
        if (e1 > 0.)
          d += e1 * e1;
      }
    }
    return d;
  }

  //! returns square of the distance between centres of this cluster and \p cl
  double centerToCenterDist2(const ExtendedBemCluster *cl) {
    double d2 = 0.0;
    const unsigned n = cl->dim();

    for (unsigned j = 0; j < n; ++j) {
      const double d = this->getcom(j) - cl->getcom(j);
      d2 += d * d;
    }
    return d2;
  }

  //! debug routine: print bounding boxes of all dofs contained in this cluster
  void printBboxes() const {
    std::cout << "Cluster with " << this->nend - this->nbeg << " dofs;\nbbox: ";
    for (int i = 0; i < Dim; ++i)
      std::cout << this->xminmax[i] << " -- " << this->xminmax[i + Dim] << "; ";
    std::cout << "\nexterior bbox: ";
    for (int i = 0; i < Dim; ++i)
      std::cout << this->extMinmax[i] << " -- " << this->extMinmax[i + Dim]
                << "; ";
    std::cout << std::endl;
  }

  unsigned int maximumBlockSize() const { return m_maximumBlockSize; }

  bool usingStrongAdmissibilityCondition() const {
    return m_strongAdmissibility;
  }

  void useStrongAdmissibilityCondition(bool value = true) {
    m_strongAdmissibility = value;
    for (int i = 0; i < this->getns(); ++i) {
      cluster *son = this->getson(i);
      if (ExtendedBemCluster *exbemson =
              dynamic_cast<ExtendedBemCluster *>(son))
        exbemson->useStrongAdmissibilityCondition(value);
    }
  }

  void clearDofPointers() {
    this->dofs = 0;
    for (int i = 0; i < this->getns(); ++i) {
      cluster *son = this->getson(i);
      if (ExtendedBemCluster *exbemson =
              dynamic_cast<ExtendedBemCluster *>(son))
        exbemson->clearDofPointers();
    }
  }

private:
  unsigned int m_maximumBlockSize;
  bool m_strongAdmissibility;
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
