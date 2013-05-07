/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with
    this software; if not, see AHMED's internet site.
*/

/*
    This file has been created by the authors of BEM++.

    Summary of changes:

    - 7/05/2013: initial version. The bbxbemcluster class template provides the
      geticom() routine needed by functions from apprx.h, but unlike to
      bemcluster it is derived from cluster_bbx rather than cluster3d_pca.
*/

#ifndef BEMCLUSTER_H
#define BEMCLUSTER_H

#include "cluster_bbx.h"

template<class T>
class bbxbemcluster : public cluster_bbx<T>
{
private:
    typedef cluster_bbx<T> Base;

public:
    // cluster with entries k <= i < l
    bbxbemcluster(T* dofs, unsigned k, unsigned l) :
        Base(dofs, k, l), icom(0) { }

    //! returns the index of the degree of freedom closest to the centroid
    unsigned geticom() const { return icom; }

    //! set the index of the degree of freedom closest to the centroid
    void seticom(unsigned icom_) { icom = icom_; }

private:
    // index closest to the centroid (needed for the starting index of ACA)
    unsigned icom;
};

#endif
