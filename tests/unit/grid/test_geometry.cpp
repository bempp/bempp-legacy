// Copyright (C) 2011 by the BEM++ Authors
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

#include "test_entity.hpp"
#include "grid/armadillo_helpers.hpp"
#include "grid/geometry.hpp"

#include <boost/test/unit_test.hpp>
#include "../num_template.hpp"

using namespace Bempp;
typedef double ctype;


BOOST_FIXTURE_TEST_SUITE(Geometry_Triangular, TriangularEntityManager)

// type()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(type_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    BOOST_CHECK_EQUAL(geo.type(), duneGeo.type());
}

// affine()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(affine_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    BOOST_CHECK_EQUAL(geo.affine(), duneGeo.affine());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(affine_is_true_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();

    BOOST_CHECK(geo.affine());
}

// cornerCount()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(cornerCount_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    BOOST_CHECK_EQUAL((int)geo.cornerCount(), duneGeo.corners());
}

// corners()

//_for_single_point

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(corners_first_corner_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    const int nTestedCorner = 0;

    arma::Mat<ctype> corners;
    geo.getCorners(corners);

    Dune::FieldVector<ctype, DuneGrid::dimensionworld> duneCorner = duneGeo.corner(nTestedCorner);

    BOOST_CHECK_EQUAL(corners.col(nTestedCorner), duneCorner);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(corners_last_corner_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    const int nTestedCorner = geo.cornerCount() - 1;

    arma::Mat<ctype> corners;
    geo.getCorners(corners);

    Dune::FieldVector<ctype, DuneGrid::dimensionworld> duneCorner = duneGeo.corner(nTestedCorner);

    BOOST_CHECK_EQUAL(corners.col(nTestedCorner), duneCorner);
}

// ____________________________________________________________________________
// local2global()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(local2global_agrees_with_Dune_for_one_point_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, 1);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int i = 0; i < dimLocal; ++i)
        local(i,0) = duneLocal[i] = 0.1 * (i + 1);

    arma::Mat<ctype> global;
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;

    geo.local2global(local, global);
    duneGlobal = duneGeo.global(duneLocal);

    BOOST_CHECK_EQUAL(global.col(0), duneGlobal);
}

// Helper function for the two following tests
template <int codim>
static void local2global_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    // Note: matrix global is uninitialised
    arma::Mat<ctype> global;
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;

    geo.local2global(local, global);
    duneGlobal = duneGeo.global(duneLocal);

    BOOST_CHECK_EQUAL(global.col(nTestedPoint), duneGlobal);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(local2global_agrees_with_Dune_for_first_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    local2global_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(local2global_agrees_with_Dune_for_last_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    local2global_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

// Helper function for the two following tests
template <int codim>
static void local2global_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    // Note: matrix global is initialised to incorrect shape
    // (local2global is expected to resize it)
    arma::Mat<ctype> global(50, 40);
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;

    geo.local2global(local, global);
    duneGlobal = duneGeo.global(duneLocal);

    BOOST_CHECK_EQUAL(global.col(nTestedPoint), duneGlobal);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(local2global_agrees_with_Dune_for_first_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    local2global_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(local2global_agrees_with_Dune_for_last_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    local2global_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

// ____________________________________________________________________________
// global2local()

// Note: we do not test here whether the calculated values actually make sense
// (they in fact do not since the global points don't necessarily lie on the
// chosen entity). Such tests properly belong to the FoamGrid's test suite. We
// mostly test the correctness of the mapping from Dune field vectors to
// Armadillo matrices.

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(global2local_agrees_with_Dune_for_one_point_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> global(dimGlobal, 1);
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;
    for (int i = 0; i < dimGlobal; ++i)
        global(i,0) = duneGlobal[i] = 0.1 * (i + 1);

    arma::Mat<ctype> local;
    geo.global2local(global, local);
    Dune::FieldVector<ctype, dimLocal> duneLocal = duneGeo.local(duneGlobal);

    BOOST_CHECK_EQUAL(local.col(0), duneLocal);
}

// Helper function for the two following tests
template <int codim>
static void global2local_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> global(dimGlobal, nPoints);
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimGlobal; ++i)
            global(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimGlobal; ++i)
        duneGlobal[i] = global(i,nTestedPoint);

    // Note: matrix local is uninitialised
    arma::Mat<ctype> local;
    geo.global2local(global, local);
    Dune::FieldVector<ctype, dimLocal> duneLocal = duneGeo.local(duneGlobal);

    BOOST_CHECK_EQUAL(local.col(nTestedPoint), duneLocal);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(global2local_agrees_with_Dune_for_first_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    global2local_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(global2local_agrees_with_Dune_for_last_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    global2local_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

// Helper function for the two following tests
template <int codim>
static void global2local_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> global(dimGlobal, nPoints);
    Dune::FieldVector<ctype, dimGlobal> duneGlobal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimGlobal; ++i)
            global(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimGlobal; ++i)
        duneGlobal[i] = global(i,nTestedPoint);

    // Note: matrix local is initialised to incorrect shape
    // (global2local is expected to resize it)
    arma::Mat<ctype> local(50, 40);
    geo.global2local(global, local);
    Dune::FieldVector<ctype, dimLocal> duneLocal = duneGeo.local(duneGlobal);

    BOOST_CHECK_EQUAL(local.col(nTestedPoint), duneLocal);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(global2local_agrees_with_Dune_for_first_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    global2local_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(global2local_agrees_with_Dune_for_last_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    global2local_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

// ____________________________________________________________________________
// integrationElement()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(integrationElement_agrees_with_Dune_for_one_point_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    //const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, 1);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int i = 0; i < dimLocal; ++i)
        local(i,0) = duneLocal[i] = 0.1 * (i + 1);

    arma::Row<ctype> intElement;
    geo.getIntegrationElements(local, intElement);
    ctype duneIntElement = duneGeo.integrationElement(duneLocal);

    BOOST_CHECK_EQUAL(intElement(0), duneIntElement);
}

// Helper function for the two following tests
template <int codim>
static void integrationElement_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    arma::Row<ctype> intElement;
    geo.getIntegrationElements(local, intElement);
    ctype duneIntElement = duneGeo.integrationElement(duneLocal);

    BOOST_CHECK_EQUAL(intElement(nTestedPoint), duneIntElement);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(integrationElement_agrees_with_Dune_for_first_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    integrationElement_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(integrationElement_agrees_with_Dune_for_last_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    integrationElement_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

// Helper function for the two following tests
template <int codim>
static void integrationElement_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    // Note: vector intElement is initialised to incorrect shape
    // (integrationElement is expected to resize it)
    arma::Row<ctype> intElement(50);
    geo.getIntegrationElements(local, intElement);
    ctype duneIntElement = duneGeo.integrationElement(duneLocal);

    BOOST_CHECK_EQUAL(intElement(nTestedPoint), duneIntElement);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(integrationElement_agrees_with_Dune_for_first_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    integrationElement_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(integrationElement_agrees_with_Dune_for_last_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    integrationElement_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

//virtual ctype volume() const {
//    return m_dune_geometry->volume();
//}

//virtual void center(arma::Col<ctype>& c) const {

// ____________________________________________________________________________
// volume()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(volume_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    BOOST_CHECK_EQUAL(geo.volume(), duneGeo.volume());
}

// ____________________________________________________________________________
// center()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(center_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int dimGlobal = DuneGrid::dimensionworld;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    const typename DuneGrid::Codim<codim>::Entity::Geometry& duneGeo = duneEp->geometry();

    arma::Col<ctype> center;
    geo.getCenter(center);
    Dune::FieldVector<ctype, dimGlobal> duneCenter = duneGeo.center();

    BOOST_CHECK_EQUAL(center, duneCenter);
}

// ____________________________________________________________________________
// jacobianTransposed()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianTransposed_agrees_with_Dune_for_one_point_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    //const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    typedef typename DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;
    const DuneGeometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, 1);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int i = 0; i < dimLocal; ++i)
        local(i,0) = duneLocal[i] = 0.1 * (i + 1);

    arma::Cube<ctype> jacobianT;

    geo.getJacobiansTransposed(local, jacobianT);
    Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::mydimension,DuneGeometry::coorddimension> duneJacobianT = duneGeo.jacobianTransposed(duneLocal);

    // work around a "bug" in Armadillo (which has some issues with handling empty cubes)
    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianT.slice(0), duneJacobianT);
    else
        BOOST_CHECK(jacobianT.is_empty());
}

// Helper function for the two following tests
template <int codim>
static void jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    typedef typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;
    const DuneGeometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    arma::Cube<ctype> jacobianT;
    geo.getJacobiansTransposed(local, jacobianT);
    Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::mydimension,DuneGeometry::coorddimension> duneJacobianT
      = duneGeo.jacobianTransposed(duneLocal);

    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianT.slice(nTestedPoint), duneJacobianT);
    else
        BOOST_CHECK(jacobianT.is_empty());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianTransposed_agrees_with_Dune_for_first_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianTransposed_agrees_with_Dune_for_last_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

// Helper function for the two following tests
template <int codim>
static void jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    typedef typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;
    const DuneGeometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    // Note: matrix jacobianT is initialised to incorrect shape
    // (jacobianTransposed is expected to resize it)
    arma::Cube<ctype> jacobianT(10,10,10);
    geo.getJacobiansTransposed(local, jacobianT);
    typedef Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::mydimension,DuneGeometry::coorddimension>
      DuneJacobianTransposedType;
    DuneJacobianTransposedType duneJacobianT = duneGeo.jacobianTransposed(duneLocal);

    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianT.slice(nTestedPoint), duneJacobianT);
    else
        BOOST_CHECK(jacobianT.is_empty());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianTransposed_agrees_with_Dune_for_first_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianTransposed_agrees_with_Dune_for_last_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    jacobianTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

// ____________________________________________________________________________
// jacobianInverseTransposed()


BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianInverseTransposed_agrees_with_Dune_for_one_point_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    //const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
            getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    typedef typename DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;
    const DuneGeometry& duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, 1);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int i = 0; i < dimLocal; ++i)
        local(i,0) = duneLocal[i] = 0.1 * (i + 1);

    arma::Cube<ctype> jacobianInvT;

    geo.getJacobianInversesTransposed(local, jacobianInvT);
    Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::coorddimension,DuneGeometry::mydimension> duneJacobianInvT
            = duneGeo.jacobianInverseTransposed(duneLocal);

    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianInvT.slice(0), duneJacobianInvT);
    else
        BOOST_CHECK(jacobianInvT.is_empty());
}

// Helper function for the two following tests
template <int codim>
static void jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    const Geometry& geo = ep->entity().geometry();
    typedef typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;

    //const typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry&
    DuneGeometry duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    arma::Cube<ctype> jacobianInvT;
    geo.getJacobianInversesTransposed(local, jacobianInvT);
    Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::coorddimension,DuneGeometry::mydimension> duneJacobianInvT
            = duneGeo.jacobianInverseTransposed(duneLocal);

    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianInvT.slice(nTestedPoint), duneJacobianInvT);
    else
        BOOST_CHECK(jacobianInvT.is_empty());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianInverseTransposed_agrees_with_Dune_for_first_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianInverseTransposed_agrees_with_Dune_for_last_of_several_points_and_uninitialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_uninitialised_output_and_codim<codim>(*this, nTestedPoint);
}

// Helper function for the two following tests
template <int codim>
static void jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim(
    BOOST_AUTO_TEST_CASE_FIXTURE& f, int nTestedPoint)
{
    //const int dimGlobal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimensionworld;
    const int dimLocal = BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = f.getPointerToSecondEntityOnLevel0<codim>();
    typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::LevelGridView::Codim<codim>::Iterator duneEp =
        f.getDunePointerToSecondEntityOnLevel0<codim>();

    typedef typename BOOST_AUTO_TEST_CASE_FIXTURE::DuneGrid::Codim<codim>::Entity::Geometry DuneGeometry;
    const Geometry& geo = ep->entity().geometry();
    DuneGeometry duneGeo = duneEp->geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    for (int i = 0; i < dimLocal; ++i)
        duneLocal[i] = local(i,nTestedPoint);

    // Note: matrix jacobianInvT is initialised to incorrect shape
    // (getJacobianInversesTransposed is expected to resize it)
    arma::Cube<ctype> jacobianInvT(10,10,10);
    geo.getJacobianInversesTransposed(local, jacobianInvT);
    Dune::FieldMatrix<typename DuneGeometry::ctype,DuneGeometry::coorddimension,DuneGeometry::mydimension>
            duneJacobianInvT = duneGeo.jacobianInverseTransposed(duneLocal);

    if (dimLocal > 0)
        BOOST_CHECK_EQUAL(jacobianInvT.slice(nTestedPoint), duneJacobianInvT);
    else
        BOOST_CHECK(jacobianInvT.is_empty());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianInverseTransposed_agrees_with_Dune_for_first_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 0;
    jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(jacobianInverseTransposed_agrees_with_Dune_for_last_of_several_points_and_initialised_output_and_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    const int nTestedPoint = 4;
    jacobianInverseTransposed_agrees_with_Dune_for_nth_of_several_points_and_initialised_output_and_codim<codim>(*this, nTestedPoint);
}

BOOST_AUTO_TEST_SUITE_END()

