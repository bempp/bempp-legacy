#define BOOST_TEST_MODULE Element

#include "element.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/*
 * Check that a triangular element's Jacobian is constant everywhere
 */
BOOST_AUTO_TEST_CASE(TriangularElement_jacobian_is_const)
{
  TriangularElement e;
  e.vertices[0].x = 0.;
  e.vertices[0].y = 0.;
  e.vertices[1].x = 1.;
  e.vertices[1].y = 0.;
  e.vertices[2].x = 0.;
  e.vertices[2].y = 1.;

  const int n_elements = 2;
  Point pts[n_elements] = {{0., 1.}, {0.2, -0.6}};
  double result[n_elements];
  e.jacobian(n_elements, pts, result);

  BOOST_CHECK_CLOSE(result[0], result[1], 1e-9);
}

/*
 * Check that a parallelogramic element's Jacobian is constant everywhere
 */
BOOST_AUTO_TEST_CASE(QuadrilateralElement_jacobian_is_const_for_parallelograms)
{
  QuadrilateralElement e;
  e.vertices[0].x = 0.;
  e.vertices[0].y = 0.;
  e.vertices[1].x = 1.;
  e.vertices[1].y = 0.;
  e.vertices[2].x = 2.;
  e.vertices[2].y = 1.;
  e.vertices[3].x = 1.;
  e.vertices[3].y = 1.;

  const int n_elements = 2;
  Point pts[n_elements] = {{0., 1.}, {0.2, -0.6}};
  double result[n_elements];
  e.jacobian(n_elements, pts, result);

  BOOST_CHECK_CLOSE(result[0], result[1], 1e-9);
}

/*
 * Check that a quadrilateral element's Jacobian is constant everywhere ->
 * THIS TEST WILL FAIL
 */
BOOST_AUTO_TEST_CASE(QuadrilateralElement_jacobian_is_const)
{
  QuadrilateralElement e;
  e.vertices[0].x = 0.;
  e.vertices[0].y = 0.;
  e.vertices[1].x = 1.;
  e.vertices[1].y = 0.;
  e.vertices[2].x = 2.;
  e.vertices[2].y = 1.;
  e.vertices[3].x = 0.5;
  e.vertices[3].y = 1.;

  const int n_elements = 2;
  Point pts[n_elements] = {{0., 1.}, {0.2, -0.6}};
  double result[n_elements];
  e.jacobian(n_elements, pts, result);

  BOOST_CHECK_CLOSE(result[0], result[1], 1e-9);
}
