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

#include "element.hpp"

#include <cstdlib>
#include <iostream>
#include <sys/time.h>

using std::endl;
using std::cout;

const int N_POINTS = 1;
const int N_ELEMENTS = 1000;
const int N_TRIALS = 100000;
const int SEED = 5;

// Returns the time needed for calculating the Jacobian 
// N_TRIALS * N_ELEMENTS times at N_POINTS
template <class Element, int NVertices>
double test_element()
{
  srand(SEED);

  Element* elements = new Element[N_ELEMENTS];
  for (int i = 0; i < N_ELEMENTS; ++i)
    for (int v = 0; v < NVertices; ++v)
    {
      elements[i].vertices[v].x = rand() / double(RAND_MAX);
      elements[i].vertices[v].y = rand() / double(RAND_MAX);
    }

  double j[N_POINTS];
  Point pts[N_POINTS];

  j[0] = 0;

  // Some dummy values
  for (int i = 0; i < N_POINTS; ++i)
  {
    pts[i].x = -1. + 0.2 * i;
    pts[i].y = 0.8 - 0.3 * i;
  }

  timeval start, end;
  gettimeofday(&start, 0);
  for (int t = 0; t < N_TRIALS; ++t)
    for (int i = 0; i < N_ELEMENTS; ++i)
      elements[i].jacobian(N_POINTS, pts, j);
  gettimeofday(&end, 0);

  delete[] elements;
  return (end.tv_sec + end.tv_usec / 1000000.) - 
    (start.tv_sec + start.tv_usec / 1000000.);
}

int main()
{
  cout << "Time for TriangularElement: " 
       << test_element<TriangularElement, 3>() << endl;
  cout << "Time for QuadrilateralElement: " 
       << test_element<QuadrilateralElement, 4>() << endl;

  return 0;
}
