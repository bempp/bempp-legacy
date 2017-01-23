// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef cuda_const_memory_hpp
#define cuda_const_memory_hpp

// Constant memory is typically 64 KB
// Use double as most memory consuming possible data type to ensure enough space
extern __constant__ double constTestQuadWeights[6];
extern __constant__ double constTrialQuadWeights[6];

extern __constant__ double constTestGeomShapeFun0[6];
extern __constant__ double constTestGeomShapeFun1[6];
extern __constant__ double constTestGeomShapeFun2[6];

extern __constant__ double constTrialGeomShapeFun0[6];
extern __constant__ double constTrialGeomShapeFun1[6];
extern __constant__ double constTrialGeomShapeFun2[6];

#endif /* cuda_const_memory_hpp */
