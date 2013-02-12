# Copyright (C) 2011-2012 by the BEM++ Authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
from bempp import lib

def evaluatePotentialOnPlane(potential, gridfun, limits, dimensions, plane="xy",
                             origin=[0,0,0], evalOps=None):
    """evaluatePotentialOnPlane(potential, gridfun, limits, dimensions, plane="xy",
                                evalOps=None, quadStrategy=None)

       Evaluate a given potential operator on a plane

       Parameters:
       -----------
       potential    : Instance of a potential operator.
       gridfun      : grid function with which to evaluate the potential
       limits       : tuple (xmin,xmax,ymin,ymax) specifying the
                      extent of the plane.
       dimensions   : tuple (xdim,ydim) specifying the number of points in
                      the (x,y) dimensions of the plane.
       plane        : either "xy", "xz" or "yz", (default: "xy"); the orientation
                      of the plane in 3d space.
       origin       : origin of the plane (default: [0,0,0])
       evalOps      : Optional EvaluationOptions object. Use default options
                      if none is given.

       Returns:
       --------
       points     : Array [p_1,p_2,,,] of points in 3d space.
       values     : Array [v_1,v_2,..] of values v_i of the potential at points p_i.
    """

    if not plane in ["xy","xz","yz"]:
        raise ValueError("'plane' must by either 'xy', 'xz', or 'yz'")

    xmin, xmax, ymin, ymax = limits
    x, y = np.mgrid[xmin:xmax:dimensions[0]*1j,
                    ymin:ymax:dimensions[1]*1j]

    points_2d = np.array([x.T.ravel(),y.T.ravel()],dtype ='d')

    dims = dimensions

    if plane=="xy":
        points = np.array([points_2d[0]+origin[0],points_2d[1]+origin[1],np.zeros(dims[0]*dims[1],dtype='d')+origin[2]])
    elif plane=="xz":
        points = np.array([points_2d[0]+origin[0],np.zeros(dims[0]*dims[1],dtype='d')+origin[1],points_2d[1]+origin[2]])
    elif plane=="yz":
        points = np.array([np.zeros(dims[0]*dims[1],dtype='d')+origin[0],points_2d[0]+origin[1],points_2d[1]+origin[2]])

    if evalOps is None:
        evalOps = lib.createEvaluationOptions()
    values = potential.evaluateAtPoints(gridfun, points, evalOps)

    return (points, values)

class RealOperator(object):
    """RealOperator(operator)

       Convert a complex operator into a real operator of twice the dimension

       Parameters:
       -----------
       operator    : A Python object that supports that provides a 'matvec'
                     method and a 'dtype' attribute.

    """



    def __init__(self,operator):
        
        self.operator = operator
        self.operator_shape = operator.shape
        self.shape = (2*self.operator_shape[0],2*self.operator_shape[1])
        if operator.dtype==np.dtype('complex128'):
            self.dtype = np.dtype('float64')
        elif operator.dtype==np.dtype('complex64'):
            self.dtype = np.dtype('float32')
        elif operator.dtype==np.dtype('float64'):
            self.dtype = np.dtype('float64')
        elif operator.dtype==np.dtype('float32'):
            self.dtype = np.dtype('float32')
        else:
            raise Exception('RealOperator.__init__: Datatype of operator not supported.')


    def matvec(self,x):

        if len(x.shape) == 1:
            if x.shape[0] != self.shape[1]:
                raise Exception("RealOperator.matvec: wrong dimension.")
        elif len(x.shape) == 2:
            if x.shape[1] != 1 or x.shape[0] != self.shape[1]:
                raise Exception("RealOperator.matvec: wrong dimension.")

        res = self.operator.matvec(x[:self.operator_shape[1]]+1j*x[self.operator_shape[1]:])
        
        if len(res.shape)==1:
            return np.hstack([np.real(res),np.imag(res)])
        else:
            return np.vstack([np.real(res),np.imag(res)])




            
        
