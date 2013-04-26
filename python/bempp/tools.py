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
                                evalOps=None)

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

def evaluatePotentialInBox(potentials, gridfuns, coefficients, limits, dimensions, evalOps=None):
    """evaluatePotentialOnPlane(potential, gridfun, limits, dimensions, 
                                evalOps=None)

       Evaluate a linear combination of given potential operators in a box

       Parameters:
       -----------
       potentials   : List of instances of a potential operators.
       gridfuns     : List grid functions with which to evaluate the potentials
       coefficients : Scalar coefficients of the potentials
       limits       : tuple (xmin,xmax,ymin,ymax,zmin,zmax) specifying the
                      extent of the plane.
       dimensions   : tuple (xdim,ydim,zdim) specifying the number of points in
                      the (x,y) dimensions of the plane.
       evalOps      : Optional EvaluationOptions object. Use default options
                      if none is given.

       Returns:
       --------
       points     : A numpy.mgrid object defining the points in 3d space.
       values     : Array [v_1,v_2,...], where each v_i has the same shape as the
                    components of the mgrid object points. Each element v_i corresponds
                    to one dimension of the return value of the potential, e.g. for Maxwell
                    the values v_1,v_2, and v_3 correspond to the x, y and z component of the
                    potential field.
    """


    pointsArray, pointsMgrid = pointsInBox(limits,dimensions,mgridObj=True)

    if evalOps is None:
        evalOps = lib.createEvaluationOptions()

    values = coefficients[0]*potentials[0].evaluateAtPoints(gridfuns[0],pointsArray,evalOps)
    for i in range(1,len(potentials)):
        values +=coefficients[i]*potentials[i].evaluateAtPoints(gridfuns[i],pointsArray,evalOps)
    valsMgrid = np.array([v.reshape(dimensions[::-1][:]) for v in values])

    return (pointsMgrid, valsMgrid)


def pointsInBox(limits, dimensions, mgridObj=False):
    """pointsInBox(limits, dimensions, mgridObj=False) 

       Define a set of points in a three dimensional box

       Parameters:
       -----------
       limits       : tuple (xmin,xmax,ymin,ymax,zmin,zmax) specifying the
                      extent of the plane.
       dimensions   : tuple (xdim,ydim,zdim) specifying the number of points in
                      the (x,y) dimensions of the plane.
       mgridObj     : If True return the points also as np.mgrid object

       Returns:
       --------
       pointsArray : An array [p0,p1,...] defining the points in 3d space 
       pointsMgrid : A numpy.mgrid object defining the points in 3d space
                    (only if mgridObj=True)
    """

    xmin, xmax, ymin, ymax, zmin, zmax = limits
    pointsMgrid = np.mgrid[xmin:xmax:dimensions[0]*1j,
                    ymin:ymax:dimensions[1]*1j,
                    zmin:zmax:dimensions[2]*1j]
    x, y, z = pointsMgrid

    pointsArray = np.array([x.ravel(),y.ravel(), z.ravel()],dtype ='d')

    if mgridObj == True:
        return (pointsArray,pointsMgrid)
    else:
        return pointsArray


def array2Mgrid(values,dimensions):
    """array2Mgrid(values,dimensions)

       Reshape an array of scalar values to a shape compatible with coordinates
       given as mgrid object.

       Parameters:
       -----------
       values      : Array [v_1,v_2,...] of scalar values
       dimensions  : Dimension (xdim,ydim,zdim) for reshaping the input array

       Returns:
       --------
       An [v_0,v_1,v_2], where each v_i corresponds to one component of the input
       values array converted into a shape compatible with coordinates given as mgrid
       object.
     """
        
    return np.array([v.reshape(dimensions[::-1][:]) for v in values])

def mgridPoints2Array(points):
    """mgridPoints2Mgrid(values)

       Reshape an mgrid array of points into an array of points that is
       compatible with input needed by BEM++ routines

       Parameters:
       -----------
       points     : Array of points coming from numpy.mgrid

       Returns:
       --------
       An array of points [p_0,p_1,...].
     """
        
    return np.array([p.ravel() for p in points])



class RealOperator(object):
    """RealOperator(operator)

       Convert a complex operator into a real operator of twice the dimension

       Parameters:
       -----------
       operator    : A Python object that supports that provides a 'matvec'
                     method and a 'dtype' attribute.

    """



    def __init__(self,operator):
        
        self.__operator = operator
        self.__operator_shape = operator.shape
        self.shape = (2*self.__operator_shape[0],2*self.__operator_shape[1])
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

    @property
    def operator(self):
        return self.__operator

    def matvec(self,x):

        if len(x.shape) == 1:
            if x.shape[0] != self.shape[1]:
                raise Exception("RealOperator.matvec: wrong dimension.")
        elif len(x.shape) == 2:
            if x.shape[1] != 1 or x.shape[0] != self.shape[1]:
                raise Exception("RealOperator.matvec: wrong dimension.")

        res = self.__operator.matvec(x[:self.__operator_shape[1]]+1j*x[self.__operator_shape[1]:])
        
        if len(res.shape)==1:
            return np.hstack([np.real(res),np.imag(res)])
        else:
            return np.vstack([np.real(res),np.imag(res)])

    def matmat(self,x):

        if x.shape[0] != self.shape[1]:
            raise Exception("RealOperator.matmat: wrong dimension")

        res = self.__operator.matmat(x[:self.__operator_shape[1]]+1j*x[self.__operator_shape[1]:])
        return np.vstack([np.real(res),np.imag(res)])

            
        
