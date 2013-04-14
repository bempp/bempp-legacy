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
from PyTrilinos import Epetra, AztecOO, Anasazi, Teuchos

Vector = Epetra.Vector


class PyTrilinosOperator(Epetra.Operator):

    def __init__(self,operator,pid=0,label=""):
        Epetra.Operator.__init__(self)


        from tools import RealOperator    

        self.__comm = Epetra.PyComm()
        self.__label = label
        self.__pid = pid

        if self.__comm.MyPID() == pid:
            if operator.dtype == 'complex128' or operator.dtype == 'complex64':
                self.__operator = RealOperator(operator)
            else:
                self.__operator = operator
        else:
            self.__operator = None




        if self.__comm.MyPID() == pid:
            num_my_elems_domain = self.__operator.shape[1]
            num_my_elems_range = self.__operator.shape[0]
        else:
            num_my_elems_domain = 0
            num_my_elems_range = 0


        self.__rangeMap = Epetra.Map(-1,num_my_elems_range,0,self.__comm)
        self.__domainMap = Epetra.Map(-1,num_my_elems_domain,0,self.__comm)


    def Label(self):
        return self.__label

    def OperatorDomainMap(self):
        return self.__domainMap

    def OperatorRangeMap(self):
        return self.__rangeMap

    def Comm(self):
        return self.__comm

    def ApplyInverse(self,x,y):
        return -1

    def HasNormInf(self):
        return False

    def NormInf(self):
        return -1

    def SetUseTranspose(self, useTranspose):
        return -1

    def UseTranspose(self):
        return 0

    def Apply(self,x,y):
        try:
            if self.__comm.MyPID() == self.__pid:
                xvec = x.ExtractView().T
                yvec = y.ExtractView().T
                if len(xvec.shape) == 1:
                    yvec[:] = self.__operator.matvec(xvec)
                else:
                    yvec[:] = self.__operator.matmat(xvec)
            return 0
        except Exception, e:
            print "Exceptin in "+self.__label+".apply:"
            print e
            return -1


class PyTrilinosInverseOperator(Epetra.Operator):

    def __init__(self,operator,pid=0,label=""):
        Epetra.Operator.__init__(self)


        from tools import RealOperator    

        self.__comm = Epetra.PyComm()
        self.__label = label
        self.__pid = pid

        if self.__comm.MyPID() == pid:
            if operator.dtype == 'complex128' or operator.dtype == 'complex64':
                self.__operator = RealOperator(operator)
            else:
                self.__operator = operator
        else:
            self.__operator = None




        if self.__comm.MyPID() == pid:
            num_my_elems_domain = self.__operator.shape[1]
            num_my_elems_range = self.__operator.shape[0]
        else:
            num_my_elems_domain = 0
            num_my_elems_range = 0


        self.__rangeMap = Epetra.Map(-1,num_my_elems_range,0,self.__comm)
        self.__domainMap = Epetra.Map(-1,num_my_elems_domain,0,self.__comm)


    def Label(self):
        return self.__label

    def OperatorDomainMap(self):
        return self.__domainMap

    def OperatorRangeMap(self):
        return self.__rangeMap

    def Comm(self):
        return self.__comm

    def Apply(self,x,y):
        return -1

    def HasNormInf(self):
        return False

    def NormInf(self):
        return -1

    def SetUseTranspose(self, useTranspose):
        return -1

    def UseTranspose(self):
        return 0

    def ApplyInverse(self,x,y):
        try:
            if self.__comm.MyPID() == self.__pid:
                xvec = x.ExtractView().T
                yvec = y.ExtractView().T
                if len(xvec.shape) == 1:
                    yvec[:] = self.__operator.matvec(xvec)
                else:
                    yvec[:] = self.__operator.matmat(xvec)
            return 0
        except Exception, e:
            print "Exceptin in "+self.__label+".apply:"
            print e
            return -1







       

        
