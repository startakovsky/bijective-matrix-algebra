r"""
ReductionMapsDicts

Reduction Mappings are defined in and inspired by the following article in Linear Algebra and its Applications:

Bijective Matrix Algebra: http://www.sciencedirect.com/science/article/pii/S0024379506000218

Essentially the inputs are two sets and two functions which obey certain properties to show that 
one set '''reduces''' another.

AUTHORS:

- Steven Tartakovsky (2012): initial version
"""
#*****************************************************************************
#       Copyright (C) 2012 Steven Tartakovsky <startakovsky@gmail.com>, 
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps
from sage.functions.other import sqrt
from sage.matrix.all import *
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from copy import copy



class ReductionMapsDict(dict):
    """
    INPUT:
     - d a square dictionary with ordered pairs as keys and ReductionMaps as values
    
        
    EXAMPLES::

    TBD

    """
    def __init__(self,dic,repr=None):
        for key in dic.keys():
            self[key] = copy(dic[key])
        self._dim = int(sqrt(len(self.keys())))
        self._repr = repr
        if self._repr is None:
            self._repr = "This is a matrix reduction object: description missing"
        elif type(self._repr)!=str:
            raise ValueError, "Fifth input must be descriptive string or empty argument."
        else:
            self._repr = "This is a matrix reduction object: " + str(self._repr)

    def __repr__(self):
        return self._repr
        
    def get_matrix_A(self):
        dim = self.get_dim()
        mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
        L = list()
        for i in range(dim):
            L.append(list())
            for j in range(dim):
                L[i].append(self[i,j].get_A())
        return mat_space(L)

    def get_matrix_B(self):
        dim = self.get_dim()
        mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
        L = list()
        for i in range(dim):
            L.append(list())
            for j in range(dim):
                L[i].append(self[i,j].get_B())
        return mat_space(L)
              
    def get_reduction_dict(self):
        return self._reduction_dic
    
    def get_dim(self):
        return self._dim
        
    def transitive(self,other):
        """
        Implement transitivity lemma for matrices
        """
        dim = self.get_dim()
        d = dict()
        for i in range(dim):
            for j in range(dim):
                d[i,j] = self.get_reduction_dict()[i,j].transitive(other.get_reduction_dict()[i,j])
        return ReductionMapsDict(d)
    
    def confluence(self,other):
        """
        Implement confluence lemma for matrices
        """
        dim = self.get_dim()
        d = dict()
        for i in range(dim):
            for j in range(dim):
                d[i,j] = self.get_reduction_dict()[i,j].confluence(other.get_reduction_dict()[i,j])
        return ReductionMapsDict(d)
     
    def print_dict(self):
        for key in sorted(self):
            print "row: " + str(key[0]) + ", column: " + str(key[1])
            print self[key]
            func = self[key].get_SRWP()
            for i in func.domain():
                print str(i) + ", " + str(i.get_sign()) + " --> " + str(func(i)) + ", " + str(func(i).get_sign())
            print "***********************************"      