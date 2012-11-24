r"""
ReductionMaps

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

from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.bijectivematrixalgebra.map_methods import fixed_points
from sage.bijectivematrixalgebra.map_methods import is_SRWP_involution
from sage.bijectivematrixalgebra.map_methods import is_SPWP_bijection
from sage.bijectivematrixalgebra.map_methods import inverse
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.structure.sage_object import SageObject
from sage.sets.finite_set_map_cy import FiniteSetEndoMap_Set
from sage.sets.finite_set_maps import FiniteSetMap_Set
from copy import copy
#from sage.symbolic.expression import Expression



class ReductionMaps(SageObject):
    """
    INPUT:
     - A a Combinatorial Scalar
     - B a Combinatorial Scalar
     - f a SRWP involution on A
     - f0 a SPWP bijection from fixed(f) to B
    
        
    EXAMPLES::

    TBD

    """
    def __init__(self, A,B,f,f0):
        """
        TBD
        """
        
        if type(A)!=CombinatorialScalarWrapper:
            raise ValueError, "The first input must be a Combinatorial Scalar Wrapper"
        elif type(B)!=CombinatorialScalarWrapper:
            raise ValueError, "The second input must be a Combinatorial Scalar Wrapper"
        elif not(bool(A.get_generating_function()==B.get_generating_function())):
            raise ValueError, "The generating functions of the scalars are not equal"
        elif type(f)!= FiniteSetMap_Set and type(f) != FiniteSetEndoMap_Set:
            raise ValueError, "The third input must be a map in FiniteSetMaps"
        elif type(f0) != FiniteSetMap_Set and type(f0) != FiniteSetEndoMap_Set:
            raise ValueError, "The fourth input must be a map in FiniteSetMaps"
        elif set(f.domain()) != set(A):
            raise ValueError, "The third input must have domain of first input"
        elif set(f0.domain()) != set(fixed_points(f)):
            raise ValueError, "The fourth input must have domain of fixed points of third the input"
        elif set(f0.image_set())!=set(B):
            raise ValueError, "The fourth input must have image of second input"
        elif not(is_SRWP_involution(f)):
            raise ValueError, "The third input must be an SRWP involution"
        elif not(is_SPWP_bijection(f0)):
            raise ValueError, "The fourth input must be an SPWP bijection"
        else:
            self._A = A
            self._B = B
            self._f = f
            self._f0 = f0

    def __eq__(self,other):
        if self.get_A() != other.get_A():
            return False
        elif self.get_B() != other.get_B():
            return False
        elif self.get_SRWP() != other.get_SRWP():
            return False
        elif self.get_SPWP() != other.get_SPWP():
            return False
        else:
            return True

    def __repr__(self):
        return "This represents the reduction of " + str(list(self._A)) + " to " + str(list(self._B)) + "."

    def get_SRWP(self):
        """
        TBD
        """
        return self._f

    def get_SPWP(self):
        """
        TBD
        """
        return self._f0

    def get_A(self):
        """
        TBD
        """
        return self._A

    def get_B(self):
        """
        TBD
        """
        return self._B

    def print_involution(self):
        func = self.get_SRWP()
        for i in func.domain():
            print str(i) + ", " + str(i.get_sign()) + " --> " + str(func(i)) + ", " + str(func(i).get_sign())

    def reverse(self):
        if self.get_A().get_size() != self.get_B().get_size():
            raise ValueError, "Reduction direction cannot be reversed unless scalars are equivalent"
        else:
            A = self.get_B()
            B = self.get_A()
            f0 = inverse(self.get_SPWP())
            dic_f = dict()
            for elm in A:
                dic_f[elm] = elm
            f = FiniteSetMaps(A).from_dict(dic_f)
            return ReductionMaps(A,B,f,f0)

    def transitive(self,other=None):
        """
        TBD
        """
        if other is None:
            raise ValueError, "Enter a reduction map as a parameter"
        elif self.get_B() != other.get_A():
            raise ValueError, "These are not good candidate sets for reduction mapping"
        else:
            dic_h = dict()
            dic_h0 = dict()
            A = self.get_A()
            B = self.get_B()
            C = other.get_B()
            f = self.get_SRWP()
            f0 = self.get_SPWP()
            g = other.get_SRWP()
            g0 = other.get_SPWP()
            for elm in A:
                if elm in fixed_points(f):
                    if f0(elm) in fixed_points(g):
                        dic_h0[elm] = g0(f0(elm))
                        dic_h[elm] = elm
                    else:
                        dic_h[elm] = inverse(f0)(g(f0(elm)))
                else:
                    dic_h[elm] = f(elm)
            h = FiniteSetMaps(A).from_dict(dic_h)
            h0 = FiniteSetMaps(CombinatorialScalarWrapper(dic_h0.keys()),C).from_dict(dic_h0)
            return ReductionMaps(A,C,h,h0)

    def confluence(self,other=None):
        """
        Returns the reduction of C to B, where the usage is:
        r1.confluence(r2), and r1 maps A to B, B fully cancelled,
        and r2 maps A to C.
        """
        if other is None:
            raise ValueError, "Enter a redution map as a parameter"
        elif self.get_A() != other.get_A():
            raise ValueError, """ "A" """ " sets have to match"
        elif not(self.get_B().is_fully_cancelled()):
            raise ValueError, """ "B" """ " on the left is not fully cancelled"
        else:
            dic_h = dict()
            dic_h0 = dict()
            A = self.get_A()
            B = self.get_B()
            C = other.get_B()
            f = self.get_SRWP()
            f0 = self.get_SPWP()
            g = other.get_SRWP()
            g0 = other.get_SPWP()
            fxd_f = fixed_points(self.get_SRWP())
            fxd_g = fixed_points(other.get_SRWP())
            for c in C:
                x = inverse(g0)(c)
                while True:
                    if x in fxd_f: #a fixed point of h
                        dic_h[c] = c
                        dic_h0[c] = f0(x)
                        break
                    else:
                        x = f(x)
                    if x in fxd_g: #not a fixed point of h
                        dic_h[c]=g0(x)
                        break
                    else:
                        x = g(x)
            h = FiniteSetMaps(C).from_dict(dic_h)
            h0 = FiniteSetMaps(CombinatorialScalarWrapper(dic_h0.values()),B).from_dict(dic_h0)
            return ReductionMaps(C,B,h,h0)
