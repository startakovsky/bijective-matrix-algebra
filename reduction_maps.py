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
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar
from sage.bijectivematrixalgebra.main import fixed_points
from sage.bijectivematrixalgebra.main import is_SRWP_involution
from sage.bijectivematrixalgebra.main import is_SPWP_bijection
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.structure.sage_object import SageObject
from sage.sets.finite_set_map_cy import FiniteSetEndoMap_Set
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
        elif not(bool(A.value.get_generating_function()==B.value.get_generating_function())):
            raise ValueError, "The generating functions of the scalars are not equal"
        elif type(f) != FiniteSetEndoMap_Set:
            raise ValueError, "The third input must be a map in FiniteSetMaps"
        elif type(f0) != FiniteSetEndoMap_Set:
            raise ValueError, "The fourth input must be a map in FiniteSetMaps"
        elif set(f.domain()) != A.value:
            raise ValueError, "The third input must have domain of first input"
        elif set(f0.domain()) != fixed_points(f):
            raise ValueError, "The fourth input must have domain of fixed points of third input"
        elif set(f0.image_set())!=B.value:
            raise ValueError, "The fourth input must have image of second input"
        elif not(is_SRWP_involution(f)):
            raise ValueError, "The third input must be an SRWP involution"
        elif not(is_SPWP_bijection(f0)):
            raise ValueError, "The fourth input must be an SPWP bijection"
        else:
            self._A = A.value
            self._B = B.value
            self._f = f
            self._f0 = f0

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

    def transitivity(self,other=None):
		"""
		TBD
		"""
		return "This needs to be built out."
		

    def confluence(self,other=None):
		"""
		TBD
		"""	
		return "This needs to be built out."