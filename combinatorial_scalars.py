r"""
Combinatorial Scalars

Combinatorial Scalars are defined in and inspired by the following article in Linear Algebra and its Applications:

Bijective Matrix Algebra: http://www.sciencedirect.com/science/article/pii/S0024379506000218

Essentially they are finite sets where each element has a weight and sign associated with it.

The sign is +1 or -1, and the weight is a monomial.

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

#TBD TBD

class CombinatorialScalar(Set):
    """
    INPUT:
     - l an iterable consisting of CombinatorialScalarElements
    
    Returns l, wrapped as a combinatorial scalar
    
    EXAMPLES::
    
        sage: 
        
    """
    
    def __init__(self, l):
        """
        TBD
        """
        self._sign_dict = dict()
        self._weight_dict = dict()
        self._generating_function = 0
        self._size = 0
        for i in self:
        	self._sign_dict[i] = i.get_sign()
        	self._weight_dict[i] = i.get_weight()
        	self._generating_function += i.get_genfunc()
        	self._size += 1

    def __repr__(self):
        return "Combinatorial Scalar of cardinality " + str(self._size) + "."

    def get_sign_funcion(self):
        """
        TBD
        """
        M = FiniteSetMaps(self,(-1,1))
        return M.from_dict(self._sign_dict())
        
	def get_weight_function(self):
		"""
		TBD
		"""
		M = FiniteSetMaps(self,set(self._weight_dict.viewvalues()))
		return M.from_dict(self._weight_dict)
		
	def get_size(self):
		"""
		TBD
		"""
		return self._size
		
	def get_generating_function(self):
		"""
		TBD
		"""
		return self._generating_function
		
	def get_list(self):
		"""
		TBD
		"""
		return "TBD"
        	
        	
        	
        
        
        
        
        
        
        
        
