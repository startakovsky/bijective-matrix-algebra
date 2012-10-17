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

#from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.all import SageObject
#from sage.sets.all import *
from sage.all import *

class CombinatorialScalar(set):
    """
    INPUT:
     - l an iterable consisting of CombinatorialObjects
    
    Returns l, wrapped as a combinatorial scalar
    
    EXAMPLES::
    
		sage: C = CombinatorialObject("Rock",1,[0,1,1,2]);C  
		Combinatorial scalar element with name Rock, sign 1, and monomial x2*x3*x4^2.
		
		sage: D = CombinatorialObject("Paper",-1,[1,2]);D    
		Combinatorial scalar element with name Paper, sign -1, and monomial x1*x2^2.
		
		sage: E = CombinatorialObject("Scissors",1,[1,2,0]);E
		Combinatorial scalar element with name Scissors, sign 1, and monomial x1*x2^2.
		
		sage: F = CombinatorialObject("Lizard",1,[1,2,0]);F  
		Combinatorial scalar element with name Lizard, sign 1, and monomial x1*x2^2.
		
		sage: G = CombinatorialScalar((C,D,E,F)); G                 
		Combinatorial Scalar of cardinality 4.
		
		sage: G.get_generating_function()                           
		x2*x3*x4^2 + x1*x2^2
		
		sage: G.get_sign_function()                                 
		map: Combinatorial scalar element with name Rock, sign 1, and monomial x2*x3*x4^2. -> 1, Combinatorial scalar element with name 			Paper, sign -1, and monomial x1*x2^2. -> -1, Combinatorial scalar element with name Scissors, sign 1, and monomial x1*x2^2. -> 1, 			Combinatorial scalar element with name Lizard, sign 1, and monomial x1*x2^2. -> 1
		
		sage: G.get_weight_function()                               
		map: Combinatorial scalar element with name Rock, sign 1, and monomial x2*x3*x4^2. -> x2*x3*x4^2, Combinatorial scalar element with 		name Paper, sign -1, and monomial x1*x2^2. -> x1*x2^2, Combinatorial scalar element with name Scissors, sign 1, and monomial x1*x2^2. 		-> x1*x2^2, Combinatorial scalar element with name Lizard, sign 1, and monomial x1*x2^2. -> x1*x2^2
		
		sage: G.get_size()                                        
		4
		
		sage: G.is_fully_cancelled()
		False
		
		sage: G.print_list()
		Combinatorial scalar element with name Rock, sign 1, and monomial x2*x3*x4^2.
		Combinatorial scalar element with name Paper, sign -1, and monomial x1*x2^2.
		Combinatorial scalar element with name Scissors, sign 1, and monomial x1*x2^2.
		Combinatorial scalar element with name Lizard, sign 1, and monomial x1*x2^2.


    """
    
    def __init__(self, l):
        """
        Initiates the object by iterating through the set of combinatorial elements and capturing all relevant information.
        """
        self._sign_dict = dict()
        self._weight_dict = dict()
        self._generating_function = 0
        self._size = 0
        for i in l:
        	self._sign_dict[i] = i.get_sign()
        	self._weight_dict[i] = i.get_weight()
        	self._generating_function += i.get_genfunc()
        	self._size += 1
        	self.add(i)

    def __repr__(self):
        return "Combinatorial Scalar of cardinality " + str(self._size) + "."

    def get_generating_function(self):
		"""
		Returns the generating function of the combinatorial scalar.
		"""
		return self._generating_function

    def get_sign_function(self):
        """
        Returns the sign function of the combinatorial scalar.
        """
        M = FiniteSetMaps(self,(-1,1))
        return M.from_dict(self._sign_dict)
        
    def get_weight_function(self):
		"""
		Returns the weight function of the combinatorial scalar.
		"""
		M = FiniteSetMaps(self,self._weight_dict.viewvalues())
		return M.from_dict(self._weight_dict)
		
    def get_size(self):
		"""
		Returns the cardinality of the combinatorial scalar.
		"""
		return self._size
		
		
    def is_fully_cancelled(self):
		"""
		Returns 'True' if the combinatorial scalar is fully cancelled; 'False' otherwise.
		"""
		positive_set = set()
		negative_set = set()
		for i in self:
			if i.get_sign()==1:
				positive_set.add(i.get_weight())
			else:
				negative_set.add(i.get_weight())
		if positive_set.intersection(negative_set).issubset(()):
			return True
		else:
			return False
	
    def print_list(self):
		"""
		Returns a detailed list of all objects in the combinatorial scalar.
		"""
		for i in self:
			print i.get_detail()
        	
        	
        	
        
        
        
        
        
        
        
        
