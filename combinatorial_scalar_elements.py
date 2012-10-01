r"""
Combinatorial Scalar Element

Combinatorial Scalar Elements are the elements of Combinatorial Scalars.

Each element has a weight and sign associated with it.

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
from sage.rings import *
from sage.bijectivematrixalgebra.main import *
from sage.structure.sage_object import SageObject

class CombinatorialScalarElement(SageObject):
    """
    INPUT:
     - element_name a string indicating the name of the element
     - sign a number in the set {1,-1}
     - monomial a list indicating the weight monomial (eg: [0,0,1,2,1] is x2*x3^2*x4)
    
    OUTPUT:
    Returns l, wrapped as a combinatorial class
    
    EXAMPLES::
    
        sage: C = CombinatorialScalarElement("Cool",-1,[4,3,0,6,2]);C
		Combinatorial scalar element with name Cool, sign -1, and monomial x1^4*x2^3*x4^6*x5^2.
		
		sage: C.get_tuple()                                          
		('Cool', -1, x1^4*x2^3*x4^6*x5^2)
		
		sage: C.get_genfunc()                                        
		-x1^4*x2^3*x4^6*x5^2
		
		sage: C.get_name()                                           
		'Cool'

    """

    		    
    def __init__(self, element_name, sign, weight):
        """
        Initiates the object and sets monomial list to variable.
        See ``CombinatorialScalarElement`` for full documentation.
        """
        self._element_name = element_name
        self._sign = sign
        self._weight_monomial = assign_weight_monomial(weight)
        
    def __repr__(self):
        return "Combinatorial scalar element with name %s, sign %d, and monomial %s." %(self._element_name, self._sign, str(self._weight_monomial))
    	
    def get_tuple(self):
    	return (self._element_name,self._sign,self._weight_monomial)
    
    def get_genfunc(self):
    	return self._sign * self._weight_monomial
    
    def get_name(self):
    	return self._element_name