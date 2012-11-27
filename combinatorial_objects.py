r"""
Combinatorial Objects

Combinatorial Objects are the elements of Combinatorial Scalars.

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
#from sage.rings import *
#from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.all import SageObject
from sage.symbolic.expression import Expression
from sage.all import *

class CombinatorialObject(SageObject):
    r"""
    INPUT:
    - object an instance of the class SageObject
    - sign a number in the set {1,-1}
    - monomial a list indicating the weight monomial (eg: [0,0,1,2,1] is x2*x3^2*x4)
    
    OUTPUT:
    Returns l, wrapped as a combinatorial class
    
    EXAMPLES::
    
        sage: C = CombinatorialObject("Cool",-1,[4,3,0,6,2]);C
		Combinatorial Object Cool, sign -1, and monomial x1^4*x2^3*x4^6*x5^2.
		
		sage: C.get_tuple()                                          
		('Cool', -1, x1^4*x2^3*x4^6*x5^2)
		
		sage: C.get_genfunc()                                        
		-x1^4*x2^3*x4^6*x5^2
		
		sage: C.get_object()                                           
		'Cool'
		
		sage: C.get_weight()
		x1^4*x2^3*x4^6*x5^2
		
		sage: C.get_sign()
		-1

    """

    		    
    def __init__(self, object, sign, weight = [0]):
    	r"""
        Initiates the object and sets monomial list to variable.
        See ``CombinatorialObject`` for full documentation.
        """
        self._object = object
        self._sign = sign
        if type(weight) == Expression:
            self._weight_monomial = weight
        else:
            self._weight_monomial = assign_weight_monomial(weight)

    def __repr__(self):
        return str(self._object)

    def __eq__(self,other):
        return self.get_object() == other.get_object() and bool(self.get_genfunc() == other.get_genfunc())

    def __hash__(self):
        return hash(self.get_tuple())
        
    def get_detail(self):
        return "Combinatorial Object %s, sign %d, and weight %s." %(self._object, self._sign, str(self._weight_monomial))

    def get_tuple(self):
    	return (self.get_object(),self.get_sign(),self.get_weight())

    def get_genfunc(self):
    	return self._sign * self._weight_monomial
    
    def get_object(self):
    	return self._object
    
    def set_object(self,obj):
        self._object = obj
        return self
    	
    def get_weight(self):
    	return self._weight_monomial
    	
    def get_sign(self):
    	return self._sign
    	
    def get_cleaned_up_version(self):
        r"""
        Returns a new cleaned up version of this object.
        Eliminates nested tuples.  Does not change actual object.
        """
        t = deepcopy(self)
        t.set_object(clean_up_object(self.get_object()))
        return t

def create_variables(n):
    r"""
    Given an integer n, this method returns the string 'x1,..,xn'
    """
    variables = ''
    for i in range(n):
        variables+= 'x' + str(i+1) + ','
    variables = variables[0:len(variables)-1]
    return variables

def assign_weight_monomial(l):
    r"""
    Returns weight monomial from list l.
    """
    l.append(0) #to handle the case of a power of x1
    v = var(create_variables(len(l)))
    monomial = 1
    for i in range(len(l)):
        monomial = monomial * v[i]**l[i]
    return monomial
	
def clean_up_object(l):
    r"""
    Returns just the objects inside nested tuples and preserves the order.
    """
    K = list()
    if type(l)==tuple:
        for i in range(len(l)):
            if type(l[i]) == tuple:
                K.extend(clean_up_object(l[i]))
            else:
                K.append(l[i])
        return tuple(K)
    else:
        return l
