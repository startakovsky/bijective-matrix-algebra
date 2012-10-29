r"""
Combinatorial Scalar Ring

See combinatorial_scalars for background.

The primary purpose of this class is to be able to treat Combinatorial Scalars as entries in a matrix.

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
from copy import deepcopy
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar

class CombinatorialScalarWrapper(RingElement):
    """
    TBD
    """
    wrapped_class = CombinatorialScalar
    def __repr__(self):
        return str(list(self.value))
    def __add__(self,other):
        return CombinatorialScalarWrapper(CombinatorialScalar(self.value.union(deepcopy(other.value))), parent = self.parent())
    def __mul__(self,other):
        new_set = set()
        for s in self.value:
            for o in other.value:
                new_set.add(CombinatorialObject((s.get_object(),o.get_object()),s.get_sign()*o.get_sign(),s.get_weight()*o.get_weight()))
        return CombinatorialScalarWrapper(CombinatorialScalar(new_set), parent = self.parent())
        
        
class CombinatorialScalarRing(Ring):
    """
    TBD
    """
    def __init__(self):
        Parent.__init__(self,category = Rings())
    def __repr__(self):
        return "Combinatorial Scalar Ring"
    def __eq__(self,other):
        return type(other) == type(self)
    Element = CombinatorialScalarWrapper