r"""
Combinatorial Scalar Ring and Ring Element.

The main goal with this is to make a ring out of Combinatorial Scalars with operations such as multiplication and addition and uses a __getattr__ special method to make this behave as a set, and add an iterator.  

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
from sage.structure.unique_representation import UniqueRepresentation


class CombinatorialScalarWrapper(RingElement):
    def __init__(self,_set):
        RingElement.__init__(self,CombinatorialScalarRing())
        self.values = CombinatorialScalar(_set)
        self.iterator = self.__iter__()
        if self.values != set():
            self.next = self.iterator.next()
    def __iter__(self):
        for this_entry in self.values:
            yield this_entry
    def get_set(self):
        return set(self.values)
    def get_scalar(self):
        return self.values
    def __repr__(self):
        return str(list(self.values))
    def __getattr__(self,attr):
        return getattr(self.values,attr)

class CombinatorialScalarWrapper(CombinatorialScalarWrapper):
    """
    Defines plus and times operations as well as special equality boolean method.
    """
    def __eq__(self,other):
        return self.get_set()==other.get_set() and self.parent() == other.parent()
    def __add__(self,other):
        return CombinatorialScalarWrapper(self.get_set().union(other.get_set()))
    def __mul__(self,other):
        new_set = set()
        for s in self:
            for o in other:
                new_set.add(CombinatorialObject((s.get_object(),o.get_object()),s.get_sign()*o.get_sign(),s.get_weight()*o.get_weight()))
        return CombinatorialScalarWrapper(new_set)
        
class CombinatorialScalarRing(Ring,UniqueRepresentation):
    """
    TBD
    """
    def __init__(self):
        Parent.__init__(self,category = Rings())
    def __repr__(self):
        return "Combinatorial Scalar Ring"
    def __eq__(self,other):
        return type(other) == type(self)
    def _zero(self):
        return CombinatorialScalarWrapper(set())
    def _one(self):
        return CombinatorialScalarWrapper([CombinatorialObject(1,1)])
    def _element_constructor_(self,i):
        return CombinatorialScalarWrapper(i)
    Element = CombinatorialScalarWrapper
