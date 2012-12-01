r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It generates the standard combinatorial interpretations of Stirling Matrices.

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


from sage.matrix.all import matrix
from sage.matrix.all import MatrixSpace
from sage.combinat.permutation import *
from sage.combinat.set_partition import *
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing

PermutationOptions(display = 'cycle')
PermutationOptions(display = 'singleton')

def _stirling1_row(row,dim):
    r = list()
    for i in range(dim):
        r.append(set())
    P = Permutations(row)
    for p in P:
        t = CombinatorialObject(p,(-1)**(row-len(p.to_cycles())))
        r[len(p.to_cycles())].add(t)
    for j in range(dim):
        r[j] = CombinatorialScalarWrapper(r[j])
    return r

def _stirling2_row(row,dim):
    r = list()
    for j in range(dim):
        r.append(set())
    for p in SetPartitions(row):
        t = CombinatorialObject(p,1)
        r[len(p)].add(t)
    for j in range(dim):
        r[j] = CombinatorialScalarWrapper(r[j])
    return r


def Stirling1Matrix(dim):
    r"""
    Returns Stirling1 Matrix whose entries are Combinatorial Scalars of signed permutations.
    """
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    l = list()
    for row in range(dim):
        l.append(_stirling1_row(row,dim))
    return mat_space(l)

def Stirling2Matrix(dim):
    r"""
    Returns Stirling2 Matrix whose entries are Combinatorial Scalars of set partitions.
    """ 
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    l = list()
    for row in range(dim):
        l.append(_stirling2_row(row,dim))
    return mat_space(l)

def find_row(elm):
    obj = elm.get_object()
    if type(obj) == Permutation_class:
        return obj.size()
    elif type(obj) == type(SetPartitions(0).random_element()):
        newset = set([0])
        for i in obj:
            newset.add(max(i))
        return max(newset)
    else:
        raise ValueError, "This algorithm expects SetPartitions or Permutations"

def find_col(elm):
    obj = elm.get_object()
    if type(obj) == Permutation_class:
        return len(obj.cycle_type())
    elif type(obj) == type(SetPartitions(0).random_element()):
        return len(obj)
    else:
        raise ValueError, "This algorithm expects SetPartitions or Permutations"
