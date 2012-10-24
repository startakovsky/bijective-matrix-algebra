r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It generates the standard combinatorial interpretations of Stirling Matrices and their product.

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

import sage.rings.rational
from sage.matrix.all import *
from sage.combinat.permutation import *
from sage.combinat.set_partition import *
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing

PermutationOptions(display = 'cycle')
PermutationOptions(display = 'singleton')

def _stirling1_row(row,dim,prnt):
    r = list()
    for i in range(dim):
        r.append(set())
    P = Permutations(row)
    for p in P:
        t = CombinatorialObject(p,(-1)**(row-len(p.to_cycles())))
        r[len(p.to_cycles())].add(t)
    for j in range(dim):
        r[j] = CombinatorialScalarWrapper(CombinatorialScalar(r[j]),parent = prnt)
    return r

def _stirling2_row(row,dim, prnt):
    r = list()
    for j in range(dim):
        r.append(set())
    for p in SetPartitions(row):
        t = CombinatorialObject(p,1)
        r[len(p)].add(t)
    for j in range(dim):
        r[j] = CombinatorialScalarWrapper(CombinatorialScalar(r[j]),parent = prnt)
    return r

def _product_row(mat1, mat2, row):
    dim = mat1.nrows()
    prnt = mat1.matrix_space().base_ring()
    r = list()
    for j in range(dim):
        C = CombinatorialScalarWrapper(CombinatorialScalar(set()), parent = prnt)
        for k in range(dim):
            C = C + (mat1[row,k]*mat2[k,j])
        r.append(C)
    return r

def matrix_multiply(mat1,mat2):
    """
    Only works for square matrices.
    """
    mat_space = mat1.matrix_space()
    dim = mat_space.nrows()
    l = list()
    for row in range(dim):
        l.append(_product_row(mat1,mat2,row))
    return mat_space(l)

def Stirling1Matrix(dim):
    """
    Returns Stirling1 Matrix whose entries are Combinatorial Scalars of signed permutations.
    """
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    prnt = mat_space.base_ring()
    l = list()
    for row in range(dim):
        l.append(_stirling1_row(row,dim,prnt))
    return mat_space(l)

def Stirling2Matrix(dim):
    """
    Returns Stirling2 Matrix whose entries are Combinatorial Scalars of set partitions.
    """ 
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    prnt = mat_space.base_ring()
    l = list()
    for row in range(dim):
        l.append(_stirling2_row(row,dim,prnt))
    return mat_space(l)

def matrix_generating_function(m):
    """
    returns the generating function of each scalar as a matrix
    """
    dimx = m.nrows()
    dimy = m.ncols()
    d = dict()
    mat_space = MatrixSpace(CombinatorialScalarRing(),dimx,dimy)
    for x in range(dimx):
        for y in range(dimy):
            d[(x,y)]=m[x,y].value.get_generating_function()
    return mat_space(QQ,d)