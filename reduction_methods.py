r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It generates the standard reductions discussed in the Loehr-Mendes paper
and pulls from the ReductionMaps and ReductionMapsDict classes.

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
#                 http://www.gnu.org/licenses/
#*****************************************************************************

from sage.bijectivematrixalgebra.matrix_methods import *
from sage.matrix.all import matrix
from sage.matrix.all import MatrixSpace
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.bijectivematrixalgebra.reduction_maps_dicts import ReductionMapsDict
from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps
from copy import copy

def _involution_dict(mat):
    """
    Returns a dictionary of arbitrary involutions on the entries of a Combinatorial Matrix.
    Note all weights must be 1 for this.  It may be extended upon in the future.
    """
    mat_gen_func = matrix_generating_function(mat)
    if mat_gen_func != mat_gen_func.parent().identity_matrix():
        raise ValueError, "Input needs to be equal to the identity."
    else:
        func = dict()
        for x in range(mat.nrows()):
            for y in range(mat.ncols()):
                if x <> y:
                    func[(x,y)] = mat[x,y].create_involution()
                else:
                    t = mat[x,y]
                    for i in t:
                        _M = FiniteSetMaps(t)
                        func[(x,y)] = _M.from_dict({i:i})
        return func

def matrix_clean_up_reduction(mat, st="standard clean up"):
    dim = mat.nrows()
    d = dict()
    B = matrix_clean_up(mat)
    A = mat
    for i in range(dim):
        for j in range(dim):
            dic_f = dict()
            dic_f0 = dict()
            for elm in A[i,j]:
                dic_f[elm] = elm
                dic_f0[elm] = elm.get_cleaned_up_version()
            f = FiniteSetMaps(A[i,j]).from_dict(dic_f)
            f0 = FiniteSetMaps(A[i,j],B[i,j]).from_dict(dic_f0)
            d[i,j] = ReductionMaps(A[i,j],B[i,j],f,f0)
    return ReductionMapsDict(d,st)

def matrix_identity_reduction(mat, st = None):
    """
    When a matrix reduces to the identity, this returns
    a ReductionMapDict of from a matrix to I.
    """
    dim = mat.nrows()
    fs = _involution_dict(mat)
    f0s = dict()
    I = identity_matrix(dim)
    for i in range(dim):
        for j in range(dim):
            if i==j:
                f0s[i,j] = FiniteSetMaps(mat[i,j],I[i,j]).from_dict({mat[i,j].get_set().pop():CombinatorialObject(1,1)})
            else:
                f0s[i,j] = FiniteSetMaps(set()).from_dict({})
    d = dict()
    for i in range(dim):
        for j in range(dim):
            d[i,j] = ReductionMaps(mat[i,j],I[i,j],fs[i,j],f0s[i,j])
    return ReductionMapsDict(d,st)
