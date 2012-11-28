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
from sage.combinat.permutation import *
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.bijectivematrixalgebra.reduction_maps_dicts import ReductionMapsDict
from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps
from copy import copy

def _involution_dict(mat):
    r"""
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

def reduction_matrix_clean_up(mat, st="standard clean up"):
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
            f = FiniteSetMaps(A[i,j],A[i,j]).from_dict(dic_f)
            f0 = FiniteSetMaps(A[i,j],B[i,j]).from_dict(dic_f0)
            d[i,j] = ReductionMaps(A[i,j],B[i,j],f,f0)
    return ReductionMapsDict(d,st)

def reduction_identity_matrix(mat, st = None):
    r"""
    When a matrix reduces to the identity, this returns
    a ReductionMapDict of from a matrix to I.
    This should be used when creating arbitrary SRWP involutions.
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
                f0s[i,j] = FiniteSetMaps(set(),set()).from_dict({})
    d = dict()
    for i in range(dim):
        for j in range(dim):
            d[i,j] = ReductionMaps(mat[i,j],I[i,j],fs[i,j],f0s[i,j])
    return ReductionMapsDict(d,st)

def reduction_lemma_40(mat, st = "lemma 40"):
    r"""
    Returns the reduction of mat = adj_A times A
    to det_A times I.
    """
    dim = mat.nrows()
    d = dict()
    B = matrix_adjoint_lemma_40(mat)
    A = mat
    for i in range(dim):
        for j in range(dim):
            dic_f = dict()
            dic_f0 = dict()
            copyset = deepcopy(A[i,j].get_set())
            if i==j:
                for elm in copyset:
                    tmp = list(elm.get_object()[0].get_object())
                    #object is tuple, elements come from the actual tuple, hence double get_object()
                    index = tmp.index(CombinatorialObject('_',1))
                    tmp[index]=elm.get_object()[1]
                    dic_f[elm] = elm
                    copyelm= deepcopy(elm)
                    dic_f0[elm] = copyelm.set_object(tuple(tmp))
                f0 = FiniteSetMaps(A[i,j],B[i,j]).from_dict(dic_f0)
            else:
                for elm in copyset:
                    tmp = list(elm.get_object()[0].get_object())
                    #object is tuple, elements come from the actual tuple, hence double get_object()
                    p = tmp[0].get_object()
                    ii = i+1
                    jj = j+1
                    q = Permutation((ii,jj))
                    pq = q*p #p composed with q
                    tmp[0] = CombinatorialObject(pq,pq.signature())
                    elm_range2 = tmp[jj]
                    tmp[jj] = elm.get_object()[1]
                    sgn = 1
                    weight = 1
                    for flop in tmp:
                        sgn *= flop.get_sign()
                        weight *= flop.get_weight()
                    elm_range1 = CombinatorialObject(tuple(tmp),sgn,weight)
                    elm_range = CombinatorialObject((elm_range1,elm_range2),elm_range1.get_sign()*elm_range2.get_sign(),elm_range1.get_weight()*elm_range2.get_weight())
                    dic_f[elm] = elm_range
                f0 = FiniteSetMaps(set()).from_dict({})
            f = FiniteSetMaps(A[i,j]).from_dict(dic_f)
            d[i,j] = ReductionMaps(A[i,j],B[i,j],f,f0)
    return ReductionMapsDict(d,st)
