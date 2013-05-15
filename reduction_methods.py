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
from sage.bijectivematrixalgebra.map_methods import fixed_points
from sage.bijectivematrixalgebra.reduction_maps_dicts import ReductionMapsDict
from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps
import sage.bijectivematrixalgebra.stirling as stirling
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
                        _M = FiniteSetMaps(t,t)
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

def reduction_identity_matrix(mat,involution_dict=None, st = None):
    r"""
    When a matrix reduces to the identity, this returns
    a ReductionMapDict of from a matrix to I.
    """
    dim = mat.nrows()
    if involution_dict is None:
        fs = _involution_dict
    else:
        fs = involution_dict
    f0s = dict()
    I = identity_matrix(dim)
    for i in range(dim):
        for j in range(dim):
            if i==j:
                tmp = CombinatorialScalarWrapper(fixed_points(fs[i,j]))
                f0s[i,j] = FiniteSetMaps(tmp,I[i,j]).from_dict({tmp.get_set().pop():CombinatorialObject(1,1)})
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
                f0 = FiniteSetMaps(set(),set()).from_dict({})
            f = FiniteSetMaps(A[i,j],A[i,j]).from_dict(dic_f)
            d[i,j] = ReductionMaps(A[i,j],B[i,j],f,f0)
    return ReductionMapsDict(d,st)

def reduction_matrix_AIB_AB(mat,st = "remove middle Identity matrix"):
    r"""
    Input AIB and output is AB, that is, we remove the middle 
    Combinatorial Object 1 from each triple and return only
    the outside two entries.
    """
    if mat.nrows()!=mat.ncols():
        raise ValueError, "Check dimensions"
    else:
        dim = mat.nrows()
        d = dict()
        for i in range(dim):
            for j in range(dim):
                dic_f0 = dict()
                dic_f = dict()
                newset = set()
                for elm in mat[i,j]:
                    newelm0 = elm.get_object()[0]
                    newelm1 = elm.get_object()[2]
                    tmp = newelm0*newelm1
                    newset.add(tmp)
                    dic_f[elm] = elm
                    dic_f0[elm] = tmp
                f = FiniteSetMaps(mat[i,j],mat[i,j]).from_dict(dic_f)
                f0 = FiniteSetMaps(mat[i,j],newset).from_dict(dic_f0)
                d[i,j] = ReductionMaps(mat[i,j],CombinatorialScalarWrapper(newset),f,f0)
        return ReductionMapsDict(d,st)

def reduction_matrix_IAB_AB(mat,st = "remove left Identity matrix"):
    r"""
    Input IAB and output is AB, that is, we remove the left 
    Combinatorial Object 1 from each triple and return only
    the other two entries.
    """
    if mat.nrows()!=mat.ncols():
        raise ValueError, "Check dimensions"
    else:
        dim = mat.nrows()
        d = dict()
        for i in range(dim):
            for j in range(dim):
                dic_f0 = dict()
                dic_f = dict()
                newset = set()
                for elm in mat[i,j]:
                    newelm1 = elm.get_object()[1]
                    newelm2 = elm.get_object()[2]
                    tmp = newelm1*newelm2
                    newset.add(tmp)
                    dic_f[elm] = elm
                    dic_f0[elm] = tmp
                f = FiniteSetMaps(mat[i,j],mat[i,j]).from_dict(dic_f)
                f0 = FiniteSetMaps(mat[i,j],newset).from_dict(dic_f0)
                d[i,j] = ReductionMaps(mat[i,j],CombinatorialScalarWrapper(newset),f,f0)
        return ReductionMapsDict(d,st)

def reduction_matrix_ABCD_to_ApBCpD(A,B,C,D,st = None):
    r"""
    returns the reduction/equivalence of the product
    of ABCD to A(BC)D.
    """
    mat = matrix_multiply(A,matrix_multiply(matrix_multiply(B,C),D))
    dim = mat.nrows()
    d = dict()
    for i in range(dim):
        for j in range(dim):
            newsetA = set()
            newsetB = set()
            dic_f = dict()
            dic_f0 = dict()
            for elm in mat[i,j]:
                tmp0 = elm.get_object()[0]
                tmp1 = elm.get_object()[1].get_object()[0]
                tmp2 = elm.get_object()[1].get_object()[1]
                tmpB = CombinatorialObject((tmp0,tmp1,tmp2),elm.get_sign(),elm.get_weight())
                newsetB.add(tmpB)
                tmp10 = tmp1.get_object()[0]
                tmp11 = tmp1.get_object()[1]
                tmpA = CombinatorialObject((tmp0,tmp10,tmp11,tmp2),elm.get_sign(),elm.get_weight())
                newsetA.add(tmpA)
                dic_f[tmpA] = tmpA
                dic_f0[tmpA] = tmpB
            f = FiniteSetMaps(newsetA,newsetA).from_dict(dic_f)
            f0 = FiniteSetMaps(newsetA,newsetB).from_dict(dic_f0)
            d[i,j] = ReductionMaps(CombinatorialScalarWrapper(newsetA),CombinatorialScalarWrapper(newsetB),f,f0)
    return ReductionMapsDict(d,st)
    
def reduction_matrix_ABCD_to_pABpCD(A,B,C,D,st = None,reduction = None):
    r"""
    returns the reduction/equivalence of the product
    of ABCD to (AB)CD.
    """
    if reduction == None:
        mat = matrix_multiply(matrix_multiply(matrix_multiply(A,B),C),D)
    else: #not used... yet
        mat = reduction.get_matrix_A()
    d = dict()
    dim = mat.nrows()
    for i in range(dim):
        for j in range(dim):
            newsetA = set()
            newsetB = set()
            dic_f = dict()
            dic_f0 = dict()
            for elm in mat[i,j]:
                tmp0 = elm.get_object()[0].get_object()[0]
                tmp01 = tmp0.get_object()[0]
                tmp02 = tmp0.get_object()[1]
                tmp1 = elm.get_object()[0].get_object()[1]
                tmp2 = elm.get_object()[1]
                tmpB = CombinatorialObject((tmp0,tmp1,tmp2),elm.get_sign(),elm.get_weight())
                newsetB.add(tmpB)
                tmpA = CombinatorialObject((tmp01,tmp02,tmp1,tmp2),elm.get_sign(),elm.get_weight())
                newsetA.add(tmpA)
                dic_f[tmpA] = tmpA
                dic_f0[tmpA] = tmpB
            f = FiniteSetMaps(newsetA,newsetA).from_dict(dic_f)
            f0 = FiniteSetMaps(newsetA,newsetB).from_dict(dic_f0)
            d[i,j] = ReductionMaps(CombinatorialScalarWrapper(newsetA),CombinatorialScalarWrapper(newsetB),f,f0)
    return ReductionMapsDict(d,st)

def reduction_lemma_28_23(mat, red_AB_to_I, st = "an application of lemma 28, reduction_23"):
    r"""
    Because only one matrix here has a nontrivial SRWP map,
    we need not apply the formal indexing given in the proof
    of lemma 28.  Simply enter a matrix adj_A(AB)A and the 
    reduction of AB to I.
    """
    d = dict()
    dim = mat.nrows()
    for i in range(dim):
        for j in range(dim):
            dic_f = dict()
            dic_f0 = dict()
            newset = set()
            for elm in mat[i,j]:
                #break tuple apart
                tmp0 = elm.get_object()[0]
                tmp1 = elm.get_object()[1]
                tmp2 = elm.get_object()[2]
                row = tmp1.get_row()
                col = tmp1.get_col()
                f_row_col = red_AB_to_I[row,col].get_SRWP()
                f0_row_col = red_AB_to_I[row,col].get_SPWP()
                #assign map
                tmp = f_row_col(tmp1)
                sign = tmp.get_sign()*tmp0.get_sign()*tmp2.get_sign()
                weight = tmp.get_weight()*tmp0.get_weight()*tmp2.get_weight()
                dic_f[elm] = CombinatorialObject((tmp0,tmp,tmp2),sign,weight)
                if tmp1 in fixed_points(f_row_col):
                    tmpfxd = CombinatorialObject((tmp0,f0_row_col(tmp1),tmp2),elm.get_sign(),elm.get_weight())
                    dic_f0[elm] = tmpfxd
                    newset.add(tmpfxd)
            f = FiniteSetMaps(mat[i,j],mat[i,j]).from_dict(dic_f)
            f0 = FiniteSetMaps(dic_f0.keys(),newset).from_dict(dic_f0)
            d[i,j] = ReductionMaps(mat[i,j],CombinatorialScalarWrapper(newset),f,f0)
    return ReductionMapsDict(d,st)

def reduction_lemma_28_68(mat, red_adjAA_to_I, st = "an application of lemma 28, reduction_68"):
    r"""
    Because only one matrix here has a nontrivial SRWP map,
    we need not apply the formal indexing given in the proof
    of lemma 28.  Simply enter a matrix (adj_AA)BA and the
    reduction of adj_AA to I.
    """
    d = dict()
    dim = mat.nrows()
    for i in range(dim):
        for j in range(dim):
            dic_f = dict()
            dic_f0 = dict()
            newset = set()
            for elm in mat[i,j]:
                #break tuple apart
                tmp0 = elm.get_object()[0]
                tmp1 = elm.get_object()[1]
                tmp2 = elm.get_object()[2]
                row = tmp0.get_row()
                col = tmp0.get_col()
                f_row_col = red_adjAA_to_I[row,col].get_SRWP()
                f0_row_col = red_adjAA_to_I[row,col].get_SPWP()
                #assign map
                tmp = f_row_col(tmp0)
                sign = tmp.get_sign()*tmp1.get_sign()*tmp2.get_sign()
                weight = tmp.get_weight()*tmp1.get_weight()*tmp2.get_weight()
                dic_f[elm] = CombinatorialObject((tmp,tmp1,tmp2),sign,weight)
                if tmp in fixed_points(f_row_col):
                    tmpfxd = CombinatorialObject((f0_row_col(tmp0),tmp1,tmp2),elm.get_sign(),elm.get_weight())
                    dic_f0[elm] = tmpfxd
                    newset.add(tmpfxd)
            f = FiniteSetMaps(mat[i,j],mat[i,j]).from_dict(dic_f)
            f0 = FiniteSetMaps(dic_f0.keys(),newset).from_dict(dic_f0)
            d[i,j] = ReductionMaps(mat[i,j],CombinatorialScalarWrapper(newset),f,f0)
    return ReductionMapsDict(d,st)
