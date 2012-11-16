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


from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.matrix.all import *
from sage.combinat.permutation import *
from sage.combinat.set_partition import *
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.bijectivematrixalgebra.reduction_maps_dicts import ReductionMapsDict
from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps
from copy import copy

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

def _product_row(mat1, mat2, row):
    dim = mat1.nrows()
    r = list()
    for j in range(dim):
        C = CombinatorialScalarWrapper(set())
        for k in range(dim):
            C = C + (mat1[row,k]*mat2[k,j])
        r.append(C)
    return r

def identity_matrix(dim):
    """
    Returns standard combinatorial identity matrix
    """
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    prnt = mat_space.base_ring()
    L = list()
    for i in range(dim):
        L.append(list())
        for j in range(dim):
            if i==j:
                L[i].append(prnt._one())
            else:
                L[i].append(prnt._zero())
    return mat_space(L)

def _involution_dict(mat):
    """
    Returns a dictionary of arbitrary involutions on the entries of a Combinatorial Matrix.
    Note all weights must be 1 for this.  It may be extended upon in the future.
    """
    mat_gen_func = matrix_generating_function(mat)
    if mat_gen_func != mat_gen_func.parent().identity_matrix():
        print "ERROR: Input needs to be equal to the identity."
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
    l = list()
    for row in range(dim):
        l.append(_stirling1_row(row,dim))
    return mat_space(l)

def Stirling2Matrix(dim):
    """
    Returns Stirling2 Matrix whose entries are Combinatorial Scalars of set partitions.
    """ 
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    l = list()
    for row in range(dim):
        l.append(_stirling2_row(row,dim))
    return mat_space(l)

def matrix_generating_function(m):
    """
    returns the generating function of each scalar as a matrix
    """
    num_var = 10
    dimx = m.nrows()
    dimy = m.ncols()
    d = dict()
    for x in range(dimx):
        for y in range(dimy):
            d[(x,y)]=m[x,y].get_generating_function()
    return matrix(PolynomialRing(ZZ,['x'+str(i) for i in range(num_var)]),d)

def matrix_remove_row_col(mat,row,col):
    """
    Return matrix with row, col removed
    """
    L = list()
    newrows = range(mat.nrows())
    newcols= range(mat.ncols())
    newrows.pop(row)
    newcols.pop(col)
    for x in newrows:
        L.append(list())
        for y in newcols:
            L[len(L)-1].append(mat[x,y])
    return matrix(mat.parent().base_ring(),len(newrows),len(newcols),L)

def matrix_determinant(mat):
    """
    Return determinant scalar
    """
    dim = mat.nrows()
    P = Permutations(dim)
    S = set()
    for p in P:
        l = list()
        sgn = CombinatorialScalarWrapper([CombinatorialObject(p.signature(),p.signature())])
        for i in range(1,dim+1):
            l.append(mat[i-1,p(i)-1])
        cp = CartesianProduct(sgn,*l)
        for i in cp:
            weight = 1
            sign = 1
            for elm in i:
                sign = sign*elm.get_sign()
                weight = weight*elm.get_weight()
            S.add(CombinatorialObject(tuple(i),sign,weight))
    return CombinatorialScalarWrapper(S)

def matrix_combinatorial_adjoint(mat):
    """
    Return Combinatorial Adjoint.
    """
    dim = mat.nrows()
    M = list()
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    prnt = mat_space.base_ring()
    #create list L of lists of sets, dimension is increased by one to mitigate index confusion
    L = list()
    for i in range(dim+1):
        L.append(list())
        for j in range(dim+1):
            L[i].append(set())
    P = Permutations(dim)
    for p in P:
        p_comb = CombinatorialScalarWrapper([CombinatorialObject(p,p.signature())])
        l = list()
        for i in range(1,dim+1):
            l.append(mat[p(i)-1,i-1])
        #This list will have empty sets, which will yield to an empty cartesian product
        #especially when the matrix input is triangular (except for the identity permutation).
        #We will now iterate through the selected entries in each column
        #and create a set of a singleton of an empty string that corresponds 
        #to the "missing" element of the tuple described in definition 39.
        for i in range(1,dim+1):
            copy_l = copy(l)
            copy_l[i-1]=CombinatorialScalarWrapper([CombinatorialObject('_',1)])
            cp = CartesianProduct(p_comb,*copy_l)
            for tupel in cp:
                tupel = tuple(tupel)
                tupel_weight = 1
                tupel_sign = 1
                for elm in tupel:
                    tupel_sign = tupel_sign*elm.get_sign()
                    tupel_weight = tupel_weight*elm.get_weight()
                L[i][p(i)].add(CombinatorialObject(tupel,tupel_sign,tupel_weight))
    #turn these sets into CombinatorialScalars
    for i in range(1,dim+1):
        l = list()
        for j in range(1,dim+1):
            l.append(CombinatorialScalarWrapper(L[i][j]))
        M.append(l)
    return mat_space(M)

def matrix_clean_up(mat):
    """
    Apply get_cleaned_up_version to each object within
    each entry and return a cleaned up matrix
    """
    dim = mat.nrows()
    L = list()
    mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
    prnt = mat_space.base_ring()
    for i in range(dim):
        L.append(list())
        for j in range(dim):
            L[i].append(CombinatorialScalarWrapper(mat[i,j].get_cleaned_up_version()))
    return mat_space(L)

def matrix_clean_up_reduction(mat):
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
    return ReductionMapsDict(d,"clean up")

def matrix_print(mat):
    print "Printing..."
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            print "row " + str(i) + ", column " + str(j) + "; " + str(mat[i,j].get_size()) + " elements"
            mat[i,j].print_list()
            print "------------------------------"

def matrix_identity_reduction(mat):
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
    return ReductionMapsDict(d," an arbitrary involution")
