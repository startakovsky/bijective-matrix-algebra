r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It creates all the foundational methods around matrices, including multiplication,
determinant, adjoint and printing the matrix.

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
from sage.matrix.all import matrix
from sage.matrix.all import MatrixSpace
from sage.combinat.permutation import *
from sage.combinat.cartesian_product import CartesianProduct
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.sets.finite_set_maps import FiniteSetMaps
from copy import copy
from copy import deepcopy


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
    Return determinant scalar, the form of which is:
    (\sigma,a_1,...,a_n) where sign \sigma is sgn(\sigma)
    and weight \sigma is 1.
    """
    dim = mat.nrows()
    P = Permutations(dim)
    S = set()
    for p in P:
        l = list()
        p_comb = CombinatorialScalarWrapper([CombinatorialObject(p,p.signature())])
        for i in range(1,dim+1):
            l.append(mat[i-1,p(i)-1])
        cp = CartesianProduct(p_comb,*l)
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

def matrix_print(mat):
    print "Printing..."
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            print "row " + str(i) + ", column " + str(j) + "; " + str(mat[i,j].get_size()) + " elements"
            mat[i,j].print_list()
            print "------------------------------"

def matrix_comparison(matA,matB):
    nrows = matA.nrows()
    ncols = matB.ncols()
    if matA.ncols()!=matB.ncols() or matA.nrows()!=matB.nrows():
        raise ValueError, "Check dimensions of input"
    else:
        for i in range(nrows):
            for j in range(ncols):
                if matA[i,j]!=matB[i,j]:
                    return False
        return True
        
def matrix_multiply_scalar(mat,scal):
    nrows = mat.nrows()
    ncols = mat.ncols()
    L = list()
    for i in range(nrows):
        L.append(list())
        for j in range(ncols):
            L[i].append(scal*mat[i,j])
    return matrix(CombinatorialScalarRing(),nrows,ncols,L)

def matrix_adjoint_lemma_40(mat):
    """
    Returns a matrix which is the target in the
    reduction of adjoint(A) times A to det(A) times I.
    """
    if mat.nrows()!=mat.ncols():
        raise ValueError, "Make sure that this is indeed a Combinatorial Adjoint matrix"
    else:
        dim = mat.nrows()
        L = list()
        mat_space = MatrixSpace(CombinatorialScalarRing(),dim)
        for i in range(dim):
            L.append(list())
            for j in range(dim):
                if i==j:
                    copyset = deepcopy(mat[i,j].get_set())
                    for elm in copyset:
                        tmp = list(elm.get_object()[0].get_object()) 
                        #object is tuple, elements come from the actual tuple, hence double get_object()
                        index = tmp.index(CombinatorialObject('_',1))
                        tmp[index]=elm.get_object()[1]
                        elm.set_object(tuple(tmp))
                else:
                    copyset = CombinatorialScalarWrapper(set())
                L[i].append(CombinatorialScalarWrapper(copyset))
        return  mat_space(L)
