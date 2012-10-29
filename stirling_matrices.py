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

from sage.rings.rational_field import RationalField
from sage.matrix.all import *
from sage.combinat.permutation import *
from sage.combinat.set_partition import *
from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from copy import copy

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
    for x in range(dimx):
        for y in range(dimy):
            d[(x,y)]=m[x,y].value.get_generating_function()
    return matrix(RationalField(),d)

def matrix_remove_row_col(mat,row,col):
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
    dim = mat.nrows()
    P = Permutations(dim)
    S = set()
    for p in P:
        l = list()
        sgn = CombinatorialScalar([CombinatorialObject(p.signature(),p.signature())])
        for i in range(1,dim+1):
            l.append(mat[i-1,p(i)-1].value)
        cp = CartesianProduct(sgn,*l)
        for i in cp:
            weight = 1
            sign = 1
            for elm in i:
                sign = sign*elm.get_sign()
                weight = weight*elm.get_weight()
            S.add(CombinatorialObject(tuple(i),sign,weight))
    return CombinatorialScalarWrapper(CombinatorialScalar(S),parent = CombinatorialScalarRing())

def matrix_combinatorial_adjoint(mat):
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
        p_comb = CombinatorialScalar([CombinatorialObject(p,p.signature())])
        l = list()
        for i in range(1,dim+1):
            l.append(mat[p(i)-1,i-1].value)
        #This list will have empty sets, which will yield to an empty cartesian product
        #especially when the matrix input is triangular (except for the identity permutation).
        #We will now iterate through the selected entries in each column
        #and create a set of a singleton of an empty string that corresponds 
        #to the "missing" element of the tuple described in definition 39.
        for i in range(1,dim+1):
            copy_l = copy(l)
            copy_l[i-1]=CombinatorialScalar([CombinatorialObject('_',1)])
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
            l.append(CombinatorialScalarWrapper(CombinatorialScalar(L[i][j]),parent=prnt))
        M.append(l)
    return mat_space(M)

def matrix_print(mat):
    print "Printing..."
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            print "row " + str(i) + ", column " + str(j) + "; " + str(mat[i,j].value.get_size()) + " elements"
            mat[i,j].value.print_list()
            print "------------------------------"