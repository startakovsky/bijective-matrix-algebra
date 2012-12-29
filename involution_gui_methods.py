r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It generates the methods needed in order to generate one's own involutions
for any matrix which reduces to I and is equivalent to a singleton along
each entry in its diagonal.

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
from sage.bijectivematrixalgebra.reduction_methods import *
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.bijectivematrixalgebra.reduction_maps_dicts import ReductionMapsDict
from sage.bijectivematrixalgebra.reduction_maps import ReductionMaps

def tab_list(y, headers = None):
    '''
    Converts a list into an html table with borders.
    '''
    s = '<table border = 1>'
    if headers:
        for q in headers:
            s = s + '<th>' + str(q) + '</th>'
    for x in y:
        s = s + '<tr>'
        for q in x:
            s = s + '<td>' + str(q) + '</td>'
        s = s + '</tr>'
    s = s + '</table>'
    return s

def total_maps(mat):
    nrows = mat.nrows()
    ncols = mat.ncols()
    total = 0
    for i in range(nrows):
        for j in range(ncols):
            total += mat[i,j].get_size()
    return total
    
def reduction_identity(mat,involution_dict, st = None):
    r"""
    When a matrix reduces to the identity, this returns
    a ReductionMapDict of from a matrix to I.
    """
    dim = mat.nrows()
    fs = involution_dict
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
