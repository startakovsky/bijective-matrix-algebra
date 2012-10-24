r"""
This file contains the methods used in the Bijective Matrix Algebra package.

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

from sage.matrix.all import *
from sage.rings.rational_field import RationalField
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.sets.set import Set
from sage.bijectivematrixalgebra.stirling_matrices import matrix_generating_function


def is_bijection(f):
	"""
	Returns True if the function is a bijection; False otherwise.
	"""
	return f.domain().cardinality() == f.image_set().cardinality()

def is_sign_preserving(f):
	"""
	Returns True if the function is sign preserving; False otherwise.
	"""
	for i in f.domain():
		if f(i).get_sign() != i.get_sign():
			return False
	return True

def is_sign_reversing(f):
	"""
	Returns True if the function is sign reversing; False otherwise.
	"""
	for i in f.domain():
		if f(i).get_sign() != -i.get_sign():
			return False
	return True

def is_weight_preserving(f):
	"""
	Returns True if the function is weight preserving; False otherwise.
	"""
	for i in f.domain():
		if f(i).get_weight() != i.get_weight():
			return False
	return True

			
def is_involution(f):
	"""
	Returns True if the function is an involution; False otherwise.
	"""
	if f.domain() != f.codomain():
		return False
	else:
		for i in f.domain():
			if f(f(i)) != i:
				return False
		return True

def involution_dict(mat):
    if matrix_generating_function(mat)!=MatrixSpace(RationalField(),mat.nrows(),mat.ncols()).identity_matrix():
        print "ERROR: Input needs to be equal to the identity."
    else:
        f = dict()
        for x in range(mat.nrows()):
            for y in range(mat.ncols()):
                if x <> y:
                    f[(x,y)] = mat[x,y].value.create_involution()
                else:
                    t = mat[x,y].value
                    for i in t:
                        _M = FiniteSetMaps(t)
                        f[(x,y)] = _M.from_dict({i:i})
        return f
        
def print_involution_dict(f):
    for key in sorted(f):
        print "---------------------------------------"
        print "row " + str(key[0]) + ", column " + str(key[1]) + ": " + str(len(f[key].domain())) + " elements."
        for elm in f[key].domain():
            print str(elm) + " --> " + str(f[key](elm))