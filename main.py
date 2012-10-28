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
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar


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

def fixed_points(f):
    """
    Returns the Combinatorial Scalar of the fixed points of a map.
    """
    S = set()
    for i in f.domain():
        if f(i) == i:
            S.add(i)
    return CombinatorialScalar(S)

def not_fixed_points(f):
    """
    Returns the Combinatorial Scalar of the non-fixed points of a map.
    """
    return CombinatorialScalar(set(f.domain()).difference(fixed_points(f)))

def restrict_map_fixed(f):
    """
    Returns the same map whose domain is restricted to its fixed points.
    """
    fxd = fixed_points(f)
    M = FiniteSetMaps(fxd)
    d = dict()
    for i in fxd:
        d[i] = f(i)
    return M.from_dict(d)

def restrict_map_not_fixed(f):
    """
    Returns the the same map whose domain is restricted to its non-fixed points.
    """
    not_fxd = not_fixed_points(f)
    M = FiniteSetMaps(fxd)
    d = dict()
    for i in not_fxd:
        d[i] = f(i)
    return M.from_dict(d)

def is_SPWP(f):
    """
    Returns True if the function is sign preserving and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(f)):
        return False
    elif not(is_sign_preserving(f)):
        return False
    else:
        return True

def is_SRWP(f):
    """
    Returns True if the function is sign reversing and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(f)):
        return False
    elif not(is_sign_reversing(restrict_map_not_fixed(f))):
        return False
    else:
        return True

def is_SRWP_involution(f):
    """
    Returns True if the function is a sign reversing, weight preserving involution; False otherwise.
    """
    return is_SRWP(f) and is_involution(f)

def is_SPWP_bijection(f):
    """
    Returns True if the function is a sign preserving, weight preserving bijection; False otherwise.
    """
    return is_SPWP and is_bijection(f)

def involution_dict(mat):
    """
    Returns a dictionary of arbitrary involutions on the entries of a Combinatorial Matrix.
    """
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
    """
    Prints out all involutions in a dict of involutions (or other maps).
    """
    for key in sorted(f):
        print "---------------------------------------"
        print "row " + str(key[0]) + ", column " + str(key[1]) + ": " + str(len(f[key].domain())) + " elements."
        for elm in f[key].domain():
            print str(elm) + " --> " + str(f[key](elm))