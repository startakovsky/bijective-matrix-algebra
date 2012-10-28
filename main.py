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


def is_bijection(func):
	"""
	Returns True if the function is a bijection; False otherwise.
	"""
	return func.domain().cardinality() == func.image_set().cardinality()

def is_sign_preserving(func):
	"""
	Returns True if the function is sign preserving; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_sign() != i.get_sign():
			return False
	return True

def is_sign_reversing(func):
	"""
	Returns True if the function is sign reversing; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_sign() != -i.get_sign():
			return False
	return True

def is_weight_preserving(func):
	"""
	Returns True if the function is weight preserving; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_weight() != i.get_weight():
			return False
	return True

			
def is_involution(func):
	"""
	Returns True if the function is an involution; False otherwise.
	"""
	if func.domain() != func.codomain():
		return False
	else:
		for i in func.domain():
			if func(func(i)) != i:
				return False
		return True

def fixed_points(func):
    """
    Returns the Combinatorial Scalar of the fixed points of a map.
    """
    S = set()
    for i in func.domain():
        if func(i) == i:
            S.add(i)
    return CombinatorialScalar(S)

def not_fixed_points(func):
    """
    Returns the Combinatorial Scalar of the non-fixed points of a map.
    """
    return CombinatorialScalar(set(func.domain()).difference(fixed_points(func)))

def restrict_map_fixed(func):
    """
    Returns the same map whose domain is restricted to its fixed points.
    """
    fxd = fixed_points(func)
    M = FiniteSetMaps(fxd)
    d = dict()
    for i in fxd:
        d[i] = func(i)
    return M.from_dict(d)

def restrict_map_not_fixed(func):
    """
    Returns the the same map whose domain is restricted to its non-fixed points.
    """
    not_fxd = not_fixed_points(func)
    M = FiniteSetMaps(not_fxd)
    d = dict()
    for i in not_fxd:
        d[i] = func(i)
    return M.from_dict(d)

def is_SPWP(func):
    """
    Returns True if the function is sign preserving and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(func)):
        return False
    elif not(is_sign_preserving(func)):
        return False
    else:
        return True

def is_SRWP(func):
    """
    Returns True if the function is sign reversing and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(func)):
        return False
    elif not(is_sign_reversing(restrict_map_not_fixed(func))):
        return False
    elif not(is_sign_preserving(restrict_map_fixed(func))):
        return False
    else:
        return True

def is_SRWP_involution(func):
    """
    Returns True if the function is a sign reversing, weight preserving involution; False otherwise.
    """
    return is_SRWP(func) and is_involution(func)

def is_SPWP_bijection(func):
    """
    Returns True if the function is a sign preserving, weight preserving bijection; False otherwise.
    """
    return is_SPWP and is_bijection(func)

def involution_dict(mat):
    """
    Returns a dictionary of arbitrary involutions on the entries of a Combinatorial Matrix.
    Note all weights must be 1 for this.  It may be extended upon in the future.
    """
    if matrix_generating_function(mat)!=MatrixSpace(RationalField(),mat.nrows(),mat.ncols()).identity_matrix():
        print "ERROR: Input needs to be equal to the identity."
    else:
        func = dict()
        for x in range(mat.nrows()):
            for y in range(mat.ncols()):
                if x <> y:
                    func[(x,y)] = mat[x,y].value.create_involution()
                else:
                    t = mat[x,y].value
                    for i in t:
                        _M = FiniteSetMaps(t)
                        func[(x,y)] = _M.from_dict({i:i})
        return func
        
def print_involution_dict(func):
    """
    Prints out all involutions from a dictionary of involutions indexed by matrix entries (or other maps).
    """
    for key in sorted(func):
        print "---------------------------------------"
        print "row " + str(key[0]) + ", column " + str(key[1]) + ": " + str(len(func[key].domain())) + " elements."
        for elm in func[key].domain():
            print str(elm) + " --> " + str(func[key](elm))