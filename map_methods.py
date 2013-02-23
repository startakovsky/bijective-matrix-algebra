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

from sage.sets.finite_set_maps import FiniteSetMaps
from sage.sets.set import Set
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper


def is_bijection(func):
	r"""
	Returns True if the function is a bijection; False otherwise.
	"""
	return func.domain().cardinality() == func.image_set().cardinality()

def is_sign_preserving(func):
	r"""
	Returns True if the function is sign preserving; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_sign() != i.get_sign():
			return False
	return True

def is_sign_reversing(func):
	r"""
	Returns True if the function is sign reversing; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_sign() != -i.get_sign():
			return False
	return True

def is_weight_preserving(func):
	r"""
	Returns True if the function is weight preserving; False otherwise.
	"""
	for i in func.domain():
		if func(i).get_weight() != i.get_weight():
			return False
	return True

			
def is_involution(func):
	r"""
	Returns True if the function is an involution; False otherwise.
	"""
	if func.domain() != func.codomain() and not(is_bijection(func)):
		return False
	else:
		for i in func.domain():
			if func(func(i)) != i:
				return False
		return True

def fixed_points(func):
    r"""
    Returns the Combinatorial Scalar of the fixed points of a map.
    """
    S = set()
    for i in func.domain():
        if func(i) == i:
            S.add(i)
    return CombinatorialScalarWrapper(S)

def not_fixed_points(func):
    r"""
    Returns the Combinatorial Scalar of the non-fixed points of a map.
    """
    return CombinatorialScalarWrapper(set(func.domain()).difference(fixed_points(func)))

def restrict_map_fixed(func):
    r"""
    Returns the same map whose domain is restricted to its fixed points.
    """
    fxd = fixed_points(func)
    M = FiniteSetMaps(fxd)
    d = dict()
    for i in fxd:
        d[i] = func(i)
    return M.from_dict(d)

def restrict_map_not_fixed(func):
    r"""
    Returns the the same map whose domain is restricted to its non-fixed points.
    """
    not_fxd = not_fixed_points(func)
    M = FiniteSetMaps(not_fxd)
    d = dict()
    for i in not_fxd:
        d[i] = func(i)
    return M.from_dict(d)

def is_SPWP(func):
    r"""
    Returns True if the function is sign preserving and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(func)):
        return False
    elif not(is_sign_preserving(func)):
        return False
    else:
        return True

def is_SRWP(func):
    r"""
    Returns True if the function is sign reversing and weight preserving; False otherwise.
    """
    if not(is_weight_preserving(func)):
        return False
    elif not(is_sign_reversing(restrict_map_not_fixed(func))):
        return False
    else:
        return True

def is_SRWP_involution(func):
    r"""
    Returns True if the function is a sign reversing, weight preserving involution; False otherwise.
    """
    return is_SRWP(func) and is_involution(func)

def is_SPWP_bijection(func):
    r"""
    Returns True if the function is a sign preserving, weight preserving bijection; False otherwise.
    """
    return is_SPWP and is_bijection(func)
            
def inverse(func):
    dic = func.fibers()
    for i in func.codomain():
        dic[i] = set(dic[i]).pop()
    return FiniteSetMaps(func.codomain(),func.domain()).from_dict(dic)
