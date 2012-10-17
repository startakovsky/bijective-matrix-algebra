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
#this was a tip from website... minimal importing

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
