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

from sage.all import *


def create_variables(n):
	"""Given an integer n, this method returns the string 'x1,..,xn' """
	variables = ''
	for i in range(n):
		variables+= 'x' + str(i+1) + ','
	variables = variables[0:len(variables)-1]
	return variables

def assign_weight_monomial(l):
	"""Returns weight monomial to list l."""
	v = var(create_variables(len(l)))
	monomial = 1
	for i in range(len(l)):
		monomial = monomial * v[i]**l[i]
	return monomial
