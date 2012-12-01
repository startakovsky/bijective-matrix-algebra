r"""

This file contains the methods used in the Bijective Matrix Algebra package.

It generates the LoehrMendes bijection as stated in Theorem 45 of the paper
by Loehr and Mendes.

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

from sage.structure.all import SageObject
from sage.bijectivematrixalgebra.matrix_methods import *
from sage.bijectivematrixalgebra.reduction_methods import *

class LoehrMendes(SageObject):
    r"""
    This class represents the process as outlined in Theorem 47
    """
    
    def __init__(self,A,B,red_AB_to_I,repr=None):
        if repr == None:
            self._repr = "description missing"
        else:
            self_repr = repr
        _dim = A.nrows()
        self._reduction_AB_to_I = red_AB_to_I
        self._A = A
        self._B = B
        self._AB = matrix_multiply(self._A,self._B)
        self._BA = matrix_multiply(self._B,self._A)
        self._adj_A = matrix_combinatorial_adjoint(self._A)
        _I = identity_matrix(_dim)
        _adj_AA = matrix_multiply(self._adj_A,self._A)
        _det_A = matrix_determinant(self._A)
        _det_AI = matrix_identity_multiply_scalar(_det_A,_dim,_dim)
        _reduction_45 = reduction_identity_matrix(_det_AI,"detAI to I or reduction_45")
        _reduction_adj_AA_to_detAI = reduction_lemma_40(_adj_AA)
        _reduction_adj_AA_to_I = _reduction_adj_AA_to_detAI.transitive(_reduction_45)
        #the setup complete

        _reduction_AB_to_I = red_AB_to_I
        
        _reduction_12 = reduction_matrix_ABCD_to_ApBCpD(self._adj_A,self._A,self._B,self._A,"reduction_12")
        _mat2 = _reduction_12.get_matrix_B()
        _reduction_23 = reduction_lemma_28_23(_mat2, self._reduction_AB_to_I)
        _reduction_13 = _reduction_12.transitive(_reduction_23)
        _mat3 = _reduction_23.get_matrix_B()
        #reduction_13 complete...
        
        
        _reduction_34 = reduction_matrix_AIB_AB(_mat3).transitive(_reduction_adj_AA_to_detAI,"reduction_34")
        _reduction_35 = _reduction_34.transitive(_reduction_45,"reduction_35")
        #reduction_35 complete...
        
        self._reduction_15 = _reduction_13.transitive(_reduction_35,"reduction_15")
        #reduction_15 and its dependencies complete
        
        
        _reduction_16 = reduction_matrix_ABCD_to_pABpCD(self._adj_A,self._A,self._B,self._A,"reduction_16")
        _mat6 = _reduction_16.get_matrix_B()
        _reduction_68 = reduction_lemma_28_68(_mat6,_reduction_adj_AA_to_I)
        _reduction_18 = _reduction_16.transitive(_reduction_68)
        _mat8= _reduction_68.get_matrix_B()
        #reduction_18 complete...
        
        _reduction_89 = reduction_matrix_IAB_AB(_mat8)
        self._reduction_19 = _reduction_18.transitive(_reduction_89,"reduction_19")
        #reduction_19 complete
        
        self._reduction_LoehrMendes = self._reduction_15.confluence(self._reduction_19, "LoehrMendes")
        #confluence lemma complete#
        
    def __repr__(self):
        return "The Loehr-Mendes Bijection: " + self._repr

    def __eq__(self,other):
        return self.get_confluence_reduction() == other.get_confluence_reduction()
    
    def get_reduction15(self):
        return self._reduction_15

    def get_reduction19(self):
        return self._reduction_19
    
    def get_confluence_reduction(self):
        return self._reduction_LoehrMendes
    
    def get_original_reduction(self):
        return self._reduction_AB_to_I

    def get_reduction(self):
        return self._reduction_LoehrMendes
    
    def get_original_left(self):
        return self._A
    
    def get_original_right(self):
        return self._B
    
    def get_leftright(self):
        return self._AB
    
    def get_rightleft(self):
        return self._BA
    
    def get_adjleft(self):
        return self._adj_A
