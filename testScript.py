from sage.bijectivematrixalgebra.combinatorial_objects import CombinatorialObject
from sage.bijectivematrixalgebra.combinatorial_scalars import CombinatorialScalar
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarRing
from sage.bijectivematrixalgebra.combinatorial_scalar_rings_and_elements import CombinatorialScalarWrapper



"""
Test Script
"""


Cg = CombinatorialObject("Rock",1);Cg 
Dg = CombinatorialObject("Paper",-1,[1]);Dg    
Eg = CombinatorialObject("Scissors",-1,[2]);Eg
Fg = CombinatorialObject("Lizard",1);Fg  
g = CombinatorialScalar((Cg,Dg,Eg,Fg)); g                


Cg.get_sign()                                          
Cg.get_object()                                          
Cg.get_weight()                                        
Cg.get_genfunc()
Cg.get_detail()                                       


g.get_generating_function()                           
g.get_sign_function()                                 
g.get_weight_function()                               
g.get_size()                                        
g.is_fully_cancelled()
g.print_list()

Ch = CombinatorialObject("Corn",1,[0,0,0,0,0,0,2]);Ch 
Dh = CombinatorialObject("Garlic",-1,[1,2]);Dh
Eh = CombinatorialObject("Cucumbers",1,[1,2,2]);Eh
Fh = CombinatorialObject("Blizzard",1,[0]);Fh  
h = CombinatorialScalar((Ch,Dh,Eh,Fh)); h

R = CombinatorialScalarRing()
G = CombinatorialScalarWrapper(g, parent=R)
H = CombinatorialScalarWrapper(h,parent=R)

G+H
print 'G*H'
print G*H

n = 5;

