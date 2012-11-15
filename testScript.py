from sage.bijectivematrixalgebra.all import *




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

R = CombinatorialScalarRing()
G = CombinatorialScalarWrapper(g)
H = CombinatorialScalarWrapper([Ch,Dh,Eh,Fh])

G+H
print 'G*H'
print G*H

dim = 4

m1 = Stirling1Matrix(dim)
m2 = Stirling2Matrix(dim)
Ss = matrix_multiply(m2,m1)
f = bma.involution_dict(Ss)
bma.print_dict_of_maps(f)
