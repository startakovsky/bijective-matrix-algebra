from sage.all import *



"""
Test Script
"""


sage: C = CombinatorialScalarElement("Rock",1,[0,1,1,2]);C  
sage: D = CombinatorialScalarElement("Paper",-1,[1,2]);D    

sage: E = CombinatorialScalarElement("Scissors",1,[1,2,0]);E
sage: F = CombinatorialScalarElement("Lizard",1,[1,2,0]);F  
sage: G = CombinatorialScalar((C,D,E,F)); G                 

sage: C.get_sign()                                          

sage: C.get_name()                                          

sage: C.get_weight()                                        

sage: C.get_genfunc()                                       

sage: G.get_generating_function()                           

sage: G.get_sign_function()                                 
sage: G.get_weight_function()                               

sage: G.get_size()                                        

sage: G.is_fully_cancelled()

sage: G.print_list()

