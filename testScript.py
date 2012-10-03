from sage.all import *



"""
Test Script
"""


C = CombinatorialScalarElement("Rock",1,[0,1,1,2]);C  
D = CombinatorialScalarElement("Paper",-1,[1,2]);D    
E = CombinatorialScalarElement("Scissors",1,[1,2,0]);E
F = CombinatorialScalarElement("Lizard",1,[1,2,0]);F  
G = CombinatorialScalar((C,D,E,F)); G                 

C.get_sign()                                          
C.get_name()                                          
C.get_weight()                                        
C.get_genfunc()                                       

G.get_generating_function()                           
G.get_sign_function()                                 
G.get_weight_function()                               
G.get_size()                                        
G.is_fully_cancelled()
G.print_list()


C1 = CombinatorialScalarElement("Rock",1,[0,1,1,2]);C1 
D1 = CombinatorialScalarElement("Paper",-1,[1,2]);D1  
E1 = CombinatorialScalarElement("Scissors",1,[1,2,0]);E1
F1 = CombinatorialScalarElement("Lizard",1,[1,2,0]);F1  
G1 = CombinatorialScalar((C1,D1,E1,F1)); G1

M = FiniteSetMaps(G,G1)
f = M.random_element()