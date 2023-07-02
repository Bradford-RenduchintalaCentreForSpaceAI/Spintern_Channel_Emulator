import matplotlib.pyplot as plt
import random 
import numpy as np


def Bit_generator(L):
    bits = random.choices([0, 1], weights = [50,50], k = L)
    return bits
    


def Int_generator(bits,Levels):
    
    #Split the bits into levels
    int_list = []
    int_value = 0 
    for i in range(0,round(len(bits)/Levels)):
        for q in range(0,Levels):
            
            if i*Levels+q >= len(bits):
                break
            int_value += bits[i*Levels+q]*2**q
        int_list.append(int_value)
        int_value = 0 
    return int_list



            
def Generate_MASK(L,f,int_list,t,levels):
    from Filters import Filters
    t_per_bit = int(len(t)/len(int_list))
    i1_mapped = []

    
    # Map the bits to the length of t
    for q in range(0, len(int_list)):    
        for i in range(0,t_per_bit):
            i1_mapped.append(int_list[q])

    # Edge case where there are remaining t entries just map the last bits to the remaining t 
    if len(i1_mapped) != len(t):
        remaining_t = len(t)-len(i1_mapped)
        for i in range(0,remaining_t):
            i1_mapped.append(i1_mapped[len(i1_mapped)-1])
    s = []
    for q in range(0,len(t)):
        s.append((i1_mapped[q])*np.cos(f*3.14*t[q]))
    
    
    return [i1_mapped, s]



def De_mod_MASK(s,ft,f,t):
    from Filters import Filters
    s1 = s*np.cos((f*3.14*t)+3.14/4)
    
    s2 = Filters().Low_pass(s1, f/2, 15, ft)
    return s2

    
        
