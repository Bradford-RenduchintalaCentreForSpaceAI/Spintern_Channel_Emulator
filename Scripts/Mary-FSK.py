import matplotlib.pyplot as plt
import random 
import numpy as np
N = 100
L = round(N/8)
f = float(10)
levels = 4
step = 1/(2.5*(f*2*3.14*(2**levels)))
graph_scaling_factor = 160000




def Bit_generator(L):
    bits = random.choices([0, 1], weights = [50,50], k = L)
    return bits
    


def Int_generator(bits,Levels):
    
    #Split the bits into levels
    cnt = 0
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
    
            
    
    
    return int_list
def Generate_MFSK(L,f,int_list,t,levels):
    w = [f*2*3.1425*i for i in range(1,(2**levels)+1)]

    t_per_bit = int(len(t)/len(int_list))

    i1_mapped = []

    possible_levels = [i for i in range(0,2**levels)]
    
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
    for i in range(0,len(t)):
        w_relative = w[i1_mapped[i]]
        s.append(np.cos(t[i]*w_relative))
    return [i1_mapped,s]
    

bits = Bit_generator(L)
in_list =Int_generator(bits, levels)
shit = Generate_MFSK(L,f,in_list,t,levels)[0]




fig1, ax1 = plt.subplots()
ax1.plot(shit[0])
#ax1.plot(t,shit[1])
plt.show()


