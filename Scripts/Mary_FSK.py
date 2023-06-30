import matplotlib.pyplot as plt
import random 
import numpy as np
N = 60

f = float(1)

levels = 2

step = 1/(4*(f*2*3.14*(2**levels)))

graph_scaling_factor = 250

t = np.arange(0,100,step)


def Bit_generator(L):
    bits = random.choices([0, 1], weights = [50,50], k = N)
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
    
            
def Generate_MFSK(L,f,int_list,t,levels):
    w = [f*2*3.1425*i for i in range(1,(2**levels)+1)]
    print(w)
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
    for i in range(0,len(t)):
        w_relative = w[i1_mapped[i]]
        s.append(np.cos(t[i]*w_relative))
    return [i1_mapped,s]
    

bits = Bit_generator(N)

ints = Int_generator(bits, levels)

MPSK = Generate_MFSK(N, f, ints, t, levels)


fig1, (sub1, sub2) = plt.subplots(1, 2, figsize=(10, 10))
sub1.plot([t[i] for i in range(0,graph_scaling_factor)], [MPSK[1][i] for i in range(0,graph_scaling_factor)])
sub2.plot(abs(np.fft.fft(MPSK[1])))
plt.show()
