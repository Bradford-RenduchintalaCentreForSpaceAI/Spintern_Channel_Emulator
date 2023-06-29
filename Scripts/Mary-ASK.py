import matplotlib.pyplot as plt
import random 
import numpy as np
N = 50

f = 100
Levels = 4

step = 1/(3*(f*2*3.14*(2**Levels)))

t = np.arange(0,1000,step)

graph_scaling_factor = len(t)

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



            
def Generate_MASK(L,f,int_list,t,levels):
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
        s.append((i1_mapped[q]+1)*np.cos(f*3.14*t[q]))
    
    
    return [i1_mapped, s]


def AWGN_Noise(s,ratio):
    for q in range(0,len(s)):
        s[q] = s[q] + ratio*(random.gauss(mu=0.0, sigma=1.0))
    return s

bits = Bit_generator(N)

ints = Int_generator(bits, Levels)

MASK = Generate_MASK(N, f, ints, t, Levels)

MASK[1] = AWGN_Noise(MASK[1], 1)



fig1, (sub1, sub2) = plt.subplots(1, 2, figsize=(10, 10))
sub1.plot([t[i] for i in range(0,graph_scaling_factor)], [MASK[1][i] for i in range(0,graph_scaling_factor)])
sub2.plot(abs(np.fft.fft(MASK[1])))
plt.show()

        
