import matplotlib.pyplot as plt
import random 
import numpy as np


def Bit_generator(L):
    # Generates a random series of 1 and 0 with the same chance and of length L
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



            
def Generate_MASK(L,f,int_list,t,levels,fs):
    from Filters import Filters
    t_per_bit = int(len(t)/len(int_list)) #Find the time per bit by dividing length of t by the total integers
    i1_mapped = [] #Empty array for mapped integers

    # Map the bits to the length of t
    for q in range(0, len(int_list)):    
        for i in range(0,t_per_bit):
            i1_mapped.append(int_list[q]) # For each integer append the value for length of time 

    # Edge case where there are remaining t entries just map the last bits to the remaining t 
    if len(i1_mapped) != len(t):
        remaining_t = len(t)-len(i1_mapped)
        for i in range(0,remaining_t):
            i1_mapped.append(i1_mapped[len(i1_mapped)-1])
    s = []
    i1_mapped_filtered = Filters().Low_pass(i1_mapped,2, 25, fs)
    for q in range(0,len(t)):
        s.append((i1_mapped_filtered[q])*np.cos(f*3.14*t[q]*2))

    
    return [i1_mapped, s]



def De_mod_MASK(s,ft,f,t):
    from Filters import Filters
    s1 = s*np.cos((f*3.14*t*2)+3.14/4)
    
    s2 = Filters().Low_pass(s1, f/2, 15, ft)
    return s2

    
def test():
    from Constalation_in_python import Constallation_phase
    N = 60

    f = 2

    levels = 2

    step = 1/(4*(f*2*3.14*(2**levels)))

    graph_scaling_factor = 10

    t = np.arange(0,100,step)
    import scipy as sci
    
    
    
    bits = Bit_generator(N)
    
    ints = Int_generator(bits,levels)
    
    i1_mapped,s = Generate_MASK(levels,f,ints,t,levels,1/step)
    

    
    c = np.exp(-1j*np.pi*t*f)
    
    [real_part, phase_part] = Constallation_phase(c, s)
    
    #Plots 
    fig1, (sub1, sub2,sub3) = plt.subplots(3,1, figsize=(10, 10))
    sub1.plot(t[:int(len(t)/graph_scaling_factor)], i1_mapped[:int(len(t)/graph_scaling_factor)])
    sub1.plot(t[:int(len(t)/graph_scaling_factor)], s[:int(len(t)/graph_scaling_factor)], linestyle='--')
    fft1 = np.abs(sci.fft.fft(i1_mapped))
    fft1 = fft1[:len(fft1)//2]
    fft2 = np.abs(sci.fft.fft(s))
    fft2 = fft2[:len(fft2)//2]
    fft_freq = sci.fft.fftfreq(len(t),step)
    fft_freq = fft_freq[:len(fft_freq)//2]
    sub2.plot(fft_freq,fft2, linestyle = 'dotted')
    sub2.plot(fft_freq,fft1)
    sub2.set_xlim(0,f*2**levels+1)
    sub2.set_ylim(0,10000)
    sub3.scatter(np.real(s),np.arctan2(np.imag(s),np.real(s)))
    sub3.grid(True)
    fig1.subplots_adjust(hspace=0.2)
    plt.show()
            
    
if __name__ == '__main__':
    test()