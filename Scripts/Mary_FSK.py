import matplotlib.pyplot as plt
import random 
import numpy as np
N = 60

f = 2

levels = 2

step = 1/(4*(f*2*3.14*(2**levels)))

graph_scaling_factor = 20

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
    
            
def Generate_MFSK(L,f,int_list,t,levels,fs):
    from Filters import Filters
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
 


   
def test():
    import scipy as sci
    
    bits = Bit_generator(N)
    
    ints = Int_generator(bits, levels)
    
    [ints_mapped, TxSignal] = Generate_MFSK(N, f, ints, t, levels ,1/step)
    
    
    #Plots 
    fig1, (sub1, sub2,sub3,sub4) = plt.subplots(4,1, figsize=(10, 10))
    sub1.plot(t[:int(len(t)/graph_scaling_factor)], ints_mapped[:int(len(t)/graph_scaling_factor)])
    sub1.plot(t[:int(len(t)/graph_scaling_factor)], TxSignal[:int(len(t)/graph_scaling_factor)], linestyle='--')
    fft1 = np.abs(sci.fft.fft(TxSignal))
    fft1 = fft1[:len(fft1)//2]
    fft2 = np.abs(sci.fft.fft(ints_mapped))
    fft2 = fft2[:len(fft2)//2]
    fft_freq = sci.fft.fftfreq(len(t),step)
    fft_freq = fft_freq[:len(fft_freq)//2]
    sub2.plot(fft_freq,fft2, linestyle = 'dotted')
    sub2.plot(fft_freq,fft1)
    sub2.set_xlim(0,f*2**levels+1)
    sub2.set_ylim(0,10000)
    sub3.scatter(np.real(ints_mapped),np.imag(ints_mapped))
    sub4.scatter(np.real(ints_mapped),np.imag(ints_mapped))
    sub3.grid(True)
    sub4.grid(True)
    fig1.subplots_adjust(hspace=0.2)
    plt.show()


test()