import matplotlib.pyplot as plt
import random 
import numpy as np
def Bit_generator(N):
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



def MPSK_Generator(bits,t,L,f):
    theta = [(2*np.pi*n)/(2**L) for n in range(0,2**L)]
    print(theta)

    t_per_bit = int(len(t)/len(bits))
    i1_mapped= []
    for i in range(0,len(bits)):
        for ii in range(0,t_per_bit):
            i1_mapped.append(bits[i])
            
    if len(i1_mapped) != len(t):
        remaining_t = len(t)-len(i1_mapped)
        for i in range(0,remaining_t):
            i1_mapped.append(i1_mapped[len(i1_mapped)-1])
    s = []      
    for i in range(len(t)):
        s.append(np.exp(-1j*2*np.pi*t[i]*f+1j*theta[i1_mapped[i]]))
        #print(theta[i1_mapped[i]])
            
    
        
    
    return [i1_mapped,s]
    

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
    
    i1_mapped,s = MPSK_Generator(ints, t,levels,f)
    

    
    c = np.exp(-1j*2*np.pi*t*f)
    
    [real_part, phase_part] = Constallation_phase(c, s)
    
    #Plots 
    fig1, (sub1, sub2,sub3,sub4) = plt.subplots(4,1, figsize=(10, 10))
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
    sub3.scatter(np.real(s),np.arctan(np.imag(s)/np.real(s)))
    sub4.scatter(real_part,phase_part)
    sub3.grid(True)
    sub4.grid(True)
    sub4.set_xlim((-4,4))
    sub4.set_ylim((-4,4))
    fig1.subplots_adjust(hspace=0.2)
    plt.show()
    
    
    
    
if __name__ == '__main__':
    test()