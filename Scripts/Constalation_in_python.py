
import matplotlib.pyplot as plt
import numpy as np
import random


 
def Constallation_phase(carrier,signal):
    


    real_part = np.sqrt(np.real(signal)**2+np.imag(signal)**2)

    phase_part = np.arctan(np.imag(signal)/np.real(signal)) - np.arctan(np.imag(carrier)/np.real(carrier))


    for i in range(len(phase_part)):
        if phase_part[i] < 0:
           phase_part[i] = np.pi-np.abs(phase_part[i])
    
    return [real_part,phase_part]           
           
     
           
     
        
     
def test():
    
    pi = np.pi
    phase = pi/2
    t = np.arange(0,2*pi,0.01 )
    s = np.cos(2*pi*t+phase)
    
    c = np.exp(-1j*2*pi*t)
    
    s1 = np.exp(-1j*2*pi*t-1j*phase)  
    [real_part,phase_part] = Constallation_phase(c, s1)
    fig = plt.figure()
    
    
    plt.scatter(real_part,phase_part)
    
    plt.xlim((-4,4))
    
    plt.ylim((-4,4))
    
    plt.grid(True)
    plt.show()
    
    print(phase_part)
    
    
    
if __name__ == '__main__':
    test()