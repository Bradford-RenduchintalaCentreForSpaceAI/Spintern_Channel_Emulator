import matplotlib.pyplot as plt
import numpy as np



 
def Constallation_phase(carrier,signal):
    import numpy as np
    
    carrier = np.array(carrier)
    signal = np.array(signal)
    real_part = (np.sqrt((np.real(carrier-signal)**2)+(np.imag(carrier-signal)**2)))

    phase_part = (np.arctan2(np.imag(carrier),np.real(carrier)) - np.arctan2(np.imag(signal),np.real(signal)))
    
    
    return [real_part*np.cos(phase_part),real_part*np.sin(phase_part)]           
           
     
           
     
        
     
def test():
    pi = np.pi
    phase = 2*pi/3
    t = np.arange(0,2*pi,0.01)

    c = np.exp(-1j*pi*t)
    
    s1 = []
    for i in range(len(t[:len(t)//2])):
        s1.append(np.exp(-1j*pi*t[i]-(1j*phase)))
        
    for i in range(len(t[len(t)//2:])):
        s1.append(np.exp(-1j*pi*t[i]-(1j*phase/2)))
    
    s1 = np.exp(-1j*pi*t-(1j*phase))  
    [real_part,phase_part] = Constallation_phase(c, s1)
    
    
    plt.figure()
    plt.plot(t,s1)
    plt.plot(t,c)
    plt.show()
    
    plt.figure(figsize = (10, 10))
    plt.title("BPSK Constelation diagram")
    plt.scatter(phase_part,real_part);
    plt.axhline(0, color = "black", linewidth = 1)
    plt.axvline(0, color = "black", linewidth = 1)
    plt.xlim((-4,4))
    plt.xlabel("Real")
    plt.ylabel("Phase")
    plt.ylim((-4,4))
    
    plt.grid(True)
    plt.show()
    
    
    
    
if __name__ == '__main__':
    test()