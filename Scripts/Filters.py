import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


class Filters():
    from scipy import signal
    
    
    def Low_pass(self,s,cutoff,N,fs):
        sos = signal.butter(N, cutoff, fs = fs,btype ='low', output = 'sos')
        return signal.sosfilt(sos, s)

        
        
    
    
    
    
    
    
def test(f):
    t =np.linspace(0,100,100)
    signal = np.sin(t*f)+np.sin(t)
    new_s = Filters().Low_pass(signal, 10,10,100)
    print(new_s,signal)
    plt.plot(t,new_s)
    plt.plot(t,signal, linestyle = '--')
    
