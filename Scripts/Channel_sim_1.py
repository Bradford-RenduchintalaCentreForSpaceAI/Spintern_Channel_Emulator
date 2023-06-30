## Imports
from Mary_ASK import Int_generator, Bit_generator, Generate_MASK
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy as sci
from Filters import Filters

# Basics Var declaration 

N = 100

f = 10

Levels = 3

step = 1/(3*(f*2*3.14))

t = np.arange(0,10,step)

graph_scaling_factor = 10


#AWGN Noise 

def AWGN_Noise(s,ratio):
    out_signal = []
    for q in range(0,len(s)):
        out_signal.append(s[q] + ratio*(random.gauss(mu=0.0, sigma=1.0)))
    return out_signal





# Transmitter 

bits = Bit_generator(N)

ints = Int_generator(bits, Levels)

[bits, TxSignal] = Generate_MASK(N, f, ints, t, Levels)


# Channel 

RxSignal = AWGN_Noise(TxSignal, 0)



# Reciever 

for i in range(len(RxSignal)):
    RxSignal[i] = RxSignal[i]-np.cos(f*3.14*t[i])

#RxSignal = Filters().Low_pass(s = RxSignal, cutoff=10, N = 10, fs = 1/step)








#Plots 
fig1, (sub1, sub2,sub3) = plt.subplots(3,1, figsize=(10, 10))
sub1.plot(t[:int(len(t)/graph_scaling_factor)], TxSignal[:int(len(t)/graph_scaling_factor)])
sub1.plot(t[:int(len(t)/graph_scaling_factor)], RxSignal[:int(len(t)/graph_scaling_factor)], linestyle='--')
sub1.set_ylim(-20, 20)
sub2.plot(abs(sci.fft.fft(TxSignal)))
sub2.plot(abs(sci.fft.fft(RxSignal)))
sub3.scatter(np.real(bits),np.imag(bits))
sub3.grid(True)
fig1.subplots_adjust(hspace=0.2)
plt.show()

