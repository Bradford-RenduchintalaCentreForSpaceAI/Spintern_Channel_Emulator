## Imports
from Mary_ASK import Int_generator, Bit_generator, Generate_MASK, De_mod_MASK
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy as sci
from Filters import Filters

# Basics Var declaration 

N = 10

f = 100

Levels = 3

step = 1/(3*(f*2*3.14))

t = np.arange(0,100,step)

graph_scaling_factor = 1


#AWGN Noise 

def AWGN_Noise(s,ratio):
    out_signal = []
    for q in range(0,len(s)):
        out_signal.append(s[q] + ratio*(random.gauss(mu=0.0, sigma=1.0)))
    return out_signal





# Transmitter 

bits = Bit_generator(N)

ints = Int_generator(bits, Levels)

[ints_mapped, TxSignal] = Generate_MASK(N, f, ints, t, Levels)


# Channel 

RxSignal = AWGN_Noise(TxSignal, 0.5)

# More to add here 

# Reciever 

RxInt = De_mod_MASK(RxSignal, 1/step, f, t)


#Plots 
fig1, (sub1, sub2,sub3) = plt.subplots(3,1, figsize=(10, 10))
sub1.plot(t[:int(len(t)/graph_scaling_factor)], ints_mapped[:int(len(t)/graph_scaling_factor)])
sub1.plot(t[:int(len(t)/graph_scaling_factor)], RxInt[:int(len(t)/graph_scaling_factor)], linestyle='--')
sub2.plot(abs(sci.fft.fft(TxSignal)))
sub2.plot(abs(sci.fft.fft(RxSignal)))
sub3.scatter(np.real(ints_mapped),np.imag(ints_mapped))
sub3.scatter(np.real(RxInt),np.imag(RxInt))
sub3.grid(True)
fig1.subplots_adjust(hspace=0.2)
plt.show()

