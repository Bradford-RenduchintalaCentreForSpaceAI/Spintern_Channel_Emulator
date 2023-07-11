## Imports
from .mary_ask import Int_generator, Bit_generator, Generate_MASK, De_mod_MASK
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy as sci
from .filters import low_pass
import time

# Basics Var declaration

N = 50

f = 200

Levels = 3

step = 1/(3*(f*2*3.14))

t = np.arange(0,1000,step)

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

[ints_mapped, TxSignal] = Generate_MASK(N, f, ints, t, Levels,1/step)


# Channel

RxSignal = AWGN_Noise(TxSignal, 0.25)



# Reciever

RxInt = De_mod_MASK(RxSignal, 1/step, f, t)


RxIntSampled = []
cnt = 0
for i in range(len(t)):
    if cnt == len(t)/100:
        RxIntSampled.append(RxInt[i])
        cnt = 0
    cnt += 1




#Plots
fig1, (sub1, sub2,sub3,sub4) = plt.subplots(4,1, figsize=(10, 10))
sub1.plot(t[:int(len(t)/graph_scaling_factor)], ints_mapped[:int(len(t)/graph_scaling_factor)])
sub1.plot(t[:int(len(t)/graph_scaling_factor)], RxInt[:int(len(t)/graph_scaling_factor)], linestyle='--')
fft1 = np.abs(sci.fft.fft(TxSignal))
fft1 = fft1[:len(fft1)//2]
fft2 = np.abs(sci.fft.fft(ints_mapped))
fft2 = fft2[:len(fft2)//2]
fft_freq = sci.fft.fftfreq(len(t),step)
fft_freq = fft_freq[:len(fft_freq)//2]
sub2.plot(fft_freq,fft2, linestyle = 'dotted')
sub2.plot(fft_freq,fft1)
sub2.set_xlim(0,f*4)
sub2.set_ylim(0,10000)
sub3.scatter(np.real(ints_mapped),np.imag(ints_mapped))
sub4.scatter(np.real(RxIntSampled),np.imag(RxIntSampled))
sub3.grid(True)
sub4.grid(True)
fig1.subplots_adjust(hspace=0.2)
plt.show()




# for i in range(len(t)):
#     print(f"t:{t[i]} Ints:{ints_mapped[i]} RxInts:{RxInt[i]}")
#     time.sleep(0.001)

