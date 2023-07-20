import matplotlib.pyplot as plt
import random 
import numpy as np
N = 10000
L = N/8
step = 0.01
wt1 = (1*10^9)*np.pi
wt2 = (2*10^9)*np.pi
graph_scaling_factor = 5
#Create time and but sequance 
t = np.arange(0, N, step)
b1 = random.choices([0,1], weights = (50,50), k = round(L)) 

ratio = 0
# Map the bit values to the time sequance
c1 = 10*np.sin(t*wt1)
c2 = 10*np.sin(t*wt2)
b1_mapped = []
for q in range(0,round(L)):
    for i in range(0,round((N/L)/step)):
        b1_mapped.append(b1[q])
# Generate Noise 

n = [ratio*random.gauss(mu=0.0, sigma=1.0) for i in t]


# Define the two signals
s = []
for i in range(0,len(t)):
    if b1_mapped[i] == 1:
        s.append(c1[i]+n[i])
    if b1_mapped[i] == 0:
        s.append(c2[i]+n[i])
        





fig1, (sub1, sub2) = plt.subplots(1, 2, figsize=(10, 10))
sub1.plot([t[i] for i in range(1,round(N/graph_scaling_factor))],[b1_mapped[i] for i in range(1,round(N/graph_scaling_factor))])
sub1.plot([t[i] for i in range(1,round(N/graph_scaling_factor))],[s[i] for i in range(1,round(N/graph_scaling_factor))], linestyle = 'dashed')
# sub2.plot(np.fft.fft(b1_mapped))
sub2.plot(np.fft.fft(s))
# sub2.plot(np.fft.fft(c1))

# sub2.legend(["Bit train","Signal", "Carrier"])
# sub2.set(ylim=(0, 1*10**5))

plt.show()

