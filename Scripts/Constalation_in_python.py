import matplotlib.pyplot as plt
import random 
import numpy as np
N = 10000
L = N/8
step = 0.01
wt = (1*10^9)*np.pi
graph_scaling_factor = 5

t = np.arange(0, N, step)
b1 = random.choices([0,1], weights = (50,50), k = int(L)) 



c = np.sin(t*wt)
cnt = 0
b1_mapped = []
for q in range(0,int(L)):
    for i in range(0,int((N/L)/step)):
        b1_mapped.append(b1[q])


s = []
for i in range(0,len(t)):
    s.append(b1_mapped[i]*c[i])
    


fig1, (sub1, sub2) = plt.subplots(1, 2, figsize=(10, 10))
sub1.plot([t[i] for i in range(1,int(N/graph_scaling_factor))],[b1_mapped[i] for i in range(1,int(N/graph_scaling_factor))])
sub1.plot([t[i] for i in range(1,int(N/graph_scaling_factor))],[s[i] for i in range(1,int(N/graph_scaling_factor))], linestyle = 'dashed')
sub2.plot(t, s)
plt.show()

