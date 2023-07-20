import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, rfft, rfftfreq
N = 100  
f =float(30**3)
step = 1/(2.5*f)
t = np.linspace(0,N,int(N/step))

s = 10*np.cos(f*t)

xf = rfft(s)[:rfft(s).size//2]
yf = rfftfreq(rfft(s).size,step)[:rfft(s).size//2]
fig2, ax2 = plt.subplots()
ax2.plot(xf,yf)


ax2.legend(["Bit train","Signal"])
plt.show()
 


fig2, ax2 = plt.subplots()
ax2.plot([t[i] for i in range(1,100)], [s[i] for i in range(1,100)])
ax2.legend(["Bit train","Signal"])
plt.show()
 

