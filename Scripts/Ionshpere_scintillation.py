import numpy as np
import matplotlib.pyplot as plt





def Scintilate(Type,t,s):
    from scipy.stats import nakagami
    # Theta estimations taken from 
    # F. D. Nunes and F. M. G. Sousa, "Generation of Nakagami correlated fading in GNSS signals affected by ionospheric scintillation,
    # " 2014 7th ESA Workshop on Satellite Navigation Technologies and European Workshop on GNSS Signals and Signal Processing (NAVITEC), 
    # Noordwijk, Netherlands, 2014, pp. 1-8, doi: 10.1109/NAVITEC.2014.7045141.
    if Type.lower() =="strong":
        S_4 = 0.8
        SD = 1.5
    elif Type.lower() == "moderate":
        S_4 = 0.5
        SD = 1/3
    elif Type.lower() == "weak":
        S_4 = 0.2
        SD = 2/15
    else:
        ValueError("Type is not a valid string try strong, weak or moderate")
    
    

    m = 1/(S_4**2)
    
    s_1 = [];noise_1 = []
    for i in range(len(t)):
        phase = (np.random.normal(0,SD))
        amp = (nakagami.rvs(m))
        print(phase)
        noise = np.exp(-1j*np.pi*(t[i]+phase))
        s_1.append(s[i]+noise)
        
        noise_1.append(np.exp(-1j*np.pi*(t[i]+np.pi/2)))
        
    return s_1



def Scinillilation_test():
    from Constalation_in_python import Constallation_phase
    A = 1
    f = 1
    t = np.linspace(0,10,10000)
    
    s = A*np.exp(-1j*np.pi*t*f)

    
    s_1 = Scintilate("weak", t, s)
    
    fig, (sub1, sub2) = plt.subplots(2,1, figsize=(10, 10))
    
    sub1.plot(t,s, label = "orig")
    sub1.plot(t,s_1, label = "noise")
    sub1.legend()
    
    
    real_part,phase_part = Constallation_phase(s, s)
    
    real_part_1,phase_part_1 = Constallation_phase(s, s_1)
    
    sub2.scatter(phase_part,real_part)
    sub2.scatter(phase_part_1,real_part_1)
    
    
    
    
    
    
    
    
    
    
if __name__ == "__main__":
    Scinillilation_test()