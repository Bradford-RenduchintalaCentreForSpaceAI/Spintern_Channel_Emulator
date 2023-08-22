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
        SD = 2.05
    elif Type.lower() == "moderate":
        S_4 = 0.5
        SD = 1/3
    elif Type.lower() == "weak":
        S_4 = 0.2
        SD = 2/15
    else:
        ValueError("Type is not a valid string try strong, weak or moderate")
    
    m = 1/(S_4**2)
    phase = (np.random.normal(0,SD)) 
    amp = (nakagami.rvs(m))
        
    noise = amp*np.exp(-1j*np.pi*(t+phase))


    return noise


    
    

def NeQuick_interface(home_lat,home_long,height_of_ground,sat_lat,sat_long,sat_height,year,month,hour):
    import os
    import subprocess
    import time
          
    if int(year)>2009:
        year = "2009"

    os.chdir("NeQuick2_P531-12")
    
    prog = subprocess.run(["neq2.exe"]
                          , capture_output= True,input=f"{home_lat},{home_long},{height_of_ground}\n{sat_lat},{sat_long},{sat_height}\n{year},{month},{hour}\nn\ny\n", encoding="utf8")

    print(prog.stderr)

    with open("slQu.dat","r") as file:
        data = file.read().split("\n")[15:]
    data1 = [[] for i in range(len(data))]
    for i in range(len(data)):
        data[i] = data[i].split(" ")
        for ii in range(len(data[i])):
            if data[i][ii] != "":
                
                try:
                    data1[i].append(float(data[i][ii]))
                except:
                    None
        
    
    
    data = data1[:len(data1)-3]
    

    os.chdir("..")
    
    return data

def Get_TEC(NeQuick_data,f):
    s = [];n_e = []
    for i in range(len(NeQuick_data)):
        s.append(NeQuick_data[i][0]*1000)
        n_e.append(NeQuick_data[i][5])
    import numpy as np
    
    B = 50E-6
    N_t = np.trapz(n_e,s)
    
    group_delay = 1.345*(N_t/((f)**2))*1E-7
    
    Faraday_rotation = (2.36E-14)*((B*N_t)/((f/1E9)**2))
    
    XPD = -20*np.log(np.tan(Faraday_rotation))
    
    return N_t, group_delay,Faraday_rotation,XPD
    
def Scinillilation_test_basic():
    from Constalation_in_python import Constallation_phase
    A = 1
    f = 1
    t = np.linspace(0,10,1000)
    
    s = A*np.exp(-1j*np.pi*t*f)

    s_1 = []
    
    for i in range(len(t)):
        s_1.append(s[i]+(Scintilate("weak", t[i], s[i])))
        

            
    
    
    fig, (sub1, sub2) = plt.subplots(2,1, figsize=(10, 10))
    
    sub1.plot(t,s, label = "orig")
    sub1.plot(t,s_1, label = "noise")
    sub1.legend()
    
    
    real_part,phase_part = Constallation_phase(s, s)
    
    real_part_1,phase_part_1 = Constallation_phase(s, s_1)
    
    sub2.scatter(real_part,phase_part)
    sub2.scatter(real_part_1,phase_part_1)
       
def Scinillilation_test_prop():
    from TLE import TLE_calc
    import numpy as np
    from datetime import datetime, timezone
    from Main import AWGN
    line_1 = "1 40128U 14050A   23230.09515268 -.00000077  00000+0  00000+0 0  9994"
    line_2 = "2 40128  49.9508 315.2575 1609708 139.7180 233.4163  1.85519304 61072"
    sat = TLE_calc(line_1, line_2)
    A = 4
    f = 1
    t = np.linspace(0,10,1000)
    
    s = A*np.exp(-1j*np.pi*t*f)
    f = 1202.025E6
    
    home_lat = 51.0643
    home_long = 0.8598
     
    d, Az, elv, pos = sat.get_pos_TOPO(100, home_lat, home_long, 0)
    
    q = sat.get_pos_ECEF(100)
    
    lat_of_sat = np.rad2deg(q["lat"])
    long_of_sat = np.rad2deg(q["long"])
    height = sat.get_height(100)
    
    
    time_now = datetime.now(tz=timezone.utc)
    
    year = time_now.strftime("%Y")
    
    hour = time_now.strftime("%H")
   
    month = time_now.strftime("%m")
    
    
    data = NeQuick_interface(home_lat,home_long,0,lat_of_sat,long_of_sat,height,year,month,hour)
    
    
    N_t, group_delay,Faraday_rotation,XPD = Get_TEC(data,f)
    
    
    s = AWGN(s, 100E6, 390, A)
    
    
    
    print(f"N_t:{N_t}\nGroup delay: {group_delay}\nFaraday rotation: {Faraday_rotation}\nXPD loss : {XPD}")
    
if __name__ == "__main__":
    Scinillilation_test_prop()