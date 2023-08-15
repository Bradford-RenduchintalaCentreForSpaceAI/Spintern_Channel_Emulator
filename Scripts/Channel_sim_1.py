from Main import get_data,Store_atten_data,Store_prop_data,Gen_signal, Get_data_current
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
from Constalation_in_python import Constallation_phase
from datetime import datetime, timezone
import numpy as np
"""Init Sat vars"""
file_path_prop = "Data\Prop_data.csv"
atmos_file_path = "Data\Atmos_data.csv"
TLE_line_1 = "1 44850U 19088A   23225.08599465 -.00000023  00000+0  00000+0 0  9991"
TLE_line_2 = "2 44850  65.0694 108.5397 0011909 255.7022 104.2347  2.13102271 28559"
home_lat = 51.0643
home_long = 0.8598
wvd = 2.5
ion_strength = 'strong'

# Store_prop_data(TLE_line_1,TLE_line_2,home_lat,home_long,0,file_path_prop)
# Store_atten_data(file_path_prop,atmos_file_path,1,wvd,0)
time_now = datetime.now(timezone.utc)
sat_data = get_data(atmos_file_path, file_path_prop)
sat_data_now = Get_data_current(sat_data, time_now)
"""SDR init vars"""
f = 1202.025E6
N = 4

#Signal Generation
t, c_clean,s_dirty,s_clean,A_tot = Gen_signal(N, ion_strength, f, sat_data_now)

#Q map stuff
clean_real, clean_phase = Constallation_phase(c_clean, s_clean)
dirty_real, dirty_phase = Constallation_phase(c_clean,s_dirty)






#Plotting signal stuff
fig1, (sub1_1,sub1_2) = plt.subplots(2,1, figsize=(10, 10))

sub1_1.scatter(clean_real,clean_phase)
sub1_1.scatter(dirty_real,dirty_phase)

sub1_2.plot(t,c_clean)
sub1_2.plot(t,s_dirty)




# Sat data comprehension
time = [];elv= []; A_g = [];freq_change = []
for i in range(len(sat_data)):
    time.append(sat_data[i][0])
    elv.append(sat_data[i][2])
    A_g.append(sat_data[i][1])
    freq_change.append(sat_data[i][5])
    

#Getting better elv angles
newelv =[]
for i in range(len(elv)):
    if np.rad2deg(elv[i])>10:
        newelv.append(elv[i])
        
time = time[:len(newelv)]
A_g = A_g[:len(newelv)]
freq_change = freq_change[:len(newelv)]



fig2, (sub2_1,sub2_2) = plt.subplots(2,1, figsize=(10, 10))

sub2_1.plot(time,np.rad2deg(newelv))
sub2_1.scatter(sat_data_now[0],np.rad2deg(sat_data_now[2]))
sub2_2.plot(time,freq_change)
sub2_2.scatter(sat_data_now[0],sat_data_now[5])
sub2_2.xaxis.set_major_formatter(mdates.DateFormatter('%D:%H'))
sub2_1.xaxis.set_major_formatter(mdates.DateFormatter('%D:%H'))



fig3, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot(np.rad2deg(newelv),A_g,freq_change)
ax.view_init(20,-120,0)