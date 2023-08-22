from Main import get_data,Store_atten_data,Store_prop_data,Gen_signal, Get_data_current, AWGN
from Atmosphere_attenuation import Free_space_loss
import matplotlib.pyplot as plt
from Constalation_in_python import Constallation_phase
from datetime import datetime
import numpy as np
import geopandas as gpd
import matplotlib.dates as mdates




"""INIT SETUP"""
file_path_prop = "Data\Prop_data.csv"
atmos_file_path = "Data\Atmos_data.csv"
line_1 = "1 37846U 11060A   23232.63214153 -.00000077  00000+0  00000+0 0  9996"
line_2 = "2 37846  57.1080   9.6970 0002304  72.9466 287.0814  1.70476050 73667"
home_lat = 51.4545
home_long = -2.5879
wvd = 2.5
ion_strength = 'weak'

# Store_prop_data(line_1,line_2,home_lat,home_long,0,file_path_prop)
# Store_atten_data(file_path_prop, atmos_file_path, 1202.025E6, wvd, home_lat, home_long, 0)
time_now = datetime.strptime("2023,08,22 11:30", "%Y,%m,%d %H:%M")
sat_data = get_data(atmos_file_path, file_path_prop)
sat_data_now = Get_data_current(sat_data, time_now)
"""SDR init vars"""
f = 1202.025E6
N = 2 
A = 4
B = 100E6
T = 400
#Signal Generation
t, c_clean,s_dirty,s_clean,A_tot,travel_time = Gen_signal(N, ion_strength, f, sat_data_now,A,B,T)
#Q map stuff
clean_real, clean_phase = Constallation_phase(c_clean, s_clean)
dirty_real, dirty_phase = Constallation_phase(c_clean,s_dirty)
# Sat data comprehension
time = [];elv= []; A_g = [];freq_change = [];lat= [];long = [];ducting = [];TEC = [];
ION_delay = [];Faraday_rot = [];XPD_loss = [];Tot_loss = [];d = []
for i in range(len(sat_data)):
    time.append(sat_data[i][0])
    elv.append(sat_data[i][2])
    A_g.append(sat_data[i][1])
    freq_change.append(sat_data[i][5])
    lat.append(sat_data[i][6])
    long.append(sat_data[i][7])
    ducting.append(sat_data[i][8])
    TEC.append(sat_data[i][9])
    ION_delay.append(sat_data[i][10])
    Faraday_rot.append(sat_data[i][11])
    XPD_loss.append(sat_data[i][12])
    d.append(sat_data[i][4])
    Tot_loss.append(A_g[i]+Free_space_loss((f*freq_change[i])/1E9, d[i])+XPD_loss[i])
#Getting better elv angles
newelv =[]
for i in range(len(elv)):
    if np.rad2deg(elv[i])>10:
        newelv.append(elv[i])     
time = time[:len(newelv)]
A_g = A_g[:len(newelv)]
freq_change = freq_change[:len(newelv)]
lat = lat[:len(newelv)]
long= long[:len(newelv)]
ducting = ducting[:len(newelv)]
TEC = TEC[:len(newelv)]
ION_delay = ION_delay[:len(newelv)]
Faraday_rot = Faraday_rot[:len(newelv)]
XPD_loss = Faraday_rot[:len(newelv)]
d= d[:len(newelv)]
Tot_loss = Tot_loss[:len(newelv)]

ducting_time = []
for i in range(len(ducting)):
    if ducting[i] == "True":
        ducting_time.append(time[i])

tot_loss_now = sat_data_now[1]+Free_space_loss((f*sat_data_now[5])/1E9, sat_data_now[4])+sat_data_now[12]

"""FIG1"""
fig1, (sub1_1,sub1_2) = plt.subplots(2,1, figsize=(10,7.5))
sub1_1.scatter(clean_real,clean_phase, label = "Without scintillation")
sub1_1.scatter(dirty_real,dirty_phase, label = "With scintillation")
sub1_1.set_xlabel("Real", size = 16)
sub1_1.set_ylabel("Imag", size = 16)

sub1_2.plot(t*1E9,c_clean, label = "Without scintillation", color = "tab:blue")
sub1_2.set_xlabel("Time (ns)",fontsize = 12)
sub1_2.set_ylabel("Amplitude",fontsize = 12)
sub1_2.tick_params(axis='x', labelcolor="tab:blue",size = 16)
sub1_2.ticklabel_format(useOffset=False)

sub1_3 = sub1_2.twiny()

sub1_3.plot((t+travel_time+sat_data_now[10])*1E9,s_dirty, label = "With scintillation and delay", color = "tab:orange")
sub1_3.set_ylabel("Amplitude",fontsize = 20)
sub1_3.set_xlabel("", color = "tab:orange",fontsize = 16)
sub1_3.tick_params(axis='x', labelcolor="tab:orange",size = 16)
sub1_3.ticklabel_format(useOffset=False)


fig1.suptitle("Signal Representation", fontsize = 20, y = 1)
sub1_1.set_title("Constallation plot of signal", size = 16)
sub1_2.set_title("Time series plot of signal",size = 16)
sub1_2.legend()
sub1_3.legend(loc = "upper left")
sub1_1.legend()
fig1.tight_layout()

"""FIG2"""
fig2, ax2 = plt.subplots(figsize=(10,5))
countries = gpd.read_file(
               gpd.datasets.get_path("naturalearth_lowres"))

countries.plot(color="green", ax = ax2)


ax2.plot(np.rad2deg(long),np.rad2deg(lat))
ax2.scatter(np.rad2deg(sat_data_now[7]),np.rad2deg(sat_data_now[6]),label = "Sattellite lat/long")
ax2.scatter(home_long,home_lat, label = "Bristol Lat/long")
ax2.set_xlabel("longitude")
ax2.set_ylabel("lattitude")
ax2.legend()
fig2.suptitle(f"Orbital window of GSAT0101 sattellite at {time_now} UTC",y= 0.95, fontsize = 16)
fig1.tight_layout()

"""FIG3"""
fig3, (sub3_1,sub3_2) = plt.subplots(1,2,figsize = (10,5))
sub3_1.plot(time,ION_delay)
sub3_2.plot(time,Tot_loss)
sub3_2.axvline(ducting_time[0], color = "tab:orange", label = "Ducting occours")
sub3_1.axvline(ducting_time[0], color = "tab:orange")
sub3_1.scatter(sat_data_now[0],sat_data_now[10])
sub3_2.scatter(sat_data_now[0],tot_loss_now)

sub3_2.set_xlabel("Time (Hour)")
sub3_1.set_ylabel("Ionspheric delay")
sub3_2.set_ylabel("Total loss (dB)")
sub3_1.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
sub3_2.xaxis.set_major_formatter(mdates.DateFormatter('%H'))

fig3.suptitle("Attenuation and group delay", fontsize = 16, y = 0.95)
fig3.tight_layout()
fig3.legend()