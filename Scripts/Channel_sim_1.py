from TLE import TLE_calc
from Doplar_shift import Get_orbital, Get_freq_change, Get_range_rate
from Atmosphere_attenuation import Total_atmos_atten_space_earth, Free_space_loss
from Mary_BSK import *
import matplotlib.pyplot as plt
import numpy as np


"""Init Varibles"""

TLE_line_1 = "1 37846U 11060A   23216.79436601 -.00000089  00000+0  00000+0 0  9992"
TLE_line_2 = "2 37846  57.1053  10.1244 0002842  67.3316 292.6967  1.70475823 73395"


home_lat = 51.0643
home_long = 0.8598

sat_data = Get_orbital(TLE_line_1, TLE_line_2, home_lat, home_long, 0, 0,stop_time_min=1)

f_orig = 1E9

f_scaled = f_orig/1E9

wvd = 2.5

d = sat_data[3]

N = 5

end_point = 10

A = 4

step = (4*(f_scaled*2*3.14))

t = np.linspace(0,end_point,end_point*int(step))

c = A*np.exp(-1j*np.pi*2*t*f_scaled)

bits = Bit_generator(N)

ints = Int_generator(bits, 1)

s = MPSK_Generator(ints, t, 1, f_scaled, A)

"""Doppler stuff"""
f_change = Get_freq_change(Get_range_rate(d))
f = f_change
for i in range(len(f_change)):
    f[i] = f_change[i]*f_orig
    f[i] = f[i]/1E9
 
time = sat_data[4]




"""Atmosphere"""


elv = sat_data[2]


height_of_sat = sat_data[6]
Ducting = False

A_g = [];A_f = []

""" Free space Loss"""

for i in range(len(height_of_sat)):
    A_f.append(Free_space_loss(f[i],d[i]))
    


plt.plot(t,s[1])
plt.plot(t,s[0])





