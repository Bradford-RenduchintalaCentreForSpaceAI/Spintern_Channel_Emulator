import numpy as np

def Rain_attenuation(lat,elv_angle,h_s,R_zero,f,Pol_tilt,p):
    #Step 1
    
    h_r = 4-0.0075*(lat-36)
    
    #Step 2
    elv_angle = (elv_angle*np.pi)/180
    lat_of_earth_stat = (lat*np.pi)/180
    
    if elv_angle < (5*np.pi)/180:
        L_s = ((2*(h_r-h_s))/(np.sqrt(np.sin(elv_angle)**2+(2*(h_r-h_s)/8500*10**3))+np.sin(elv_angle)))
    else:
        L_s = (h_r-h_s)/np.sin(elv_angle)
        
    # Step 3
    L_g = L_s*np.sin(elv_angle)
    
    
    #Step 4
    
    # Coeffiecients for line fitting 
    k_h_coef = [[-5.3398,-0.35351,-0.23789,-0.94158],[-0.10008,1.26970, 0.86036,0.64552],[1.13098,0.45400,0.15354
                                                                                      ,0.16817],-0.18961,0.71147]
    k_v_coef = [[-3.80595,-3.44965,-0.39902,0.50167],[0.56934,-0.22911,0.73042,1.07319],[0.81061,0.51059,0.11899,0.27195],-0.16398,0.63297]
    
    a_h_coef = [[-0.14318,0.29591,0.32177,-5.37610,16.1721],[1.82442,0.77564,0.63773,-0.96230,-3.29980],[-0.55187
                                                                                           ,0.19822,0.13164,1.47828,3.43990],0.67849,-1.95537]
    a_v_coef = [[-0.07771,0.56727,-0.20238,-48.2991,48.5833],[2.33840,0.95545,1.14520,0.791669,0.791459],[-0.76284,0.54039,0.26809,0.116226,0.116479],-0.053739,0.83433]
    
    # K_h Caculation
    
    log_k_h = 0.0
    for j in range(4):
        term = k_h_coef[0][j] * np.exp(-((np.log10(f) - k_h_coef[1][j]) / k_h_coef[2][j]) ** 2)
        log_k_h += term
    log_k_h += k_h_coef[3] * np.log10(f) +k_h_coef[4]
    
    # K_v Caculation
    
    log_k_v = 0.0
    for j in range(4):
        term = k_v_coef[0][j] * np.exp(-((np.log10(f) - k_v_coef[1][j]) / k_v_coef[2][j]) ** 2)
        log_k_v += term
    log_k_v += k_v_coef[3] * np.log10(f) +k_v_coef[4]
    
    #a_h calc
    
    a_h = 0.0
    
    for j in range(5):
        term = a_h_coef[0][j] * np.exp(-((np.log10(f) - a_h_coef[1][j]) / a_h_coef[2][j]) ** 2)
        a_h += term
    a_h += a_h_coef[3] * np.log10(f) +a_h_coef[4]
    
    #a_h calc

    a_v = 0.0
    
    for j in range(5):
        term = a_v_coef[0][j] * np.exp(-((np.log10(f) - a_v_coef[1][j]) / a_v_coef[2][j]) ** 2)
        a_v += term
    a_v += a_v_coef[3] * np.log10(f) +a_v_coef[4]
    
    #Convert log values 
    k_h = 10**(log_k_h)
    k_v = 10**(log_k_v)
    
    k = (k_h+k_v+(k_h-k_v)*(np.cos(elv_angle)**2)*np.cos(2*Pol_tilt))/2
    a = (k_h*a_h+k_v*a_v+(k_h*a_h-k_v*a_v)*(np.cos(elv_angle)**2)*np.cos(2*Pol_tilt))/2
    
    #Step 5
    
    Gamma_r = k*(R_zero)**a
    
    
    #Step 6
    
    Hoz_reduction_factor = 1/(1+0.78*np.sqrt((L_g*Gamma_r)/f)-0.38*(1-np.exp(-2*L_g)))
    
    #Step 7
    
    zeta = np.arctan((h_r-h_s)/(L_g*Hoz_reduction_factor))
    
    if zeta>elv_angle:
        L_r = (L_g*Hoz_reduction_factor)/(np.cos(elv_angle))
    else:
        L_r = (h_r-h_s)/(np.sin(elv_angle))
    
    if np.abs(lat_of_earth_stat) < (36*np.pi)/180:
        pearson = ((36*np.pi)/180)-np.abs(lat_of_earth_stat)
    else:
        pearson = 0
    
    print((L_r*Gamma_r)/f**2)
    
    term = (31*(1-np.exp(-(elv_angle)/(1+pearson))*(np.sqrt(L_r*Gamma_r)/f**2)-0.45))
    vert_adjust_ang = 1/(1+np.sqrt(np.sin(elv_angle))*(term))
    
    
    #Step 8
    L_e = L_r*vert_adjust_ang
    
    #Step 9
    A_zero = Gamma_r*L_e
    
    #Step 10
    
    if p>=1 or np.abs(lat_of_earth_stat)>= (36*np.pi)/180:
        beta = 0
    elif p<1 and np.abs(lat_of_earth_stat) < (36*np.pi)/180 and np.abs(elv_angle)>= (25*np.pi)/180:
        beta = -0.005*(np.abs(lat_of_earth_stat)-36)
    else:
        beta = -0.005*(np.abs(lat_of_earth_stat)-36)+1.8-4.25*np.sin(elv_angle)
   
    term = -(0.655+0.033*np.log(p)-0.045*np.log(A_zero)-beta*(1-p)*np.sin(elv_angle))
    
    A_p = A_zero*(p/0.01)**term
    
    return A_p
        
def Gamma_const(f,pressure,p_w_v_den,t):
    import numpy as np
    # From ITU P.676-5 Annex 2
    theta = 300/(t)
    e = p_w_v_den
    p = pressure
    
    f_o_o = [50.474238, 50.987749, 51.50335, 52.02141, 52.542394, 53.066907, 53.595749, 54.13, 54.671159, 55.221367, 
             55.783802, 56.264775, 56.363389, 56.968206, 57.612484, 58.323877, 58.44659, 59.164207, 59.590983, 
             60.306061, 60.434776, 61.15056, 61.800154, 62.411215, 62.48626, 62.997977, 63.568518, 64.127767, 
             64.678903, 65.224071, 65.764772, 66.302091, 66.83683, 67.369598, 67.900867, 68.431005, 68.960311, 
             118.750343, 368.49835, 424.763124, 487.24937, 715.39315, 773.839675, 834.14533]
   
    a_1 =[0.94, 2.46, 6.08, 14.14, 31.02, 64.1, 124.7, 228.0, 391.8, 631.6, 953.5, 548.9, 1344.0, 1763.0, 2141.0,
            2386.0, 1457.0, 2404.0, 2112.0, 2124.0, 2461.0, 2504.0, 2298.0, 1933.0, 1517.0, 1503.0, 1087.0, 733.5, 
            463.5, 274.8, 153.0, 80.09, 39.46, 18.32, 8.01, 3.3, 1.28, 945.0, 67.9, 638.0, 235.0, 99.6, 671.0, 180.0]
    
    a_2 = [9.694, 8.694, 7.744, 6.844, 6.004, 5.224, 4.484, 3.814, 3.194, 2.624, 2.119, 0.015, 1.66, 1.26, 0.915,
             0.626, 0.084, 0.391, 0.212, 0.212, 0.391, 0.626, 0.915, 1.26, 0.083, 1.665, 2.115, 2.62, 3.195, 3.815, 
             4.485, 5.225, 6.005, 6.845, 7.745, 8.695, 9.695, 0.009, 0.049, 0.044, 0.049, 0.145, 0.13, 0.147]
    
    a_3 = [8.6, 8.7, 8.9, 9.2, 9.4, 9.7, 10.0, 10.2, 10.5, 10.79, 11.1, 16.46, 11.44, 11.81, 12.21, 12.66, 14.49, 
             13.19, 13.6, 13.82, 12.97, 12.48, 12.07, 11.71, 14.68, 11.39, 11.08, 10.78, 10.5, 10.2, 10.0, 9.7, 9.4,
             9.2, 8.9, 8.7, 8.6, 16.3, 19.2, 19.16, 19.2, 18.1, 18.1, 18.1]
    a_4 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.6, 
             0.6, 0.6]
    
    a_5 = [1.6, 1.4, 1.165, 0.883, 0.579, 0.252, -0.066, -0.314, -0.706, -1.151, -0.92, 2.881, -0.596, -0.556, -2.414
             ,-2.635, 6.848, -6.032, 8.266, -7.17, 5.664, 1.731, 1.738, -0.048, -4.29, 0.134, 0.541, 0.814, 0.415, 0.069
             ,-0.143, -0.428, -0.726, -1.002, -1.255, -1.5, -1.7, -0.247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    a_6 = [5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.514, 5.025, -0.069, 4.75, 4.104, 3.536,
             2.686, -0.647, 1.858, -1.413, 0.916, -2.323, -3.039, -3.797, -4.277, 0.238, -4.86, -5.079, -5.525,
             -5.52, -5.52, -5.52, -5.52, -5.52, -5.52, -5.52, -5.52, -5.52, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    
    f_o_w = [22.23508, 67.81396, 119.995941, 183.310074, 321.225644, 325.152919, 336.187, 380.197372, 390.134508,
             437.346667, 439.150812, 443.018295, 448.001075, 470.888947, 474.689127, 488.491133, 503.568532, 
             504.482692, 556.936002, 620.700807, 658.0065, 752.033227, 841.073593, 859.865, 899.407, 902.555, 
             906.205524, 916.171582, 970.315022, 987.926764]
    
    b_1 = [0.109, 0.0011, 0.0007, 2.3, 0.0464, 1.54, 0.001, 11.9, 0.0044, 0.0637, 0.921, 0.194, 10.6, 0.33, 1.28,
           0.253, 0.0374, 0.0125, 510.0, 5.09, 0.274, 250.0, 0.013, 0.133, 0.055, 0.038, 0.183, 8.56, 9.16, 138.0]
    b_2 = [2.143, 8.735, 8.356, 0.668, 6.181, 1.54, 9.829, 1.048, 7.35, 5.05, 3.596, 5.05, 1.405, 3.599, 2.381, 2.853, 
           6.733, 6.733, 0.159, 2.2, 7.82, 0.396, 8.18, 7.989, 7.917, 8.432, 5.111, 1.442, 1.92, 0.258]
    b_3 = [28.11, 28.58, 29.48, 28.13, 23.03, 27.83, 26.93, 28.73, 21.52, 18.45, 21.0, 18.6, 26.32, 21.52, 
           23.55, 26.02, 16.12, 16.12, 32.1, 24.38, 32.1, 30.6, 15.9, 30.6, 29.85, 28.65, 24.08, 26.7, 25.5, 
           29.85]
    b_4 = [0.69, 0.69, 0.7, 0.64, 0.67, 0.68, 0.69, 0.69, 0.63, 0.6, 0.63, 0.6, 0.66, 0.66, 0.65, 0.69, 0.61, 
           0.61, 0.69, 0.71, 0.69, 0.68, 0.33, 0.68, 0.68, 0.7, 0.7, 0.7, 0.64, 0.68]
    b_5 =[4.8, 4.93, 4.78, 5.3, 4.69, 4.85, 4.74, 5.38, 4.81, 4.23, 4.29, 4.23, 4.84, 4.57, 4.65, 5.04, 3.98,
          4.01, 4.11, 4.68, 4.14, 4.09, 5.76, 4.09, 4.53, 5.1, 4.7, 4.78, 4.94, 4.55]
    b_6 = [1.0, 0.82, 0.79, 0.85, 0.54, 0.74, 0.61, 0.84, 0.55, 0.48, 0.52, 0.5, 0.67, 0.65, 0.64, 0.72,
           0.43, 0.45, 1.0, 0.68, 1.0, 0.84, 0.45, 0.84, 0.9, 0.95, 0.53, 0.78, 0.67, 0.9]
    
    
    
    term_1 = ((3.57*(theta**7.5)*e)+0.113*p)
    
    N_dash_dash_w = f*(term_1)*(10E-7)*e*(theta**3)
    
    
    
    d = (5.6E-4)*(p+1.1*e)*theta
    
    term_1 = (6.14E-5/(d*(1+(f/d)**2)))
    term_2 = (1-(1.2E-5*(f**1.5)))
    
    N_dash_dash_d = f*p*(theta**2)*(term_1+1.4E-12*term_2*p*(theta**(1.5)))
    F_i_o = 0;s_i_o = 0
    for i in range(len(f_o_o)):
        delta_f = (a_3[i]*10E-4)*(p*(theta**(0.8-a_4[i]))+1.1*e*theta)
        delta = ((a_5[i]+(a_6[i]*theta))*10E-4)*p*(theta**(0.8))
        term_1 = ((delta_f-delta*(f_o_o[i]-f))/(((f_o_o[i]-f)**2)+delta_f**2))
        term_2 = ((delta_f-delta*(f_o_o[i]+f))/(((f_o_o[i]+f)**2)+delta_f**2))
        F_i_o += ((f/f_o_o[i])*(term_1+term_2))
        s_i_o += a_1[i]*10E-7*e*(theta**3)*np.exp(a_2[i]*(1-theta))
    F_i_w = 0;s_i_w = 0  
    for i in range(len(f_o_w)):
        delta_f = b_3[i]*10E-3*(p*(theta**b_4[i])+b_5[i]*e*(theta**b_6[i]))
        term_1 = ((delta_f)/(((f_o_o[i]-f)**2)+delta_f**2))
        term_2 = ((delta_f)/(((f_o_o[i]+f)**2)+delta_f**2))
        F_i_w += ((f/f_o_o[i])*(term_1+term_2))
        s_i_w +=(b_1[i]*10E-1*e*(theta**(3.5))*np.exp(b_2[i]*(1-theta)))
    
    

    Coef = (s_i_o+s_i_w)*(F_i_o+F_i_w)
    
    
    N_dash_dash = np.imag(Coef+N_dash_dash_d+N_dash_dash_w)
    
    #1082547444113268.4
    return 0.1820*f*N_dash_dash
    


def Gamma_w_o_test():
    import matplotlib.pyplot as plt
    
    f = np.linspace(1,350,10000)
    pressure  = 1013.25
    p_w_v_den = 7.5
    t = 15-273.7
    
    gamma_o = []
    for i in range(len(f)):
        gamma_o.append(Gamma_const(f[i],pressure,p_w_v_den,t))
        print(gamma_o[i])
    plt.figure(figsize=(10,14.1))
    plt.plot(f,gamma_o)
    # plt.xscale("log")
    plt.yscale("log")
    plt.grid(True,"minor")
    plt.xlabel("Frequency f (GHz)")
    plt.ylabel("Specific Attenuation dB/Km")
    plt.title("Specific attenuation for 1013 hPa 15°C and a water vapour of 7.5 g/m³")
    # plt.ylim((10E-4,10E2))
    plt.show()
        
def Earth_atmosphere_model(h):
    #from ITU P.835-6
    import numpy as np
    h = (6356.766*h)/(6356.766+h)
    a_o = 95.571899
    a_1 = -4.011801
    a_2 = 6.424731E-2
    a_3 = -4.789660E-4
    a_4 = 1.340543E-6
    
    
    if h < 11:
        T = 288.15-6.5*h
        p = 1013.25*((288.15)/(288.15-6.5*h))**(-34.1632/6.5)
        
    if h >11 and h<=20:
        T = 216.65
        p = 226.3226*np.exp(-34.1632*(h-11)/216.65)
        
    if h >20 and h<=32:
        T = 216.65+(h-20)
        p = 54.74980*((216.65)/(216.65+(h-20)))**(34.1632)
        
    if h>32 and h<=47:
        T = 228.65+2.8*(h-32)
        p = 8.680422*((216.65)/(216.65+(h-20)))**(34.1632/2.8)
        
    if h > 47 and h <= 51:
        T = 270.65-2.8*(h-51)
        p = 1.109106*np.exp(-34.1632*(h-47)/270.65)
        
    if h>51 and h <=71:
        T = 270.65 - 2.8*(h-51)
        p = 0.6694167*((270.65)/(270.65-2.8*(h-51)))**(-34.1632/2.8)
        
    if h>71 and h<=84.852:
        T = 214.65-2*(h-71)
        p = 0.03956649*((214.65)/((214.65-2.0*(h-71))))**(-34.1632/2.0)
    if h>84.852 and h<=86:
        T = 186.8673
        p = np.exp(a_o+(a_1*h)+(a_2*h**2)+(a_3*h**3)+(a_4*h**4))
        
    if h>86 and h<=91:
        T = 186.8673
        p = np.exp(a_o+(a_1*h)+(a_2*h**2)+(a_3*h**3)+(a_4*h**4))
        
    if h>91 and h<=100:
        T = 263.1905-76.3232*(1-((h-91)/(19.9429))**2)**0.5
        p = np.exp(a_o+(a_1*h)+(a_2*h**2)+(a_3*h**3)+(a_4*h**4))
        
    if h>100:
        T = p = 0
    
    return T,p
        
def Earth_atmosphere_model_test():
    import matplotlib.pyplot as plt
    import numpy as np
    
    h = np.linspace(10,99,100)
    
    T = []
    p = []

    for i in range(len(h)):
        T_h, p_h = Earth_atmosphere_model(i)
        T.append(T_h)
        p.append(p_h)
    
    
    fig1, (sub1, sub2) = plt.subplots(2,1, figsize=(10, 10))
    
    sub1.plot(T,h)
    sub2.plot(p,h)
    sub1.grid(True);sub2.grid(True,"minor")
    sub2.set_xscale("log")
    sub2.set_xlabel("pressure")
    sub1.set_xlabel("tempreture")
    sub1.set_ylabel("Height")
    sub2.set_ylabel("Height")

def Dry_water_height(f):
    if f<1 or f>350:
        ValueError("Frequency too low/hig=h")
        
    if f>=1 and f<=56.7:
        h_o = 5.386-(3.32734E-2*f)+(1.87185E-3*f**2)-(3.52087E-5*f**3)+((83.26)/((f-60)**2+1.2))
    if f>56.7 and f<=63.3:
        h_o = 10
    if f>63.3 and f<98.5:
        term_1 = 0.039581-(1.19751E-3*f)+(9.14810E-6*f**2)
        term_2 = 1-(0.028687*f)+(2.07858E-4*f**2)
        h_o = f*(term_1/term_2)+(90.6)/((f-60)**2)
    if f>=98.5 and f<=350:
        h_o = 5.542-(1.76414E-3*f)+(3.05354E-6*f**2)+(6.815)/((f-118.75)**2+0.321)
        
    
    return h_o
    
def Wet_water_height(f):
    if f<1 or f>350:
        ValueError("Frequency too low/hig=h")
    term_1 = 1.161/((f-22.23)**2+2.91)
    term_2 = 3.33/((f-183.3)**2+4.58)
    term_3 = 1.90/((f-325.1)**2+3.34)
    
    
    h_w = 1.65*(1+term_1+term_2+term_3)
    
    return h_w

def water_vapour_pressure(T,h):
    #from ITU-R P.835-6
    import numpy as np
    ro = 7.5*np.exp(-h/2)
    
    e = (ro*T)/216.7
    
    return e,ro
   
def Zenith_attenuation_test():
    import numpy as np
    import matplotlib.pyplot as plt
    f = np.arange(1,350,0.05)
    h = np.linspace(0,20,5)
    T = [];p = [];ro = [];Zenith = [[] for i in range(len(h))]
    for i in range(len(h)):
        T_current, p_current = Earth_atmosphere_model(h[i])
        T.append(T_current-273.7);p.append(p_current)
        e_current, ro_current = water_vapour_pressure(T_current, h[i]);ro.append(ro_current)
        for ii in range(len(f)):
            Gamma_const_current = Gamma_const(f[ii], p[i], ro[i], T[i])
            h_o = Dry_water_height(f[ii])
            h_w = Dry_water_height(f[ii])
            Zenith[i].append(Gamma_const_current[0]*h_o+Gamma_const_current[1]*h_w)
    
    
    plt.figure(figsize=(10,10))
    for i in range(len(h)):
        plt.plot(f,Zenith[i],label = f'{h[i]}Km')
    plt.yscale("log")
    plt.legend()
    plt.xlim((50, 70))
    plt.ylim((10**(-2), 10**(3)))
    plt.grid(True,"minor")
    plt.show()
    
            
    
    
    



if __name__ == "__main__":
    Gamma_w_o_test()
    
    
    