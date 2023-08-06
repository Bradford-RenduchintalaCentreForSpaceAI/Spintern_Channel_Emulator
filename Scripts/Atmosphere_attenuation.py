import numpy as np
    
def Gamma_const(f,t,pressure,e):
    """
    

    Parameters
    ----------
    f : float
        Frequency in GHZ.
    t : float
        Temp in kelvin.
    pressure : f,loat
        Dry air pressure in hPA.
    e : float
        water vapour pressure hPa.

    Returns
    -------
    float
        Specific gaseous attenuation in dB/Km.

    """
    import numpy as np
    # From ITU P.676-13
    theta = 300/(t)
    p = pressure-e
    
    
    """Table 1"""
    f_o_o = [50.474214, 50.987745, 51.50336, 52.021429, 52.542418, 53.066934, 53.595775, 54.130025, 54.67118, 55.221384,
             55.783815, 56.264774, 56.363399, 56.968211, 57.612486, 58.323877, 58.446588, 59.164204, 59.590983,
             60.306056, 60.434778, 61.150562, 61.800158, 62.41122, 62.486253, 62.997984, 63.568526, 64.127775,
             64.67891, 65.224078, 65.764779, 66.302096, 66.836834, 67.369601, 67.900868, 68.431006, 68.960312,
             118.750334, 368.498246, 424.76302, 487.249273, 715.392902, 773.83949, 834.145546]
   
    a_1 =[0.975, 2.529, 6.193, 14.32, 31.24, 64.29, 124.6, 227.3, 389.7, 627.1, 945.3, 543.4, 1331.8, 1746.6,
          2120.1, 2363.7, 1442.1, 2379.9, 2090.7, 2103.4, 2438.0, 2479.5, 2275.9, 1915.4, 1503.0, 1490.2, 1078.0,
          728.7, 461.3, 274.0, 153.0, 80.4, 39.8, 18.56, 8.172, 3.397, 1.334, 940.3, 67.4, 637.7, 237.4, 98.1, 572.3,
          183.1]
   
    a_2 = [9.651, 8.653, 7.709, 6.819, 5.983, 5.201, 4.474, 3.8, 3.182, 2.618, 2.109, 0.014, 1.654, 1.255, 0.91, 0.621,
           0.083, 0.387, 0.207, 0.207, 0.386, 0.621, 0.91, 1.255, 0.083, 1.654, 2.108, 2.617, 3.181, 3.8, 4.473, 5.2,
           5.982, 6.818, 7.708, 8.652, 9.65, 0.01, 0.048, 0.044, 0.049, 0.145, 0.141, 0.145]
   
    a_3 = [6.69, 7.17, 7.64, 8.11, 8.58, 9.06, 9.55, 9.96, 10.37, 10.89, 11.34, 17.03, 11.89, 12.23, 12.62, 12.95,
           14.91, 13.53, 14.08, 14.15, 13.39, 12.92, 12.63, 12.17, 15.13, 11.74, 11.34, 10.88, 10.38, 9.96, 9.55,
           9.06, 8.58, 8.11, 7.64, 7.17, 6.69, 16.64, 16.4, 16.4, 16.0, 16.0, 16.2, 14.7]
    a_4 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           0.0, 0.0]
   
    a_5 = [2.566, 2.246, 1.947, 1.667, 1.388, 1.349, 2.227, 3.17, 3.558, 2.56, -1.172, 3.525, -2.378, -3.545,
           -5.416, -1.932, 6.768, -6.561, 6.957, -6.395, 6.342, 1.014, 5.014, 3.029, -4.499, 1.856, 0.658, -3.036,
           -3.968, -3.528, -2.548, -1.66, -1.68, -1.956, -2.216, -2.492, -2.773, -0.439, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    a_6 = [6.85, 6.8, 6.729, 6.64, 6.526, 6.206, 5.085, 3.75, 2.654, 2.952, 6.135, -0.978, 6.547, 6.451, 6.056, 0.436,
           -1.273, 2.309, -0.776, 0.699, -2.825, -0.584, -6.619, -6.759, 0.844, -6.675, -6.139, -2.895, -2.59, -3.68,
           -5.002, -6.091, -6.393, -6.475, -6.545, -6.6, -6.65, 0.079, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   
   
   
    """Table 2"""
    f_o_w = [22.23508, 67.80396, 119.99594, 183.310087, 321.22563, 325.152888, 336.227764, 380.197353, 390.134508,
             437.346667, 439.150807, 443.018343, 448.001085, 470.888999, 474.689092, 488.490108, 503.568532,
             504.482692, 547.67644, 552.02096, 556.935985, 620.700807, 645.766085, 658.00528, 752.033113, 841.051732,
             859.965698, 899.303175, 902.611085, 906.205957, 916.171582, 923.112692, 970.315022, 987.926764, 1780.0]
   
    b_1 = [0.1079, 0.0011, 0.0007, 2.273, 0.047, 1.514, 0.001, 11.67, 0.0045, 0.0632, 0.9098, 0.192, 10.41, 0.3254,
           1.26, 0.2529, 0.0372, 0.0124, 0.9785, 0.184, 497.0, 5.015, 0.0067, 0.2732, 243.4, 0.0134, 0.1325, 0.0547,
           0.0386, 0.1836, 8.4, 0.0079, 9.009, 134.6, 17506.0]
   
    b_2 = [2.144, 8.732, 8.353, 0.668, 6.179, 1.541, 9.825, 1.048, 7.347, 5.048, 3.595, 5.048, 1.405, 3.597, 2.379,
           2.852, 6.731, 6.731, 0.158, 0.158, 0.159, 2.391, 8.633, 7.816, 0.396, 8.177, 8.055, 7.914, 8.429, 5.11,
           1.441, 10.293, 1.919, 0.257, 0.952]
   
    b_3 = [26.38, 28.58, 29.48, 29.06, 24.04, 28.23, 26.93, 28.11, 21.52, 18.45, 20.07, 15.55, 25.64, 21.34, 23.2,
           25.86, 16.12, 16.12, 26.0, 26.0, 30.86, 24.38, 18.0, 32.1, 30.86, 15.9, 30.6, 29.85, 28.65, 24.08, 26.73,
           29.0, 25.5, 29.85, 196.3]
   
    b_4 = [0.76, 0.69, 0.7, 0.77, 0.67, 0.64, 0.69, 0.54, 0.63, 0.6, 0.63, 0.6, 0.66, 0.66, 0.65, 0.69, 0.61,
           0.61, 0.7, 0.7, 0.69, 0.71, 0.6, 0.69, 0.68, 0.33, 0.68, 0.68, 0.7, 0.7, 0.7, 0.7, 0.64, 0.68, 2.0]
   
    b_5 =[5.087, 4.93, 4.78, 5.022, 4.398, 4.893, 4.74, 5.063, 4.81, 4.23, 4.483, 5.083, 5.028, 4.506, 4.804,
          5.201, 3.98, 4.01, 4.5, 4.5, 4.552, 4.856, 4.0, 4.14, 4.352, 5.76, 4.09, 4.53, 5.1, 4.7, 5.15, 5.0,
          4.94, 4.55, 24.15]
   
    b_6 = [1.0, 0.82, 0.79, 0.85, 0.54, 0.74, 0.61, 0.89, 0.55, 0.48, 0.52, 0.5, 0.67, 0.65, 0.64, 0.72, 0.43, 0.45,
           1.0, 1.0, 1.0, 0.68, 0.5, 1.0, 0.84, 0.45, 0.84, 0.9, 0.95, 0.53, 0.78, 0.8, 0.67, 0.9, 5.0]

    d= (5.6E-4)*(p+e)*(theta**0.8) #Eq 9
   
   
   
    term_1 = d*(1+((f/d)**2))
    term_2 = (1.4E-12)*p*theta**1.5
    term_3 = 1+((1.9E-5)*f**1.5)
   
   
    N_d_dash = f*p*(theta**2)*(((6.14E-5)/term_1)+(term_2/term_3)) #Eq 9
   
    N_o = 0
    for i in range(len(f_o_o)):
        term_1 = (a_5[i]+(a_6[i]*theta))*1E-4
        term_2 = (p+e)*(theta**0.8)
        delta_o = (term_1)*term_2 #Eq 7
        term_1=term_2=0
        term_1 = p*theta**(0.8-a_4[i])
        term_2 = 1.1*e*theta
        delta_f_o = (a_3[i]*1E-4)*(term_1+term_2);term_1=term_2=0
        delta_f_o = np.sqrt((delta_f_o**2)+(2.25E-6)) #Eq 6a
        term_1 = (delta_f_o-(delta_o*(f_o_o[i]-f)))
        term_2 = ((f_o_o[i]-f)**2)+(delta_f_o**2)
        term_3 = (delta_f_o-(delta_o*(f_o_o[i]+f)))
        term_4 = ((f_o_o[i]+f)**2)+(delta_f_o**2)
        F_i_o = (f/f_o_o[i])*((term_1/term_2)+(term_3/term_4)) #Eq 5
        S_i_o = ((a_1[i]*10**(-7))*p*(theta**3)*np.exp(a_2[i]*(1-theta))) #Eq 3
        N_o += (S_i_o*F_i_o)
       
   
    N_o += N_d_dash

    N_w = 0
    for i in range(len(f_o_w)):
        term_1 = b_3[i]*1E-4
        term_2 = (p*(theta**(b_4[i])))
        term_3 = (b_5[i]*e*theta**(b_6[i]))
        delta_f_w = term_1*(term_2+term_3);term_1=term_2=term_3= 0 #Eq 5
        term_1 = 0.535*delta_f_w
        term_2 = 0.217*(delta_f_w**2)
        term_3 = ((2.1316E-12)*(f_o_w[i]**2))/(theta);term_1=term_2=term_3= 0
        term_1 = (delta_f_w)
        term_2 = ((f_o_w[i]-f)**2)+(delta_f_o**2)
        term_3 = (delta_f_w)
        term_4 = ((f_o_w[i]+f)**2)+(delta_f_w**2)
        F_i_w = (f/f_o_w[i])*((term_1/term_2)+(term_3/term_4)) #Eq 5
        term_1 = b_1[i]*1E-1
        term_2 = e*(theta**3.5)
        term_3 = np.exp(b_2[i]*(1-theta))
        S_i_w = term_1*term_2*term_3 #Eq 3
        N_w += S_i_w*F_i_w
   
    return (N_w+N_o)*0.1820*f

def Gamma_w_o_test():
    import matplotlib.pyplot as plt
    f = np.arange(0,1000,1)
    pressure  = 1013.25
    t = 288.15
    wvd = 7.5
    print("Start")
    gamma_o_wet = [];gamma_o_dry = []
    for i in range(len(f)):
        gamma_o_wet.append(Gamma_const(f[i],pressure,t,wvd))
        gamma_o_dry.append(Gamma_const(f[i],pressure,t,0))
    print("Done")
    plt.figure(figsize=(14,10))
    plt.plot(f,gamma_o_wet,label = "wet")
    plt.plot(f,gamma_o_dry,label = "dry")
    plt.yscale("log");plt.legend()
    plt.grid(True,which = 'minor',axis = 'both')
    plt.xlabel("Frequency f (GHz)")
    plt.ylabel("Specific Attenuation dB/Km")
    plt.title("Specific attenuation for 1013 hPa 15°C and a water vapour of 7.5 g/m³")
    plt.ylim((10E-4,10E5))
    plt.xlim((0,1000))
    plt.show()
  
def Zenith_attenuation_test():
    import numpy as np
    import matplotlib.pyplot as plt
    from Earth_atmospher import Earth_atmosphere_model, water_vapour_pressure
   
    f = np.arange(50,70,0.01)
    h = np.linspace(0,20,5)
    print("Start")
    Zenith = [[] for i in range(len(h))]
    for i in range(len(h)):
        cur_atmos = Earth_atmosphere_model(h[i])
        ro,e = water_vapour_pressure(cur_atmos[0], h[i],7.5)
        for ii in range(len(f)):
            Zenith[i].append(Gamma_const(f[ii], cur_atmos[0], cur_atmos[1], e))
            print(f"{round(((((ii/len(f))*i+1)/(len(h)))*100),2)}%")
       
       
    plt.figure(figsize=(20,20))
    for i in range(len(h)):
        plt.plot(f,Zenith[i],label = f'{h[i]}Km')
    plt.yscale("log")
    plt.legend()
    plt.xlim((50, 70))
    plt.ylim((10**(-3), 10**(2)))
    plt.grid(True,"minor")
    plt.show()

def Apparent_elv(P,T,e):
    """
    

    Parameters
    ----------
    P : float
        Dry air pressure.
    T : float
        Tempreture in Kelvin.
    e : float
        Water vapour pressure in hPA.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    #from TU-R P.453-14
    N = (77.6*(P/T))-(5.6*(e/T))+((3.75E5)*(e/T**2))
    return 1+N*1E-7

def Total_atmos_atten_earth_space(Height_of_ground,Height_of_sat,elv,f,wvd,space_to_earth):
    """
    

    Parameters
    ----------
    Height_of_ground : flat
        Height of ground station in Km.
    Height_of_sat : float
        Height of space station in Km.
    elv : float
        Elevation of satallite in rad.
    f : float
        Frequency in Ghz.
    wvd : float
        water vapour pressure at ground level hPA.


    Returns
    -------
    float
        Attenuation in dB.

    """
    import numpy as np
    from Earth_atmospher import Earth_atmosphere_model, water_vapour_pressure
   
    H_s = Height_of_sat
    H_g = Height_of_ground
    h = np.linspace(H_g,H_s,100)
    
    
    A_g_h = [];elv_array= []
    
    T_orig,P_tot_orig = Earth_atmosphere_model(H_g)
    e_orig = water_vapour_pressure(T_orig, H_g, wvd)[1]
    P_d_orig = P_tot_orig-e_orig
    
    n_orig= Apparent_elv(P_d_orig, T_orig, e_orig)
    orig_term = (6371+H_g)*n_orig
   
    for i in range(len(h)):
        T,P = Earth_atmosphere_model(h[i])
        ro,e = water_vapour_pressure(T, h[i], wvd)
        nq = Apparent_elv(P, T, e) 
        now_term = (6371+h[i])*nq
        cos_elv_apparent = (orig_term/now_term)*np.cos(elv)
        elv_array.append(np.arccos(cos_elv_apparent))
        term_1 = np.sqrt(1-(cos_elv_apparent**2))
        gamma = Gamma_const(f, T, P, e)
        A_g_h.append(gamma/term_1)
       
    A_g = np.trapz(A_g_h,x = h)
    return A_g
   
def Total_atmos_atten_space_earth(Height_of_ground,Height_of_sat,elv,f,wvd):
    """
    

    Parameters
    ----------
    Height_of_ground : flat
        Height of ground station in Km.
    Height_of_sat : float
        Height of space station in Km.
    elv : float
        Elevation of satallite in rad.
    f : float
        Frequency in Ghz.
    wvd : float
        water vapour pressure at ground level hPA.


    Returns
    -------
    A_g: float
        Attenuation in dB.
    ducting : bool
        Boolean saying if ducting is occouring or not.

    """
    import numpy as np
    from Earth_atmospher import Earth_atmosphere_model, water_vapour_pressure
   
    H_s = Height_of_sat
    H_g = Height_of_ground
    R_e = 6378.1
    ducting = False
    
    if H_s>=100:
        H_s1 = 100
    else:
        H_s1 = H_s
        
    
    T_s,P_s = Earth_atmosphere_model(H_s1)
    e_s = water_vapour_pressure(T_s, H_s1, wvd)[1]
    n_s = Apparent_elv(P_s, T_s, e_s)
    
    
    h = [];h_pointless = []
    for i in range(922):
        if i == 0:    
            h_pointless.append(0.0001*np.exp(i/100))
        else:
            h_pointless.append(h_pointless[i-1] + 0.0001*np.exp(i/100))
        
        
        if h_pointless[i]>=H_g:
            h.append(h_pointless[i])
    
    
    
    H_g = min(h)
    
    H_s = max(h[:len(h)-1])
    r_s = R_e+H_s
    
    
    A_g = 0
    for i in range(len(h)-1):
        T_n,P_n = Earth_atmosphere_model(h[i])
        e_n = water_vapour_pressure(T_n, h[i], wvd)[1]
        n_n = Apparent_elv(P_n, T_n, e_n)
        r_n_1 = R_e+h[i+1]
        r_n = R_e +h[i]
        
        #ducting check
        
        duct_term = np.cos(elv)*((R_e+H_s)*n_s)/((R_e+h[i])*n_n)
        
        if duct_term >= 1:
            A_g += 0
            ducting = True
        
        else:    
            cos_term = ((n_s/n_n)*r_s*np.cos(elv))**2
            
            term_1 = r_n_1**2
            term_2 = r_n**2
            
            term_3 = np.sqrt(term_1-cos_term)
            term_4 = np.sqrt(term_2-cos_term)
            
            lns = term_3 -term_4
            
            gamma = Gamma_const(f, T_n, P_n, e_n)
            
            
            A_g += gamma*lns
            
    return A_g, ducting
    
def Total_atmos_atten_test():
    import matplotlib.pyplot as plt
    Height_of_sat = 100
    Height_of_ground = 1
    f = 30
    step = 50
    elv_earth = np.linspace(0,np.pi/2,step)
    elv_space = np.linspace(-np.pi/2,np.deg2rad(-9.946),step)
    A_g_1 = {"12.5":[],"7.5":[],"2.5":[]}
    A_g_2 = {"12.5":[],"7.5":[],"2.5":[]}
    for i in range(len(elv_earth)):
        A_g_1["12.5"].append(Total_atmos_atten_earth_space(Height_of_ground, Height_of_sat, elv_earth[i], f, 12.5,False))
        A_g_1["7.5"].append(Total_atmos_atten_earth_space(Height_of_ground, Height_of_sat, elv_earth[i], f, 7.5,False))
        A_g_1["2.5"].append(Total_atmos_atten_earth_space(Height_of_ground, Height_of_sat, elv_earth[i], f, 2.5,False))
        
        A_g_2["12.5"].append(Total_atmos_atten_space_earth(Height_of_ground,Height_of_sat,elv_space[i],f,12.5)[0])
        A_g_2["7.5"].append(Total_atmos_atten_space_earth(Height_of_ground,Height_of_sat,elv_space[i],f,7.5)[0])
        A_g_2["2.5"].append(Total_atmos_atten_space_earth(Height_of_ground,Height_of_sat,elv_space[i],f,2.5)[0])
        print(f"{((i+1)/len(elv_earth))*100}%")
    
    fig1, (sub1, sub2) = plt.subplots(2,1, figsize=(10, 10))
    sub1.plot(np.rad2deg(elv_earth),A_g_1["2.5"], label = "2.5",linestyle = '-')
    sub1.plot(np.rad2deg(elv_earth),A_g_1["7.5"], label = "7.5",linestyle = "dashdot")
    sub1.plot(np.rad2deg(elv_earth),A_g_1["12.5"], label = "12.5",linestyle = '--')
    sub2.plot(np.rad2deg(elv_space),A_g_2["2.5"], label = "2.5",linestyle = '-')
    sub2.plot(np.rad2deg(elv_space),A_g_2["7.5"], label = "7.5",linestyle = "dashdot")
    sub2.plot(np.rad2deg(elv_space),A_g_2["12.5"], label = "12.5",linestyle = '--')
    sub1.set_title("Earth-Space path")
    sub2.set_title("space-Earth path")
    sub1.set_yscale("log")
    sub2.set_yscale("log")
    sub1.legend()
    sub2.legend()
    sub1.set_ylim(((0.1,100)))
    sub2.set_ylim(((0.1,100)))
    sub1.set_xlim(((-10,90)))
    sub2.set_xlim(((-90,0)))
    plt.show()

def elv_angle_test():
    from Earth_atmospher import Earth_atmosphere_model, water_vapour_pressure
    import numpy as np
    import matplotlib.pyplot as plt
    T_orig,P_tot_orig = Earth_atmosphere_model(0)
    e_orig = water_vapour_pressure(T_orig, 0, 7.5)[1]
    P_d_orig = P_tot_orig-e_orig
    
    n_orig= Apparent_elv(P_d_orig, T_orig, e_orig)
    orig_term = (6371)*n_orig
    h = np.linspace(1,100,1000)
    thi =[]
    for i in range(len(h)):
        T,P = Earth_atmosphere_model(h[i])
        ro,e = water_vapour_pressure(T, h[i], 7.5)
        nq = Apparent_elv(P, T, e)
        now_term = (6371+h[i])*nq
        thi.append(np.rad2deg(np.arccos((orig_term/now_term)*np.cos(np.deg2rad(-90)))))
    plt.figure()
    plt.plot(h,thi)
    print(np.trapz(thi,h)) 
          
if __name__ == "__main__":
    Total_atmos_atten_test()