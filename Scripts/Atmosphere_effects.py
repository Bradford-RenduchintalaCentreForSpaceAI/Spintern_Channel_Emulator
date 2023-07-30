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
       
def Gamma_const(f,t,pressure,e):
    import numpy as np
    # From ITU P.676-5 Annex 2
    T = t
    theta = 300/(t)
   
   
    p = pressure-e
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

    d= (5.6E-4)*(p+e)*(theta**0.8)
   
   
   
    term_1 = d*(1+((f/d)**2))
    term_2 = (1.4E-12)*p*theta**1.5
    term_3 = 1+((1.9E-5)*f**1.5)
   
   
    N_d_dash = f*p*(theta**2)*(((6.14E-5)/term_1)+(term_2/term_3))
   
    N_o = 0
    for i in range(len(f_o_o)):
        term_1 = (a_5[i]+(a_6[i]*theta))*1E-4
        term_2 = (p+e)*(theta**0.8)
        delta_o = (term_1)*term_2
        term_1=term_2=0
        term_1 = p*theta**(0.8-a_4[i])
        term_2 = 1.1*e*theta
        delta_f_o = (a_3[i]*1E-4)*(term_1+term_2);term_1=term_2=0
        delta_f_o = np.sqrt((delta_f_o**2)+(2.25E-6))
        term_1 = (delta_f_o-(delta_o*(f_o_o[i]-f)))
        term_2 = ((f_o_o[i]-f)**2)+(delta_f_o**2)
        term_3 = (delta_f_o-(delta_o*(f_o_o[i]+f)))
        term_4 = ((f_o_o[i]+f)**2)+(delta_f_o**2)
        F_i_o = (f/f_o_o[i])*((term_1/term_2)+(term_3/term_4))
        S_i_o = ((a_1[i]*10**(-7))*p*(theta**3)*np.exp(a_2[i]*(1-theta)))
        N_o += (S_i_o*F_i_o)
       
   
    N_o += N_d_dash
    #N_o = np.imag(N_o)
    N_w = 0
    for i in range(len(f_o_w)):
        term_1 = b_3[i]*1E-4
        term_2 = (p*(theta**(b_4[i])))
        term_3 = (b_5[i]*e*theta**(b_6[i]))
        delta_f_w = term_1*(term_2+term_3);term_1=term_2=term_3= 0
        term_1 = 0.535*delta_f_w
        term_2 = 0.217*(delta_f_w**2)
        term_3 = ((2.1316E-12)*(f_o_w[i]**2))/(theta);term_1=term_2=term_3= 0
        term_1 = (delta_f_w)
        term_2 = ((f_o_w[i]-f)**2)+(delta_f_o**2)
        term_3 = (delta_f_w)
        term_4 = ((f_o_w[i]+f)**2)+(delta_f_w**2)
        F_i_w = (f/f_o_w[i])*((term_1/term_2)+(term_3/term_4))
        term_1 = b_1[i]*1E-1
        term_2 = e*(theta**3.5)
        term_3 = np.exp(b_2[i]*(1-theta))
        S_i_w = term_1*term_2*term_3
        N_w += S_i_w*F_i_w
    #N_w = np.imag(N_w)
       
       
       
       
   
   
    #1082547444113268.4
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
        raise(ValueError("Alltitude too high"))
   
    return T,p #T(k) #p(hPa)
       
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

def water_vapour_pressure(T,h,p_o):
    #from ITU-R P.835-6
    import numpy as np
    ro = p_o*np.exp(-h/2) #g/m^3
   
    e = (ro*T)/216.7 #hPa
   
    return e,ro
   
def Zenith_attenuation_test():
    import numpy as np
    import matplotlib.pyplot as plt
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
   
def Total_atmos_atten(Height_of_ground,Height_of_sat,elv,f,wvd):
    import numpy as np
    H_s = Height_of_sat
    H_g = Height_of_ground
    A_g_h = [];n = []
    h = np.linspace(H_g, H_s,1000)
   
    T_orig,P_tot_orig = Earth_atmosphere_model(H_g)
    e_orig = water_vapour_pressure(T_orig, H_g, wvd)[1]
    P_d_orig = P_tot_orig-e_orig
    N_orig = (77.6*(P_d_orig/T_orig))+(72*(e_orig/T_orig))+((3.75E5)*(e_orig/T_orig**2))
   
    n_orig = 1+N_orig*1E-6
   
    orig_term = (6371+H_g)*n_orig
   
    for i in range(len(h)):
        T,P = Earth_atmosphere_model(h[i])
        ro,e = water_vapour_pressure(T, h[i], wvd)
        P_d = P-e
        N = (77.6*(P_d/T))+(72*(e/T)+((3.75E5)*(e/T**2)))
        n =  1+N*1E-6
        now_term = (6371+h[i])*n
        cos_elv_apparent = (orig_term/now_term)*np.cos(elv)
        term_1 = np.sqrt(1-cos_elv_apparent**2)
        gamma = Gamma_const(f, T, P, e)
        A_g_h.append(gamma/term_1)
       
    A_g = np.trapz(A_g_h,x = h)
    return A_g
   
def Total_atmos_atten_test():
    import matplotlib.pyplot as plt
    import os
    Height_of_sat = 100
    Height_of_ground = 1.0
    f = 30
    elv = np.linspace(-np.pi/2,np.deg2rad(-10),10)
   
   
    A_g = []
    for i in range(len(elv)):
        A_g.append(Total_atmos_atten(Height_of_ground, Height_of_sat, elv[i], f,2.5))
        print(f"{(i/len(elv))*100}%")


    plt.figure(figsize=(14,10))
    plt.plot(np.rad2deg(elv),A_g)
    plt.yscale("log")
    plt.ylim((0.1,100))
    plt.xlim((-90,0))
    plt.grid(True,"minor")
    plt.show()
   
if __name__ == "__main__":
    Zenith_attenuation_test()
    # x = input("input ")
    # if x == "1":    
    #     Zenith_attenuation_test()
    # else:
    #     Total_atmos_atten_test()