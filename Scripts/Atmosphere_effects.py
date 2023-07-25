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
    
    
def Gamma_w_o(f,pressure,p_w_v_den,t):
    # From ITU P.676-5 Annex 2
    
    if f >350:
        raise ValueError("Frequency too high")
    
    """ Gamma_o Calc"""
    r_p = pressure/1013
    r_t = 288/(273+t)
    
    gamma_o_dash_54 = 2.128*(r_p**1.4954)*(r_t**-1.6032)*np.exp(-2.5280*(1-r_t)) # 22e
    gamma_o_54 = 2.136*(r_p**1.4975)*(r_t**-1.5852)*np.exp(-2.5196*(1-r_t)) # 22f
    gamma_o_57 = 9.984*(r_p**0.9313)*(r_t**2.6732)*np.exp(0.8563*(1-r_t)) # 22g
    gamma_o_60 = 15.42*(r_p**0.8595)*(r_t**3.6178)*np.exp(1.1521*(1-r_t)) # 22h
    gamma_o_63 = 10.63*(r_p**0.9298)*(r_t**2.3284)*np.exp(0.6287*(1-r_t)) # 22i
    gamma_o_66 = 1.944*(r_p**1.6673)*(r_t**-3.3583)*np.exp(-4.1612*(1-r_t)) # 22j
    gamma_o_66_dash = 1.935*(r_p**1.6657)*(r_t**-3.3714)*np.exp(-4.1643*(1-r_t)) #22k
    
    eta_1 = 6.7665*(r_p**-0.5050)*(r_t**0.5106)*np.exp(1.5663*(1-r_t))-1 #22n
    eta_2 = 27.8843*(r_p**-0.4908 )*(r_t**0.8491)*np.exp(0.5496*(1-r_t))-1 #22o
    
    a = np.log(eta_2/eta_1)/np.log(3.5)#22l
    b = (4**a)/eta_1 #22m
    
    xi_1 = 6.9575*(r_p**-0.3461 )*(r_t**0.2535)*np.exp(1.3766*(1-r_t))-1 # 22r
    xi_2 = 42.1309*(r_p**-0.3068)*(r_t**1.2023)*np.exp(2.5147*(1-r_t))-1 # 22s
    
    c = np.log(xi_2/xi_1)/np.log(3.5) #22p
    d = (4**c)/xi_1 # 22q
    
    if f <= 60:
        N = 0
    else:
        N = -15
        
    if f <= 54:
        # 22a
        term_1 = (7.34*(r_p**2)*(r_t**3))/((f**2)+0.36*(r_p**2)*(r_t**2))
        term_2 = (0.3429*b*gamma_o_dash_54)/(((54-f)**a)+b)
        gamma_o = (term_1+term_2)*(f**2)*1E-3 #22a
        
    elif f > 54 and f < 66:
        term_1 = ((54**-N)*np.log(gamma_o_54)*(f-57)*(f-60)*(f-63)*(f-66))/1944
        term_2 = ((-57**-N)*np.log(gamma_o_57)*(f-54)*(f-60)*(f-63)*(f-66))/486
        term_3 = ((60**-N)*np.log(gamma_o_60)*(f-54)*(f-57)*(f-63)*(f-66))/324
        term_4 = ((-63**-N)*np.log(gamma_o_63)*(f-54)*(f-57)*(f-60)*(f-66))/486
        term_5 = ((66**-N)*np.log(gamma_o_66)*(f-54)*(f-57)*(f-60)*(f-63))/1944
        gamma_o = np.exp((term_1+term_2+term_3+term_4+term_5)*(f**N)) #22b
        
    elif f >= 66 and f < 120:
        term_1 = (0.2296*d*gamma_o_66_dash)/(((f-66)**c)+d)
        term_2 = (0.286*(r_p**2)*(r_t**3.8))/(((f-118.75)**2)+2.97*(r_p**2)*(r_t**1.6))
        gamma_o = ((term_1+term_2)*(f**2))*1E-3##22c
    elif f>= 120 and f <= 350:
        term_1 = 3.02*1E-4*(r_p**2)*(r_t**3.5)
        term_2 = (1.5827*(r_p**2)*(r_t**3))/((f-66)**2)
        term_3 = (0.286*(r_p**2)*(r_t**3.8))/(((f-118.75)**2)+2.97*(r_p**2)*(r_t**1.6))
        gamma_o = (term_1+term_2+term_3)*(f**2)*1E-3#22d
    
        
    """Gamma_w Calc"""
    
    xi_w_1 = 0.9544*r_p*(r_t**0.69)+0.0061*p_w_v_den #23b
    xi_w_2 = 0.95*r_p*(r_t**0.64)+0.0067*p_w_v_den #23c
    xi_w_3 = 0.9561*r_p*(r_t**0.67)+0.0059*p_w_v_den #23d
    xi_w_4 = 0.9543*r_p*(r_t**0.68)+0.0061*p_w_v_den #23e
    xi_w_5 = 0.955*r_p*(r_t**0.68)+0.006*p_w_v_den #23f
    
    g_22 = 1+((f-22.235)**2)/((f+22.235)**2) # 23g
    g_557 = 1+((f-557)**2)/((f+557)**2)
    g_752 = 1+((f-752)**2)/((f+752)**2)
    
    
    #23a
    term_1 = (3.13E-2)*r_p*(r_t**2)
    term_2 = (1.76E-3)*p_w_v_den*(r_t**8.5)
    term_3 = (3.84*xi_w_1*g_22*np.exp(2.23*(1-r_t)))/(((f-22.235)**2)+9.42*(xi_w_1**2))
    term_4 = (10.48*xi_w_2*np.exp(0.7*(1-r_t)))/(((f-183.31)**2)+9.48*(xi_w_2**2))
    term_5 = (0.078*xi_w_3*np.exp(6.4385*(1-r_t)))/(((f-321.226)**2)+6.29*(xi_w_3**2))
    term_6 = (3.76*xi_w_4*np.exp(1.6*(1-r_t)))/(((f-325.153)**2)+9.22*(xi_w_4**2))
    term_7 = (26.36*xi_w_5*np.exp(1.09*(1-r_t)))/(((f-380)**2))
    term_8 = (17.87*xi_w_5*np.exp(1.46*(1-r_t)))/(((f-448)**2))
    term_9 = (883.7*xi_w_5*g_557*np.exp(0.17*(1-r_t)))/(((f-557)**2))
    term_10 = (302.6*xi_w_5*g_752*np.exp(0.41*(1-r_t)))/(((f-752)**2))
    
    
    gamma_w = (term_1+term_2+(r_t**2.5)*(term_3+term_4+term_5+term_6+term_7+term_8+term_9+term_10)*(f**2))*p_w_v_den*(1E-4)    
    
    
    return gamma_o,gamma_w


def Gamma_w_o_test():
    import matplotlib.pyplot as plt
    
    f = np.linspace(1,350,10000)
    pressure  = 1013.25
    p_w_v_den = 7.5
    t = 15
    
    gamma_o = []
    gamma_w = []
    gamma_sum = []
    for i in range(len(f)):
        gamma_o.append(Gamma_w_o(f[i],pressure,p_w_v_den,t)[0])
        gamma_w.append(Gamma_w_o(f[i],pressure,p_w_v_den,t)[1])
        gamma_sum.append(gamma_w[i]+gamma_o[i])
    plt.figure(figsize=(10,14.1))
    plt.plot(f,gamma_w)
    plt.plot(f,gamma_o)
    plt.plot(f,(gamma_sum), linestyle = '--')
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True,"minor")
    plt.xlabel("Frequency f (GHz)")
    plt.ylabel("Specific Attenuation dB/Km")
    plt.legend(["H20", "Dry air","tot"])
    plt.title("Specific attenuation for 1013 hPa 15°C and a water vapour of 7.5 g/m³")
    plt.ylim((10E-4,10E2))
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

def water_height(f):
    if f<1 or f>350:
        ValueError("Frequency too low")
        
    if f>=1 and f<=56:
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
    
def water_vapour_pressure(T,h):
    #from ITU-R P.835-6
    import numpy as np
    ro = 7.5*np.exp(-h/2)
    
    e = (ro*T)/216.7
    
    return e
    
if __name__ == "__main__":
    Earth_atmosphere_model_test()
    
    
    