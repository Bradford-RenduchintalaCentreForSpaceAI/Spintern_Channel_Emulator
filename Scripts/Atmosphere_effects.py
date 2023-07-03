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
    
    
    
    
    
    