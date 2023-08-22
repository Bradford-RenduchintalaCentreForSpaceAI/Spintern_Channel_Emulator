def Store_prop_data(line_1,line_2,home_lat,home_long,height,file_name):
    from Doplar_shift import Get_orbital, Get_freq_change, Get_range_rate
    import csv
    from Ionshpere_scintillation import Get_TEC, NeQuick_interface
    print("Calculating trajectory")
    sat_data = Get_orbital(line_1, line_2, home_lat, home_long, height, 0)
    
    d = sat_data[3]
    range_rate = Get_range_rate(d)
    freq_change= Get_freq_change(range_rate)
    data = []
    for i in range(len(sat_data[4])):
        time = sat_data[4][i].strftime("%m/%d/%Y %H:%M:%S:%f")
        elv = sat_data[2][i]
        height_of_sat = sat_data[6][i]
        d = sat_data[3][i]
        lat = sat_data[7][i]
        long = sat_data[8][i]
        data.append([time,elv,height_of_sat,d,freq_change[i],lat,long])
    
        
    
        
    print("Writing to file")
    with open(file_name, 'w', encoding='UTF8', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Time","elv","Height","Slant Distance","Freq Change","lat","long"])
        writer.writerows(data)
    print("Done")
   
def Store_atten_data(prop_file_path,atmos_file_path,f_orig,wvd,home_lat, home_long,ground_height):
    from Atmosphere_attenuation import Total_atmos_atten_space_earth
    import csv
    from datetime import datetime, timedelta
    import numpy as np
    from Ionshpere_scintillation import NeQuick_interface, Get_TEC
    with open(prop_file_path, "r",newline = "") as file:
        csvreader  = csv.reader(file)
        data = []
        for row in csvreader:
            data.append(row)
            
    A_g = [];time = [];N_t = []; group_delay = [];Faraday_rotation= [];XPD = []
    data = data[1:]
    for i in range(len(data)):
        # print(data[i][0])
        time = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
        if i == 0:
            A_g1 = Total_atmos_atten_space_earth(ground_height, float(data[i][2]), float(data[i][1]), (f_orig*float(data[i][4])/1E9), wvd)
            time_since_calc = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
            Ne_Quick_data = NeQuick_interface(home_lat, home_long, ground_height, np.rad2deg(float(data[i][5])), np.rad2deg(float(data[i][6])), float(data[i][2]), time.strftime("%Y"),
                              time.strftime("%m"), time.strftime("%H"))
            N_t1, group_delay1,Faraday_rotation1,XPD1 = Get_TEC(Ne_Quick_data, f_orig*float(data[i][4]))
        elif time_since_calc + timedelta(0,10)<time:
            A_g1 = Total_atmos_atten_space_earth(ground_height, float(data[i][2]), float(data[i][1]), f_orig*float(data[i][4])/1E9, wvd)
            time_since_calc = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
            Ne_Quick_data = NeQuick_interface(home_lat, home_long, ground_height, np.rad2deg(float(data[i][5])), np.rad2deg(float(data[i][6])), float(data[i][2]), time.strftime("%Y"),
                              time.strftime("%m"), time.strftime("%H"))
            N_t1, group_delay1,Faraday_rotation1,XPD1 = Get_TEC(Ne_Quick_data, f_orig*float(data[i][4]))
           
            
        
        
        
        A_g.append([time,A_g1[0],A_g1[1],N_t1, group_delay1,Faraday_rotation1,XPD1])
        
        print(f"Calculating Atmosphere Data: {round(((i+1)/(len(data)))*100,2)}%")
            
    
    print("Writing to file")
    with open(atmos_file_path, 'w', encoding='UTF8', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Time","A_g","ducting","N_t", "group delay","Faraday rotation","XPD loss"])
        writer.writerows(A_g)
    print("Done")
    
def get_data(atmos_file_path,file_path_prop):
    import csv
    from datetime import datetime, timedelta
    with open(file_path_prop, "r",newline = "") as file:
        csvreader  = csv.reader(file)
        propdata = []
        for row in csvreader:
            propdata.append(row)
    
    with open(atmos_file_path, "r",newline = "") as file:
       csvreader  = csv.reader(file)
       atmosdata = []
       for row in csvreader:
           atmosdata.append(row)
           
    propdata = propdata[1:]
    atmosdata = atmosdata[1:]           
    data = []       
    for i in range(len(propdata)):
        data.append([datetime.strptime(propdata[i][0],"%m/%d/%Y %H:%M:%S:%f"),float(atmosdata[i][1]),float(propdata[i][1]),
                     float(propdata[i][2]),float(propdata[i][3]),float(propdata[i][4]),float(propdata[i][5]),float(propdata[i][6]),
                     str(atmosdata[i][2]),float(atmosdata[i][3]),float(atmosdata[i][4]),float(atmosdata[i][5]),float(atmosdata[i][6])])
           
    return data
     
def Get_data_current(sat_data,time_now):
    import pytz
    i = 0
    try:
        if time_now<sat_data[i][0].replace(tzinfo=pytz.utc):
            raise IndexError("Generate new data as sattellite out of date")
        try:
            while sat_data[i][0].replace(tzinfo=pytz.utc)<time_now:
                i +=1
            return sat_data[i]
    
        except IndexError:
            raise IndexError("Generate new data as sattellite out of date")
    except:
        if time_now<sat_data[i][0]:
            raise IndexError("Generate new data as sattellite out of date")
        try:
            while sat_data[i][0]<time_now:
                i +=1
            return sat_data[i]
    
        except IndexError:
            raise IndexError("Generate new data as sattellite out of date")
        

def AWGN(s,B,T,amp):
    import numpy as np
    k = 1.38E-23
    N = k*T*B
    SNR = amp/N
    N_o = amp/SNR
    N = []
    for i in range(len(s)):
        AWGN = np.random.normal(0,np.sqrt(N_o/2))
        N.append(AWGN)
        s[i] = s[i]+AWGN
    return s
        
def Gen_signal(N,ion_strength,f,sat_data_now,A,B,T):
    from Mary_BSK import Bit_generator,MPSK_Generator
    from Ionshpere_scintillation import Scintilate
    import numpy as np
    from Atmosphere_attenuation import Free_space_loss
    f_scaled = f/1E9
    step = 1/(N*10*(2*np.pi*f_scaled))
    bits = Bit_generator(N)
    t_scaled = np.arange(0,2*N/f_scaled,step)
    t = t_scaled/1E9
    i1_mapped,s_clean,c_clean = MPSK_Generator(bits, t_scaled, 1, f_scaled, A)
    c = 2997924580
    time_stamp = sat_data_now[0]
    A_g= sat_data_now[1]
    elv = sat_data_now[2]
    H_s = sat_data_now[3]
    d = sat_data_now[4]
    f_change = sat_data_now[5]
    lat = sat_data_now[6]
    long = sat_data_now[7]
    f_doppler = f_scaled*f_change
    i1_mapped_1,s_dirty,c_dirty = MPSK_Generator(bits, t_scaled, 1, f_doppler, A)
    
    
    
    
    s_dirty = AWGN(s_dirty, B, T, A)
    
    travel_time = (d*1E3)/c
    
    A_f = Free_space_loss(f_doppler, d)
    
    A_tot = A_f+A_g
    

    for i in range(len(s_dirty)):
        s_dirty[i] = s_dirty[i]+Scintilate(ion_strength, t_scaled[i], s_dirty[i])
        
    
    return t,c_clean,s_dirty,s_clean, A_tot,travel_time

    