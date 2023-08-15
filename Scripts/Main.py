def Store_prop_data(line_1,line_2,home_lat,home_long,height,file_name):
    from Doplar_shift import Get_orbital, Get_freq_change, Get_range_rate
    from TLE import TLE_calc
    import csv
    from datetime import datetime
    print("Calculating trajectory")
    sat_data = Get_orbital(line_1, line_2, home_lat, home_long, height, 0)
    print("Done")
    
    d = sat_data[3]
    range_rate = Get_range_rate(d)
    freq_change= Get_freq_change(range_rate)
    data = []
    for i in range(len(sat_data[4])):
        time = sat_data[4][i].strftime("%m/%d/%Y %H:%M:%S:%f")
        elv = sat_data[2][i]
        height_of_sat = sat_data[6][i]
        d = sat_data[3][i]
        data.append([time,elv,height_of_sat,d,freq_change[i]])
    
        
    
        
        
    with open(file_name, 'w', encoding='UTF8', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Time","elv","Height","Slant Distance","Freq Change"])
        writer.writerows(data)
   
def Store_atten_data(prop_file_path,atmos_file_path,f_orig,wvd,ground_height):
    from Atmosphere_attenuation import Total_atmos_atten_space_earth
    import csv
    from datetime import datetime, timedelta
    with open(prop_file_path, "r",newline = "") as file:
        csvreader  = csv.reader(file)
        data = []
        for row in csvreader:
            data.append(row)
            
    A_g = [];time = []
    data = data[1:]
    for i in range(len(data)):
        # print(data[i][0])
        time = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
        if i == 0:
            A_g1 = Total_atmos_atten_space_earth(ground_height, float(data[i][2]), -float(data[i][1]), f_orig*float(data[i][4]), wvd)
            time_since_calc = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
        elif time_since_calc + timedelta(0,10)<time:
            A_g1 = Total_atmos_atten_space_earth(ground_height, float(data[i][2]), -float(data[i][1]), f_orig*float(data[i][4]), wvd)
            time_since_calc = (datetime.strptime(data[i][0], "%m/%d/%Y %H:%M:%S:%f"))
        
        A_g.append([time,A_g1[0],A_g1[1]])
            
        print(f"Calculating Atmosphere Data: {round(((i+1)/(len(data)))*100,2)}% at time {data[i][0]}")
            
    
    
    with open(atmos_file_path, 'w', encoding='UTF8', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Time","A_g","ducting"])
        writer.writerows(A_g)

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
                     float(propdata[i][2]),float(propdata[i][3]),float(propdata[i][4])])
           
    return data
 
    
def Get_data_current(sat_data,time_now):
    import pytz
    i = 0
    if time_now<sat_data[i][0].replace(tzinfo=pytz.utc):
        raise IndexError("Generate new data as sattellite out of date")
    try:
        while sat_data[i][0].replace(tzinfo=pytz.utc)<time_now:
            i +=1
        return sat_data[i]
    except IndexError:
        raise IndexError("Generate new data as sattellite out of date")

def Gen_signal(N,ion_strength,f,sat_data_now):
    from Mary_BSK import Bit_generator,MPSK_Generator
    from Ionshpere_scintillation import Scintilate
    import numpy as np
    from Atmosphere_attenuation import Free_space_loss
    f_scaled = f/1E9
    A = 4
    step = 1/(N*10*(2*np.pi*f_scaled))
    bits = Bit_generator(N)
    t_scaled = np.arange(0,2*N/f_scaled,step)
    t = t_scaled/1E9
    i1_mapped,s_clean,c_clean = MPSK_Generator(bits, t_scaled, 1, f_scaled, A)
    
    
    
    time_stamp, A_g, elv,H_s,d,f_change = sat_data_now
    f_doppler = f_scaled*f_change
    i1_mapped_1,s_dirty,c_dirty = MPSK_Generator(bits, t_scaled, 1, f_doppler, A)
    
    A_f = Free_space_loss(f_doppler, d)
    A_tot = A_f+A_g
    for i in range(len(s_dirty)):
        s_dirty[i] = s_dirty[i]+Scintilate(ion_strength, t_scaled[i], s_dirty[i])
        
    
    return t,c_clean,s_dirty,s_clean, A_tot