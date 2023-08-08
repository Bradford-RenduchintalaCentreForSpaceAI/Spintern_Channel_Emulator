def Get_orbital(line_1,line_2,lat,long,height,start_elv, acc = 100,stop_time_min = 60):
    """
    Predicts the next orbital pass for a sattellite given a location and height 

    Parameters
    ----------
    line_1 : string
        First line of TLE data.
    line_2 : string
        Second line of TLE data.
    lat : float
        lattitude of ground station in degrees.
    long : float
        Longitude of ground station in degrees.
    acc : int
        accuracy of calcs.
    height : int or float
        height of ground station in m.
    start_elv : float
        elevation angle to start observation in degrees.

    Returns
    -------
    list
        Datetime object of next closet window and how long until.
    list
        Duration of window and end date of window .
    elv : list
        Elevation data of the pass.
    d : list
        Slant range data of the pass.
    time_window : list
        time data of pass.

    """
    from TLE import TLE_calc
    import numpy as np
    import datetime
    from datetime import timezone
    
    sat = TLE_calc(line_1, line_2,False)
    d,Az1, elv1,q = sat.get_pos_TOPO(acc,lat,long,height)
    height_of_sat = []
    sec = 0
    step = 0.1
    stop_time_min = stop_time_min*60
    

    time_now = datetime.datetime.now(timezone.utc)
    while (np.rad2deg(elv1) < start_elv):
        sat = TLE_calc(line_1, line_2,True,time_now.year,time_now.month,time_now.day,time_now.hour,time_now.minute,
                      time_now.second,time_now.microsecond)
        d,Az1, elv1,q = sat.get_pos_TOPO(acc,lat,long,height)
        sec = sec+step 
        time_now = datetime.datetime.now(timezone.utc)+datetime.timedelta(0,sec)
    
        
    time_till_window = datetime.timedelta(seconds=sec)

    date_at_start = time_now

    current_elv = elv1
    time_window = []
    elv = [];Az = []
    d = []
    sec = 0
    while np.rad2deg(current_elv)>start_elv  :
        time_now = date_at_start+datetime.timedelta(0,sec)
        sat = TLE_calc(line_1, line_2,True,time_now.year,time_now.month,time_now.day,time_now.hour,time_now.minute,
                      time_now.second,time_now.microsecond)
        d_current,Az1, current_elv,q = sat.get_pos_TOPO(acc,lat,long,height)
        elv.append(current_elv);d.append(d_current);time_window.append(date_at_start+datetime.timedelta(seconds=sec));Az.append(Az1)
        height_of_sat.append(sat.get_height(acc))
        sec = sec+step
        if sec>stop_time_min:
            break
        
        
    
    
    time_till_end_of_window = datetime.timedelta(seconds=sec)
    date_at_end = date_at_start+datetime.timedelta(seconds=sec)
    return [time_till_window,date_at_start],[time_till_end_of_window,date_at_end], elv,d,time_window, Az,height_of_sat

def Oribtal_time_check():
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    line_1 = "1 37846U 11060A   23216.79436601 -.00000089  00000+0  00000+0 0  9992"
    line_2 = '2 37846  57.1053  10.1244 0002842  67.3316 292.6967  1.70475823 73395'
    lat = 51.00653
    long = 0.85587541
    acc = 100
    height = 0
    [time_till_window,date_at_start],[time_till_end_of_window,date_at_end], elv,d,time, Az= Get_orbital(line_1,line_2,
                                                                                                lat,long,acc,height,0)
    print(f"""
          Time till: {time_till_window} Elv start: {np.rad2deg(elv[0])}
          Duration: {time_till_end_of_window} Elv end: {np.rad2deg(elv[len(elv)-1])}
          Start time: {date_at_start} end time: {date_at_end}
          Start Slant: {np.rad2deg(d[0])} End Slant: {np.rad2deg(d[len(d)-1])}
          Max Elv: {max(np.rad2deg(elv))}
          """)

    fig1, (sub1, sub2) = plt.subplots(2,1, figsize=(10, 10))
    sub1.plot(time,np.rad2deg(elv))
    sub2.plot(time,d)
    sub1.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    sub2.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    sub1.grid(True)
    sub2.grid(True)
    sub1.set_ylabel("Degrees")
    sub2.set_ylabel("Range m")
    sub2.set_xlabel("Time")
    sub1.set_title("Elevation angle")
    sub2.set_title("Slant Range")
    plt.show()
    
    
    import csv
    
    with open("Test_data/slant_data.csv",mode = 'w',newline = '') as f:
        csv_writer = csv.writer(f)
        
        csv_writer.writerow(["Time","d","elv"])
        
        for i in range(len(time)):
            csv_writer.writerow([time[i],d[i],elv[i]])
        
    f.close()

def Get_range_rate(d):
    """
    Gets the range rate change

    Parameters
    ----------
    d : list
        Range rate in Km*-.

    Returns
    -------
    range_rate : list
        range rate change in Km.

    """
    import numpy as np
    range_rate = np.diff(d,1)
    
    for i in range(len(range_rate)):
        if i == 0:
            ii = 0
            while range_rate[i] == 0:
                range_rate[i] = range_rate[ii]
                ii += 1
        else:
            if range_rate[i] == 0:
                range_rate[i] = range_rate[i-1]
    

    
    range_rate = range_rate.tolist()
    
    range_rate.append(range_rate[len(range_rate)-1])
    
    return range_rate
    
def Get_freq_change(range_rate):
    freq_change = []
    for i in range(len(range_rate)):
        freq_change.append((3E8)/(3E8+(range_rate[i]*1000)))
        
    return freq_change
        
def Doppler_shift_test():
    import datetime
    import csv
    import matplotlib.pyplot as plt

    
    time = [];d_old = [];elv = []
    
    with open("Test_data/slant_data.csv","r") as f:
        csvreader= csv.reader(f,delimiter=',')
        
        for row in csvreader:
            time.append(row[0]);d_old.append(row[1]);elv.append(row[2])
    d = []
    for i in range(1,len(d_old)):
        d.append(float(d_old[i]))
        elv[i] = float(elv[i])
        time[i] = datetime.datetime.strptime(time[i],"%Y-%m-%d %H:%M:%S.%f%z")

    range_rate = Get_range_rate(d)


    plt.plot((3E8)/(3E8+(range_rate*1000)))
    
def All_togther_now_test():
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import numpy as np
    # Data for the ISS and lat/long and height for Ham Street kent 
    line_1 = "1 25544U 98067A   23203.12476910  .00011930  00000-0  21722-3 0  9991"
    line_2 = '2 25544  51.6411 155.4137 0000523  72.9869 282.0755 15.49927938407206'
    lat = 51.00653
    long = 0.85587541
    height = 0
    start_elv = 0
    Pass_info = Get_orbital(line_1, line_2, lat, long, height, start_elv)
    Start_date = Pass_info[0]
    elv = Pass_info[2]
    d = Pass_info[3]
    time_window = Pass_info[4]
    Az = Pass_info[5]

    range_rate = Get_range_rate(d)

    f_change = Get_freq_change(range_rate)

    
    
    
    fig1, (sub1, sub2,sub3) = plt.subplots(3,1, figsize=(10, 10))
    sub1.plot(time_window,np.rad2deg(elv))
    sub2.plot(time_window,f_change)
    sub3.plot(time_window,np.rad2deg(Az))

    
    
    sub1.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    sub2.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    sub3.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    sub1.grid(True,"minor");sub2.grid(True,"minor");sub3.grid(True,"minor")
    
    
    
    
    
    sub1.set_ylabel("Degrees")
    sub2.set_ylabel("Change in freq")
    sub3.set_ylabel("Degrees")
    sub3.set_xlabel("Time")
    sub1.set_title("Elevation angle")
    sub2.set_title("Frequency change")
    sub3.set_title("Azmith Angle (wrong)")
    format1 = "%Y/%m/%d %H:%M:%S"
    fig1.suptitle((f"Data for next ISS pass in {Start_date[0]} at {Start_date[1].strftime(format1)} (UTC)"),
                  verticalalignment = 'center', fontsize = 'xx-large', y= 0.93)
    plt.show()
    
if __name__ == "__main__":
    Oribtal_time_check()
