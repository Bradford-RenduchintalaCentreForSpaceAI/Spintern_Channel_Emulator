def Get_orbital(line_1,line_2,lat,long,acc,height,start_elv):
    from TLE import TLE_calc
    import numpy as np
    import datetime
    from datetime import timezone
    
    init = datetime.datetime.now(timezone.utc)
    sat = TLE_calc(line_1, line_2,False)
    d,Az1, elv1,q = sat.get_pos_TOPO(acc,lat,long,height)
    sec = 0
    step = 0.1
    
    elv = []
    Az = []
    
    
    while (np.rad2deg(elv1) < start_elv):
        time_now = datetime.datetime.now(timezone.utc)+datetime.timedelta(0,sec)
        sat = TLE_calc(line_1, line_2,True,time_now.year,time_now.month,time_now.day,time_now.hour,time_now.minute,
                      time_now.second,time_now.microsecond)
        d,Az1, elv1,q = sat.get_pos_TOPO(acc,lat,long,height)
        sec = sec+step 
        
    elv.append(elv1);Az.append(Az1)
    time_till_window = datetime.timedelta(seconds=sec)
        
    date_at_start = init+datetime.timedelta(seconds=sec)
    
    current_elv = elv[0]
    time_window = [];time_window.append(date_at_start)
    sec = 0
    
    while np.rad2deg(current_elv)+0.5>1:
        time_now = date_at_start+datetime.timedelta(0,sec)
        sat = TLE_calc(line_1, line_2,True,time_now.year,time_now.month,time_now.day,time_now.hour,time_now.minute,
                      time_now.second,time_now.microsecond)
        d,Az1, current_elv,q = sat.get_pos_TOPO(acc,lat,long,height)
        elv.append(current_elv);Az.append(Az1);time_window.append(date_at_start+datetime.timedelta(seconds=sec))
        sec = sec+step
        
        
    
    
    time_till_end_of_window = datetime.timedelta(seconds=sec)
    date_at_end = date_at_start+datetime.timedelta(seconds=sec)
    return [time_till_window,date_at_start],[time_till_end_of_window,date_at_end], elv,Az,time_window




def test():
    import numpy as np
    import matplotlib.pyplot as plt
    line_1 = "1 25544U 98067A   23202.15851711  .00013480  00000-0  24456-3 0  9998"
    line_2 = '2 25544  51.6395 160.2062 0000289  47.0684 313.0329 15.49904890407063'
    lat = 51.00653
    long = 0.85587541
    acc = 100
    height = 0

    
    
    [time_till_window,date_at_start],[time_till_end_of_window,date_at_end], elv,Az,time = Get_orbital(line_1,line_2,
                                                                                                lat,long,acc,height,1)
    
    
    print(f"""
          Time till: {time_till_window} Elv start: {np.rad2deg(elv[0])}
          Duration: {time_till_end_of_window} Elv end: {np.rad2deg(elv[len(elv)-1])}
          Start time: {date_at_start} end time: {date_at_end}
          Start Az: {np.rad2deg(Az[0])} End Az: {np.rad2deg(Az[len(Az)-1])}
          
          
          """)
          
    fig = plt.figure()
    sub1 = plt.subplot(121)
    sub2 = plt.subplot(122)
    sub1.plot(time,np.rad2deg(elv))
    sub2.plot(time,np.rad2deg(Az))
    plt.show()
if __name__ == "__main__":
    test()
        