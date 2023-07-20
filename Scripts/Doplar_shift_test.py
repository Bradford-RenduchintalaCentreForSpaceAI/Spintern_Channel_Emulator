def Get_orbital():
    from TLE import TLE_calc
    import numpy as np
    import datetime
    from datetime import timezone
    line_1 = "1 25544U 98067A   23201.19134634  .00010982  00000-0  20095-3 0  9998"
    line_2 = '2 25544  51.6392 164.9972 0000314  45.9683 314.1332 15.49879208406919'
    lat = 51.0643
    long = 0.8598
    t_delta = 1
    acc = 100
    now = datetime.datetime.now(timezone.utc) 

    sat = TLE_calc(line_1, line_2,True,now.year,now.month,now.day,now.hour,now.minute,
                  now.second,now.microsecond)
    
    d,Az, elv,q = sat.get_pos_TOPO(acc,lat,long)
    sec = 0
    while np.rad2deg(Az) and np.rad2deg(elv) < 0:
        now = datetime.datetime.now(timezone.utc)+datetime.timedelta(0,sec)
        sat = TLE_calc(line_1, line_2,True,now.year,now.month,now.day,now.hour,now.minute,
                      now.second,now.microsecond)
        d,Az, elv,q = sat.get_pos_TOPO(acc,lat,long)
        sec = sec+1
        
    return datetime.timedelta(seconds=sec)

def test():
    print(f"Time till pass {Get_orbital()}")



if __name__ == "__main__":
    test()
        