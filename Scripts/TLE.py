class TLE_calc():
    def __init__(self,TLE_line_1,TLE_line_2):
        from datetime import datetime
        from sgp4.api import jday,Satrec
        self.TLE_line_1 = TLE_line_1
        self.TLE_line_2 = TLE_line_2
        #Get relevent time
        year = datetime.now().year
        month = datetime.now().month
        day = datetime.now().day
        hr = datetime.now().hour
        minute = datetime.now().minute
        sec = datetime.now().second
        self.jd, self.fr = jday(year, month, day, hr, minute, sec)
        self.sat = Satrec.twoline2rv(self.TLE_line_1,self.TLE_line_2)
        
        
        
    def get_speed(self):
        
        import numpy as np
            
        #Use SPG4 to get velocity 
        
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        speed = np.linalg.norm(v)
        return speed
    
    
    def get_pos(self):
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        return r
    
    
    


Iss_line_1 = '1 25544U 98067A   23192.44888453  .00009587  00000-0  17775-3 0  9990'
Iss_line_2 = '2 25544  51.6413 208.2874 0000213 101.0412  51.1708 15.49669613405559'

line_1 = '1 56226U 23084BK  23191.10910518  .00005567  00000+0  30866-3 0  9997'
line_2 = '2 56226  97.5132 307.3349 0012453 121.1638 239.0816 15.14117734  4126'


sat = TLE_calc(line_1,line_2)
print(sat.get_speed())

position = sat.get_pos()



data = []
import time
for i in range(0,100):
    data.append(sat.get_pos())
    time.sleep(1)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()


ax = fig.add_subplot(projection='3d')





plt.show()




