import sgp4


class TLE_calc():
    def __init__(self,TLE_line_1,TLE_line_2):
        from datetime import datetime, timezone
        import sgp4
        self.TLE_line_1 = TLE_line_1
        self.TLE_line_2 = TLE_line_2
        #Get relevent time
        self.year_now = datetime.now(timezone.utc).year
        self.month_now = datetime.now(timezone.utc).month
        self.day_now = datetime.now(timezone.utc).day
        self.hr_now = datetime.now(timezone.utc).hour
        self.minute_now = datetime.now(timezone.utc).minute
        self.sec_now = datetime.now(timezone.utc).second
        self.jd, self.fr = sgp4.api.jday(self.year_now, self.month_now, self.day_now, self.hr_now,
                                         self.minute_now, self.sec_now)
        self.sat = sgp4.api.Satrec.twoline2rv(self.TLE_line_1,self.TLE_line_2)
        
        
        
    def get_speed(self):
        
        import numpy as np
            
        #Use SPG4 to get velocity 
        
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        speed = np.linalg.norm(v)
        return speed
    
    
    def get_pos(self):
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        return r
    
    def get_m_true(self):
        M_at_gather_tle = self.sat.mo
        
        m_change = self.sat.mdot
        
        year_tle = self.sat.epochyr
        tle_days_frac = self.sat.epochdays
        
        
        month_tle, day_tle, hour_tle, min_tle, sec_tle = sgp4.api.days2mdhms(year_tle, tle_days_frac)

        year_tle = year_tle+2000
        year_delta = self.year_now-year_tle
        
        if year_delta > 1:
            raise ValueError("TLE out of date")
    
        month_delta = self.month_now-month_tle
        
        if month_delta > 1:
            raise ValueError("TLE out of date")
        
        day_delta = (self.day_now-day_tle)*86400
        if day_delta > 24*86400:
            raise ValueError("TLE out of date")
        
        
        hour_delta = (self.hr_now-hour_tle)*3600

        min_delta = (self.minute_now-min_tle)*60
        
        second_delta = self.sec_now-sec_tle
        
        time_delta = (day_delta+hour_delta+min_delta+second_delta)/60
        
        m_delta  =m_change*time_delta
        
        return m_delta+M_at_gather_tle
        
    def get_true_anom(self, accuarcy):
        from scipy.special import jv as bassel_1st
        import numpy as np
        M = self.get_m_true()
        e = self.sat.ecco
        
        sum_term = 0
        for n in range(1,accuarcy+1):
            sum_term += (bassel_1st(n,n*e)/n)*np.sin(n*M)
        E = M+2*sum_term
        
        
        Q = np.sqrt((1+e)/(1-e))*np.tan(E/2)
        v = 2*np.arctan(Q)
        
        return v
            
        
    def get_ascension(self, accuarcy):
         import numpy as np
         i = self.sat.inclo
         omega = self.sat.argpo
         v = self.get_true_anom(accuarcy)
         Ohm = self.sat.nodeo
         
         alpha = (np.arctan(np.cos(i)*np.tan(omega+v)))-Ohm
         
         return alpha
         
         
         
    def get_declanation(self,accuarcy):
        import numpy as np
        i = self.sat.inclo
        omega = self.sat.argpo
        v = self.get_true_anom(accuarcy)
        
        delta = np.arcsin(np.sin(i)*np.sin(omega+v))
        return delta
        
    
    def test(self, accuarcy):
        print(self.sat.mo)


Iss_line_1 = '1 25544U 98067A   23195.46614675  .00010182  00000-0  18788-3 0  9996'
Iss_line_2 = '2 25544  51.6412 193.3466 0000287  80.3851 355.9960 15.49731575406021'

line_1 = '1 56226U 23084BK  23191.10910518  .00005567  00000+0  30866-3 0  9997'
line_2 = '2 56226  97.5132 307.3349 0012453 121.1638 239.0816 15.14117734  4126'


sat = TLE_calc(Iss_line_1,Iss_line_2)

acc = 10000
position = sat.get_pos()

v = sat.get_true_anom(acc)

alpha = sat.get_ascension(acc)

delta = sat.get_declanation(acc)




print(alpha,delta*180/3.141)
    



