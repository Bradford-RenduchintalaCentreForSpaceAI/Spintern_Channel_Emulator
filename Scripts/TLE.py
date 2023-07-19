


class TLE_calc():
    """
    TLE propergator that gathers data from given TLE values to begin intiate a sattelite using TLE_calc(line_1,line_2)
    
    Required packages numpy, datetime, sgo4, scipy, julian date
    
    It uses sgp4 to compute a lot of varibles and returns them ti relevent data 
    
        algorithms based on Chapter 1 of 
        %0 Book
        %T Doppler Applications in LEO Satellite Communication Systems
        %@ 1475783906
        %I Springer Publishing Company, Incorporated
        %A Irfan Ali
        %A Pierino G. Bonanni
        %A Naofal Al-Dhahir
        %A John E. Hershey
        %D 2013
    
    All data has been tested with current ISS data as shown in test() was 4° off annoyingly 
    """
    def __init__(self,TLE_line_1,TLE_line_2, use_time = False,year = None,month = None,day = None,hour = None,minute = None,sec = None,milli = None):
        """
        Parameters
        ----------
        TLE_line_1 : string
            First line of TLE data.
        TLE_line_2 : string
            Second line of TLE.

        Returns NONE
        -------
        """
        from sgp4.earth_gravity import wgs72
        from datetime import datetime, timezone
        from sgp4 import api
        from sgp4.io import twoline2rv as io2l
        self.TLE_line_1 = TLE_line_1
        self.TLE_line_2 = TLE_line_2
        if use_time == True:
            # Use given time
            self.year_now = self.year
            self.month_now = self.month
            self.day_now = self.day
            self.hr_now = self.hour
            self.minute_now = self.minute
            self.sec_now = self.second
            self.milli_now = self.milli
            self.jd, self.fr = api.jday(self.year_now, self.month_now, self.day_now, self.hr_now,
                                                 self.minute_now, self.sec_now)
        
        
        else:
            #Get Time now
            self.year_now = datetime.now(timezone.utc).year
            self.month_now = datetime.now(timezone.utc).month
            self.day_now = datetime.now(timezone.utc).day
            self.hr_now = datetime.now(timezone.utc).hour
            self.minute_now = datetime.now(timezone.utc).minute
            self.sec_now = datetime.now(timezone.utc).second
            self.milli_now = datetime.now(timezone.utc).microsecond
            self.jd, self.fr = api.jday(self.year_now, self.month_now, self.day_now, self.hr_now,
                                             self.minute_now, self.sec_now)
            self.sat = api.Satrec.twoline2rv(self.TLE_line_1,self.TLE_line_2)
        
        assert io2l(self.TLE_line_1, self.TLE_line_2,wgs72)
        
    def get_speed(self):
        """
        Returns
        -------
        speed : float
            Calculates the speed of the ISS by taking the magnitude of the velocity vector
            returned by sgp4.

        """
        import numpy as np
            
        #Use SPG4 to get velocity 
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        
        speed = np.linalg.norm(v) #is just sqrt(x^2+y^2+z^2)
        
        return {"speed":speed, "velocity":v}
    
    def get_pos_TEME(self):
        from sgp4.api import SGP4_ERRORS
        """

        Returns
        -------
        r : float array
            Returns the position array (x y z) of the sattellite in idiosyncratic True Equator 
            Mean Equinox coordinate frame.
        """
        e, r, v = self.sat.sgp4(self.jd,self.fr)
        if e != 0:
            print(SGP4_ERRORS[e])
        return {"x":r[0],"y":r[1],"z":r[2]}
    
    def get_m_true(self):
        """
        Raises
        ------
        ValueError
            Raises a value error if TLE data is older then 24 days.

        Returns
        -------
        TYPE float
            The new mean anomaly given the change in time from EPOCH of tle and current 
            time in UTC.

        """
        import sgp4.api
        M_at_gather_tle = self.sat.mo #Gather current mean anomaly
        
        m_change = self.sat.mdot #Gather the change in mean anomaly in rad/min
        
        #Gather TLE time data
        year_tle = self.sat.epochyr #Year 
        tle_days_frac = self.sat.epochdays # Days
        
        month_tle, day_tle, hour_tle, min_tle, sec_tle = sgp4.api.days2mdhms(year_tle, 
                                                                             tle_days_frac) #finds the other relevent data from TLE data

        year_tle = year_tle+2000 #TLE year is in format yy add 2000 to get to yyyy (goes out of date in 3000th century)
        
        year_delta = self.year_now-year_tle #Find the change in year 
        
        month_delta = self.month_now-month_tle #Find change in mounth
        
        day_delta = (self.day_now-day_tle)*86400 #Find day change in seconds
        
        if day_delta > (24*86400) or month_delta > 1 or year_delta>1: #Error if TLE data is out of date
            raise ValueError("TLE out of date")
        
        
        hour_delta = (self.hr_now-hour_tle)*3600 #Find hour change in seconds

        min_delta = (self.minute_now-min_tle)*60 #Find minute change in seconds
        
        second_delta = self.sec_now-sec_tle #Find seconds change
        
        time_delta = (day_delta+hour_delta+min_delta+second_delta)/60 # Get the time change in mins
        
        m_delta  =m_change*time_delta #Get the total change in mean anomaly by multiply change in min
        # by the momentum change
        
        return m_delta+M_at_gather_tle
        
    def get_true_anom(self, accuarcy):
        """

        Parameters
        ----------
        accuarcy : int
            The maximum count of the bassel coeffecient.

        Returns
        -------
        v : float
            true anomaly.

        """
        from scipy.special import jv as bassel_1st
        import numpy as np
        M = self.get_m_true() #Get current mean anomaly
        e = self.sat.ecco #Get inclination
        
        #Calculates the Eccentric Anomaly
        sum_term = 0 
        for n in range(1,accuarcy+1):
            sum_term += (bassel_1st(n,n*e)/n)*np.sin(n*M)
            
        E = M+2*sum_term
        
        #Caclulate true anomaly
        Q = np.sqrt((1+e)/(1-e))*np.tan(E/2)
        v = 2*np.arctan(Q)
        
        return v
            
    def get_ascension(self, accuarcy):
        """

        Parameters
        ----------
        accuarcy : int
            The maximum count of the bassel coeffecient.

        Returns
        -------
        alpha : float
             ascension in raduins.

        """
        import numpy as np
        i = self.sat.inclo #Get inclination
        omega = self.sat.argpo #Get Argument of perigee
        v = self.get_true_anom(accuarcy) #Get true annomanly
        Ohm = self.sat.nodeo #Get right ascending node
        
        alpha = (np.arctan(np.cos(i)*np.tan(omega+v)))-Ohm #calc acension
        
        return alpha
            
    def get_declanation(self,accuarcy):
        """

        Parameters
        ----------
        accuarcy : int
            The maximum count of the bassel coeffecient.

        Returns
        -------
        delta : float
            declanation in rads.

        """
        import numpy as np
        i = self.sat.inclo #Get inclination
        omega = self.sat.argpo # Get right ascending node
        v = self.get_true_anom(accuarcy) #Get true anomanly 
        
        delta = np.arcsin(np.sin(i)*np.sin(omega+v)) #Calc acension
        return delta

    def GMST(self):
        """

        Returns
        -------
        GMST : float
            Is the greenwich mean sideral time.
        Algorithm based off https://www.astrogreg.com/snippets/greenwichMeanSiderealTime.html

        """
        import numpy as np
        import juliandate as jd
        #Uses the julian date package to caculate the julian date
        jdq= jd.from_gregorian(self.year_now, self.month_now, self.day_now
                               ,self.hr_now,self.minute_now,self.sec_now,self.milli_now)
        
        jdq = self.jd+self.fr
    
        #finds centuries since j2000
        t = (jdq-2451545.0)/36525
        
        #finally calcs GMST
        
        temp = - 6.2e-6 * t * t * t + 0.093104 * t * t  + (876600.0 * 3600.0 + 8640184.812866) * t + 67310.54841
               
        
        return np.remainder(np.deg2rad(temp)/240,2*np.pi)

    def get_pos_ECEF(self,accuarcy):
        """


        Parameters
        ----------
        accuarcy : int
            The maximum count of the bassel coeffecient.

        Returns
        -------
        dict
            Positiion dictionary describing positi.

        """
        import numpy as np
        GMST = self.GMST()

        pos = self.get_pos_TEME()
        
        
        
       
        
        #Transforms position into the ECF frame
        x_bar = pos["x"]*np.cos(GMST)+pos["y"]*np.sin(GMST) 
        y_bar = -pos["x"]*np.sin(GMST)+pos["y"]*np.cos(GMST) 
        z_bar = pos["z"]
        
        long = np.arctan2(y_bar,x_bar)
        
        lat = np.arctan(z_bar/np.sqrt(x_bar**2+y_bar**2))
        return {"x":x_bar,"y":y_bar,"z":z_bar,"long":long,"lat":lat}
    
    def get_height(self,accuarcy):
        """
        

        Returns
        
        float
            Gives the height of sattelite in km.

        """
        import numpy as np
        pos = self.get_pos_ECEF(accuarcy) #get position
        
        return np.sqrt((pos["x"]**2)+(pos["y"]**2)+(pos["z"]**2))-6378.1 #Just computes sqrt(x^2+y^2+z^2)-earth raduis

    def get_slant_range(self, accuarcy,lat_of_ground,long_of_ground):
        """
        

        Parameters
        ----------
       accuarcy : int
           The maximum count of the bassel coeffecient.
        lat_of_ground : float
            lattitude of the ground station in degrees.
        long_of_ground : TYPE
            longitude of ground station in degrees.

        Returns
        -------
        d : float
            slant range.
        Az : flaot
            Azmith of sattellite in rad.
        theta : float
            eleveation angle in rad.

        """
        import numpy as np
        pi = np.pi # For personal taste
        
        #Convert lattitude to radians
        lat_of_ground = lat_of_ground*pi/180 
        long_of_ground = long_of_ground*pi/180

        delta = self.get_declanation(accuarcy)

        r = self.get_height(accuarcy)+6378.1 #Get height of sattelite above earth 
        
        R_e = 6378.1 #raduis of earth
        
        #Get long and lat of sattellite and convert to radians
        pos = self.get_pos_ECEF(accuarcy)
        lat_of_sat = pos["lat"]
        #Algorithm in the book
        
        H = -(lat_of_sat-lat_of_ground)
        
        cos_gamma = (np.sin(long_of_ground)*np.sin(delta))+(np.cos(long_of_ground)*np.cos(delta)*np.cos(H))
        
        gamma = np.arccos(cos_gamma)
        
        d = np.sqrt((R_e**2)+(r**2)-(2*R_e*cos_gamma))
        
        Az = (np.cos(delta)*np.sin(H))/np.sin(gamma)
        
        Az = np.arcsin(Az)
        
        theta = (r/d)*np.sin(gamma)
        
        theta = np.arccos(theta)
        
        
        return d, Az, theta
    def get_pass(self,accuarcy):
        import numpy as np
        
        pos = self.get_pos_ECEF(accuarcy)
        lat_of_sat = pos["lat"]
        long_of_sat = pos["long"]
        n = self.sat.no_kozai
        T = (2*np.pi)/n
        
        print(T)
        
       
        
        
        
    


        

def test():
    """
    Test file for TLE data that takes most current TLE data for ISS from https://live.ariss.org/tle/ 
    and return relevent parameters
    """
    import urllib.request
    import numpy as np
    # try:
    #     f = urllib.request.urlopen('https://live.ariss.org/iss.txt') #Gather TLE data
    #     url_text = f.read(200).decode('utf-8')
    #     url_text.split(" ")
    #     line_1_2 = url_text[13:] # Format the data to get both line 1 and 2 of TLE data 
    #     line_1 = line_1_2[0:71]
    #     line_2 = line_1_2[71:]
    # except:    
    line_1 = "1 25544U 98067A   23200.04569411  .00013707  00000-0  24898-3 0  9999"
    line_2 = '2 25544  51.6408 170.6667 0000390  75.0699   8.6204 15.49853169406738'
    sat = TLE_calc(line_1,line_2,False) # Initiate a sattellite using the TLE data provided 
    
    acc = 10000 # Set accuracy of computation to 100 see class TLE

    alpha = sat.get_ascension(acc) # Gather relevent data 
    delta = np.rad2deg(sat.get_declanation(acc))
    pos = sat.get_pos_ECEF(acc)
    x1,y1,z1 = sat.get_speed()["velocity"]
    height = sat.get_height(acc)
    d,Az, theta = sat.get_slant_range(acc, 51.0643, 0.8598) # Ground station at the location of Ham street Ashford Kent
    
    # print(f"""
    #       Lattitude:{np.rad2deg(pos["lat"])}, Longitude:{np.rad2deg(pos["long"])}
    #       Pos in ECEF (xyz): {pos["x"], pos["y"], pos["z"]}
    #       Slant length: {d}
    #       Azmuth: {np.rad2deg(Az)} 
    #       Elevation angle: {np.rad2deg(theta)}
    #       Height: {height}
    #       """)
    
    period = sat.get_pass(acc)



    
    
        
if __name__ == "__main__":
    test()
        
        
 

