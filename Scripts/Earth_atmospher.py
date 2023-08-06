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
        raise(ValueError("Alltitude too high"))
   
    return T,p #T(k) #p(hPa)

def water_vapour_pressure(T,h,p_o):
    #from ITU-R P.835-6
    import numpy as np
    ro = p_o*np.exp(-h/2) #g/m^3
   
    e = (ro*T)/216.7 #hPa
   
    return e,ro
          
def Earth_atmosphere_model_test():
    import matplotlib.pyplot as plt
    import numpy as np
    from Earth_atmospher import Earth_atmosphere_model
   
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

if __name__ == "__main__":
    Earth_atmosphere_model_test()