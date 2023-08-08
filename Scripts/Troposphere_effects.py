def Tropo_less_then_5(elv):
    import numpy as np
    from Atmosphere_attenuation import Apparent_elv
    if elv>np.deg2rad(5):
        ValueError("Elevation angle not below 5°")
        
    
    