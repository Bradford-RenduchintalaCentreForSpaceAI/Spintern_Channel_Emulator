def test():
    from TLE import TLE_calc
    import numpy as np
    import time
    line_1 = "1 25544U 98067A   23200.04569411  .00013707  00000-0  24898-3 0  9999"
    line_2 = '2 25544  51.6408 170.6667 0000390  75.0699   8.6204 15.49853169406738'

    sat = TLE_calc(line_1, line_2, True)
    
    print(sat)
    


if __name__ == "__main__":
    test()
        