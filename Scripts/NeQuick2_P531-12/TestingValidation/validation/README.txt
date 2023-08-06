===============================
README.txt:
Purpose: to describe format of the validation input/output files.

===============================
Input files for different stations distributed around the world for doy 100 of years 2000, 2004 and 2007 
(Folder "in")
Format:
*Filename: 
	<StationCode>_ne_in.dat

contents:

year  month Time(UTC hours)  lon(sta) lat(sta) heigth(sta) lon(sat) lat(sat) heigth(sat)

with lat, lon in degrees, height in meters and assuming geographic coordinates of WGS84 reference system.

===============================
Output reference associated files (Folder "out_ref")

Format:
*Filename:
        <StationCode>_ne_out_<Flux>.dat

contents:

year  month Time(UTC hours)  lon(sta) lat(sta) heigth(sta) lon(sat) lat(sat) heigth(sat) stec

with lat, lon in degrees, height in meters and assuming geographic coordinates of WGS84 reference system.
stec in TECU

