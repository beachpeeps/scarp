function checkLAStimeOutput(filename)
c = lasdata(filename);
tGPS = get_gps_time(c);
[tUTC,tLST] = GPStoUTC(tGPS,filename);
disp(['Local start time is ' datestr(min(tLST)) '= ' datestr(min(tUTC)) ' UTC'])