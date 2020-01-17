function h = plot_lidarData(nfile,cvar,filedir, filename,THETA,XO,YO)


c = lasdata([filedir filename(nfile).name]);
tSecFromSunday = get_gps_time(c);

[tUTC,tLST] = GPStoUTC(tSecFromSunday,filename(nfile).name);
[T, x, y, z, amp, sortedInd] = sortLASobject(tUTC,c);
% rotate data
[XR, YR] = xyRotate(x,y,THETA,XO,YO);
%
ind = 1:1000:length(z);
hold on
h = scatter3(XR(ind),YR(ind),z(ind),10,cvar);