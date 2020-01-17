function [tUTC,tLST] = GPStoUTC(tSecFromSunday,filename,leapseconds)

% get the sunday start time for the filename

% figure out filename date format
dstart = regexp(filename,'[0-9]{6}');

if dstart == 1
    filedate = datetime(filename(1:8),'InputFormat','yyyyMMdd');
elseif dstart ~=1
    dd = [dstart; dstart+5]';
    dstr = filename([dd(1,1):dd(1,2) dd(2,1):dd(2,2)]);
    filedate = datetime(dstr,'InputFormat','yyyyMMddHHmmss');
end

sunstart = dateshift(filedate,'dayofweek','sunday','previous');

% account for leapseconds
if nargin<3 && year(sunstart) == 2019
    leapseconds = 18;
elseif nargin<3 && year(sunstart) ~=2019
    disp('Please enter in leapseconds for your time period')
end

% tGPS is seconds since the sunday start time
secondsinoneday = 60*60*24;
sunstart = datenum(sunstart);
t_in_days = (tSecFromSunday-leapseconds)./secondsinoneday;
tUTC = sunstart+t_in_days;

% get local time
tUTCdt = datetime(tUTC(1),'ConvertFrom','datenum','TimeZone','Etc/UTC');
tLSTt = tUTCdt; 
tLSTt.TimeZone = 'America/Los_Angeles';

% get datenum of local time
tLST = tUTC + seconds(tzoffset(tLSTt))./secondsinoneday;



