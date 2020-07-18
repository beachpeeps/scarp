function [tUTC,tLST] = GPStoUTC_midnightStart(tSecFromMidnight,filename,leapseconds)

% get the sunday start time for the filename
% returns also the time in local time for San Diego
% Needs some changes to account for local time in date of filename

% figure out filename date format
dstart = regexp(filename,'[0-9]{6}');

if isempty(dstart)
    disp('No time in the filename, please enter what day the data was collected in datetime format')
    prompt = 'Date in string format (datetime(yyyy,mm,dd))?)';
    filedate = input(prompt);
elseif dstart(1) == 1
    filedate = datetime(filename(1:8),'InputFormat','yyyyMMdd');
elseif dstart ~=1
    dd = [dstart; dstart+7]';
    dstr = filename([dd(1,1):dd(1,2) dd(2,1):dd(2,1)+5]);
    filedate = datetime(dstr,'InputFormat','yyyyMMddHHmmss');
end

% sunstart = dateshift(filedate,'dayofweek','sunday','previous');
midnightstart = dateshift(filedate,'start','day');


% tGPS is seconds since the sunday start time
secondsinoneday = 60*60*24;
midnightstart = datenum(midnightstart);
t_in_days = tSecFromMidnight./secondsinoneday;
tUTC = midnightstart+t_in_days;

% get local time
tUTCdt = datetime(tUTC(1),'ConvertFrom','datenum','TimeZone','Etc/UTC');
tLSTt = tUTCdt; 
tLSTt.TimeZone = 'America/Los_Angeles';

% get datenum of local time
tLST = tUTC + seconds(tzoffset(tLSTt))./secondsinoneday;



