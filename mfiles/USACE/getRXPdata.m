function data=getRXPdata(filename,pts2return,varargin)
%% Reads RXP files into matlab
%Example Syntax
% data=getRXPdata('test.rxp','allpoints','x','y','z','tgps','a');
% data=getRXPdata('test.zip',{'last',100},'x','y','z','r');
% data=getRXPdata({'test.rxp','test2.rxp','test3.rxp'},'allpoints','x','y','z');
% data=getRXPdata({'test.zip','frame'},'allpoints','x','y','z');
%
%Output
% data(mxn)
%     m is the number of inputs into varargin
%     n is the number of points returned
% ZipNames
%     if it's a zipfile output the files used
%Inputs
% filename
%     Option 1: a string of the filename
%     Option 2: a cell with strings of the filenames (all data combined
%     into data)
%     Option 3: a zip file as a string with RXPs inside of it
%           3b: a cell with a zip file as cell{1} and the fileID parameter as cell{2}
% pts2return
% 		 Option Arguments are:
% 		   'allshots' returns a point even when no shots returned
% 			 useful to return the Beam Direction with 'allshots'
% 		   'allpoints' returns every point that was returned
% 		   'first' returns only the first points returned
% 		   'last' returns only the last points returned
% 			 ** if only one point is returned it is will be written
% 			 ** with both the 'first' and 'last' tag
%        To return only a subset of the points, pass the parameters in as a
%        cell, where cell{2} is the number of points
% varargin
% 		 Optional Arguments are:
% 		  'a': amplitude (dB)
% 		  'b': background_radiation
% 		  'd': deviation
% 		  'R': echo_range
% 		  'f': facet number
% 		  'h': is_high_power
% 		  'p': is_pps_locked
% 		  'e': is_rising_edge
% 		  's': is_sw
% 		  'r': reflectance (dB)
% 		  'tgps': GPS time (can't be used with 'tint')
% 		  'tint': Internal time (can't be used with 'tgps')
% 		  'x': x position of point in the SOCS (m)
% 		  'y': y position of point in the SOCS (m)
% 		  'z': z position of point in the SOCS (m)
% 		  '0': beam_origin[0] (x distance from phase center) (m)
% 		  '1': beam_origin[1] (y distance from phase center) (m)
% 		  '2': beam_origin[2] (z distance from phase center) (m)
% 		  '3': beam_direction[0] (x component of unit direction vector)
% 		  '4': beam_direction[1] (y component of unit direction vector)
% 		  '5': beam_direction[2] (z component of unit direction vector)
% 		  'n': number of return
% 		  'N': number of returns for given pulse
%
% 		 *NOTE.  if using 'beam_direction' and 'range', you must take
% 		 * into account the 'beam_origin' offset to get to the SOCS
%
% 		 *SOCS = Scanners Own Coordinate System
if iscell(filename) && strcmp(filename{1}(end-2:end),'zip') %zip file with ID
    fname=filename{1};
    fnameid=filename{2};
    unzip(fname,'tempZipRXP');
    fnames=dir(['tempZipRXP\*' fnameid '*.rxp']);
    if numel(fnames)>0
        for i=1:numel(fnames)
            ZipNames{i}=['tempZipRXP\' fnames(i).name];
        end
        data=getRXPdat(ZipNames{1},pts2return,varargin);
        for i=2:numel(ZipNames)
            data=[data getRXPdat(ZipNames{i},pts2return,varargin)];
        end
    else
        warning(['No ' fnameid ' data in ' fname]);
        data=nan;
    end
    rmdir('tempZipRXP','s');
elseif iscell(filename) %cell of filenames (not zipped)
    data=getRXPdat(filename{1},pts2return,varargin);
    for i=2:numel(filename)
        data=[data getRXPdat(filename{i},pts2return,varargin)];
    end
elseif strcmp(filename(end-2:end),'zip') %zip file without ID
    fname=filename;
    fnameid=[];
    unzip(fname,'tempZipRXP');
    fnames=dir(['tempZipRXP\*' fnameid '*.rxp']);
    if numel(fnames)>0
        for i=1:numel(fnames)
            ZipNames{i}=['tempZipRXP\' fnames(i).name];
        end
        data=getRXPdat(ZipNames{1},pts2return,varargin);
        for i=2:numel(ZipNames)
            data=[data getRXPdat(ZipNames{i},pts2return,varargin)];
        end
    else
        warning(['No ' fnameid ' data in ' fname]);
        data=nan;
    end
    rmdir('tempZipRXP','s');
else %single filename
    data=getRXPdat(filename,pts2return,varargin);
end
end
function data=getRXPdat(filename,pts2return,varargin)
% example code
% data=getRXPdata('testscan.rxp','allpoints','x','y','z','tgps');

if numel(varargin)==1 && iscell(varargin{1})
    varargin=varargin{1};
end

% This is awful
if ~isempty(varargin) && ~isempty(strfind(varargin{length(varargin)}, 'RXPconvert'))
    exe = varargin{length(varargin)};
    varargin(length(varargin)) = [];
else
    exe = 'RXPconvert.exe';
end


if iscell(pts2return)
    Num2return=pts2return{2};
    pts2return=pts2return{1};
else
    Num2return=nan;
end

[~,name,~] = fileparts(filename);

tempname=['temp_' name '.bin'];

fid=fopen(tempname,'w+t');
fclose(fid);

params=[exe ' ' filename ' ' tempname ' ' pts2return] ;
doGPS2datenum=false;
doGPS2datenumNum=nan;
for i=1:numel(varargin)
    if strcmp(varargin{i},'tgps_datenum')
        varargin{i}='tgps';
        doGPS2datenum=true;
        doGPS2datenumNum=i;
    end
    params=[params ' ' varargin{i}];
end
fprintf('%s\n',params);
fprintf('writing data to temp file\n');

[status, results]=system(params);

if status==1
    results
    error('File doesnt exist or might be corrupt');
end
fprintf('reading data from temp file\n');
if exist(tempname,'file')~=0
    fid = fopen(tempname);
    if isnan(Num2return)
        data=fread(fid,'*double');
    else
        data=fread(fid,Num2return*numel(varargin),'*double');
    end
    fclose(fid);
else
    error('File either doesnt exist of cant be converted');
end
data=reshape(data,numel(varargin),numel(data)/numel(varargin));


if doGPS2datenum
    [~,fname,~]=fileparts(filename);
    fnametime=datenum(fname(1:16),'yyyymmdd-HHMM-ss');
    data(doGPS2datenumNum,:)=GPStime2datenum(fnametime,data(doGPS2datenumNum,:));
end

delete(tempname);

end
function tgpsdatenum=GPStime2datenum(filenamedatenum,tgps)
%% tgps=datenum(now)-(ceil(datenum(now))-weekday(now));

%convert tgps to datenum (days)
tgps=tgps/60/60/24;
%convert filenamedatenum to gps time
filetimefromweekstart=filenamedatenum-(ceil(filenamedatenum)-weekday(filenamedatenum));
%check for huge discrepencies of greater than 12 hours
% but also not greater than 6.25 days. (that is a wrap around end of week
% thing)
timeDiff=abs(nanmean(tgps(:))-filetimefromweekstart)>0.5;
if timeDiff>0.5 && timeDiff<6.25
   tgps=tgps+floor(filetimefromweekstart); 
end

%calculate the datenum of the start of the filename week
timeofweekstartcurrent=ceil(filenamedatenum)-weekday(filenamedatenum);
%calculate datenum assuming middle of week
tgpsdatenum=timeofweekstartcurrent+tgps;

%Case 1: gps is ahead of filename time by a max of 3.5 days.  If filetime is
%a saturday and gps time is a sunday, this correction is needed. Otherwise
%the tofweekstart is 7 days backwards in time
tgpsdatenum(tgps<1.75 & filetimefromweekstart>5.25)=tgpsdatenum(tgps<1.75 & filetimefromweekstart>5.25)+7;

%Case 2: filename time is ahead of gps time by a max of 3.5 days.  If filetime is
%a sunday and gps time is a saturday, this correction is needed. Otherwise
%the tofweekstart is 7 days ahead in time
tgpsdatenum(tgps>5.25 & filetimefromweekstart<1.75)=tgpsdatenum(tgps>5.25 & filetimefromweekstart<1.75)-7;

% If the Clock drifts by more than 1.75 days in either direction... there
% will be issues.  If this starts happening, code can be written smarter so
% that it can be good within 3.5 days.  or some other offsets can be used.

end