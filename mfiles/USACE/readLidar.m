function [data, varargout] = readLidar(filenames, varargin)
%% Reads RXP files into matlab
%
%   Required inputs:
%       - filenames: a string or a cell array strings of paths to rxp
%                    file(s) or zip files that contain rxps. Zip files will
%                    be unzipped and any rxp files inside will be loaded.
%                    Unzipped files will be deleted once the data has been
%                    loaded
%
%   Optional inputs (passed in as name, value pairs):
%       - 'returns': which returns are read from file
%                    default = 'all'
%       - 'types':   cell array with type of point data to return
%                    default = {'x', 'y', 'z', 'r', 'tgps'}
%
%                    The following types can be used for various GPS data:
%                       - 'tgps': 
%                               Always tries to fill in gaps with GPS
%                               data. Will use filename datetime if
%                               necessary
%                       - 'tgps_datenum':
%                               Same as tgps but will convert values to
%                               datenums using filename
%                       - 'tint': 
%                               Lidar internal time
%                       - 'tfilename': 
%                               Generate gps time using filename and lidar 
%                               internal time
%                       - 'tgps_interp': 
%                               Interpolates between good gps data, but 
%                               will return nans if no good gps is
%                               available
%                       - 'tgps_interp_datenum':
%                               Same as tinterp but will convert values to
%                               datenums using filename
%                       - 'tgps_raw':
%                               Returns tgps from the scanner without
%                               trying to fill in nans
%                       - 'tgps_raw_datenum':
%                               Same as tgpsraw but will convert values to
%                               datenums using filename
%
%       - 'exe':     executable to use to read RXP file
%                    default = 'RXPconvert.exe'
%       - 'rmat':    rotation matrix applied to xyz values
%                    default = eye(4)
%       - 'dbg':     print debugging output by including any string with
%                    this option
%
%   Example uses:
%
%       readRXP('/path/to/file.rxp')
%       readRXP('file.rxp', 'data', {'x', 'y', 'z'})
%       readRXP({'file1.rxp', 'file2.rxp'}, 'exe', '/path/to/RXPconvert')
%       readRXP('file.rxp', 'returns', 'all', 'data', {'x', 'y', 'z'})
%
%




%% Parse input arguments
p = inputParser;

% Required inputs
addRequired(p, 'filenames', @isValidPath);

% Default values for optional inputs
defaultReturns = 'allpoints';
defaultData = {'x', 'y', 'z', 'r', 'tgps'};
defaultRmat = eye(4);
defaultExe = 'RXPConvert.exe';

% Optional inputs
addParameter(p, 'returns', defaultReturns, @isValidReturns);
addParameter(p, 'types', defaultData, @isValidTypes);
addParameter(p, 'rmat', defaultRmat, @isValidRMat);
addParameter(p, 'dbg', 'none', @ischar);
addParameter(p, 'exe', defaultExe, @isValidExe);

% Set parsing properties
p.KeepUnmatched = true;

% Parse
parse(p, filenames, varargin{:});

% Put into local variables for easier coding
returns = p.Results.returns;
types = p.Results.types;
exe = p.Results.exe;
rmat = p.Results.rmat;
isdbg = ~strcmp(p.Results.dbg, 'none');

% Debug
if isdbg
    disp(['File(s): ', filenames]);
    disp(['Returns: ', returns]);
    disp('Data types: ');
    disp(types);
    disp(['Executable: ', exe]);
    disp('Rotation Matrix: ');
    disp(rmat);
end


%% Build list of files to be loaded

files = {};     % List of rxp files to load data from
zipdirs = {};   % List of unzipped folders that will need to be deleted

if iscell(filenames)
    
    % Loop through all filenames, checking extensions for rxp and zip
    for f=1:numel(filenames)
        
        % Check for rxp or 3dd
        if strcmp(filenames{f}(end-2:end), 'rxp') || strcmp(filenames{f}(end-2:end), '3dd')
            files = [files filenames{f}];
        
        % Check for zip
        elseif strcmp(filenames{f}(end-2:end), 'zip')
            
            % Unzip the file
            folder = [filenames{f}(1:end-4) '-unzipped'];
            unzipped = unzip(filenames{f}, folder);
            zipdirs = [zipdirs folder];
            
            % Check for rxps or 3dds
            for i=1:numel(unzipped)
                if strcmp(unzipped{i}(end-2:end), 'rxp') || strcmp(unzipped{i}(end-2:end), '3dd')
                    files = [files unzipped{i}];
                end
            end
            
        end
            
        
    end
    
elseif strcmp(filenames(end-2:end), 'zip')
    
    % Only a zip file included, unzip it and look for rxps
    folder = [filenames(1:end-4) '-unzipped'];
    unzipped = unzip(filenames, folder);
    zipdirs = [zipdirs folder];
    
    % Check for rxps or 3dds
    for i=1:numel(unzipped)
        if strcmp(unzipped{i}(end-2:end), 'rxp') ||  strcmp(unzipped{i}(end-2:end), '3dd') 
            files = [files unzipped{i}];
        end
    end
    
elseif strcmp(filenames(end-2:end), 'rxp') || strcmp(filenames(end-2:end), '3dd')
    
    % A single rxp was included
    files = [files filenames];
    
end

%% Check that exe and file types are same
[~, exeName, ~] = fileparts(exe);
for f = files
    if ~strcmpi(exeName(1:3), f{:}(end-2:end))
        error(['Mismatched exe ' exeName(1:3) ' and file type ' ...
            f{:}(end-2:end) ', for file ' f{:} ]);
    end
end


%% Load data from file(s)

data = [];  % This will contain the full dataset

% Loop through each file
for i=1:numel(files)
    
    %% Load the data
    
    %   tgps_raw - get raw gps time
    %   tgps_raw_datenum - apply datenum to tgps_raw
    %   tgps_interp - interpolate tgps_raw
    %   tgps_interp_datenum - apply datenum to tgps_interp
    %   tgps - if tgps_interp still has all nans use tfilename
    %   tgps_datenum - apply datenum to tgps
    %
    %   tint - get raw internal time
    %   tfilename - use filename + tint to make tgps
    %   tfilename_datenum - apply datenum to tfilename
    
    % Build temporary filename to hold binary data
    [filepath, filename, ~] = fileparts(files{i});
    tmpfile = fullfile(filepath, [filename '.bin']);
    
    % Start building the command
    cmd = [exe ' ' files{i} ' ' tmpfile];
    
    % Check whether RXP for types and time data
    if strcmpi(exeName(1:3), 'RXP')
        cmd = [cmd ' ' returns];
    
        % Put gps options into categories
        tgps_types = {'tgps_raw', 'tgps_raw_datenum', 'tgps_interp', ...
                      'tgps_interp_datenum', 'tgps', 'tgps_datenum'};
        tint_types = {'tint', 'tfilename', 'tfilename_datenum'};

        % Add the data types
        for t=1:numel(types)
            if ismember(types{t}, tgps_types)
                cmd = [cmd ' tgps'];
            elseif ismember(types{t}, tint_types)
                cmd = [cmd ' tint'];
            else
                cmd = [cmd ' ' types{t}];
            end
        end
    elseif strcmpi(exeName(1:3), '3dd')
        types = {'x', 'y', 'z', 'r', 'tgps_raw'};
    end
    
    
    % Run command -- Linux needs LD_LIBRARY_PATH cleared from Matlab libs
    opSystemCmd = '';
    if ~ismember(pwd, 'C:\')
        opSystemCmd = 'env LD_LIBRARY_PATH='''' ';
    end
    [status, results] = system([opSystemCmd cmd]);
    
    % Check for errors
    if status==1
        disp(['Command: ' opSystemCmd cmd]);
        disp(results);
        error('File doesn''t exist or might be corrupt');
    end
    
    % Read data from temporary file
    if exist(tmpfile, 'file')
        fid = fopen(tmpfile);
        fdata = fread(fid, '*double');
        fclose(fid);
        delete(tmpfile);
        
    else
        
        error('File either doesn''t exist or can''t be converted');
        
    end
    
    %% Reshape the data
    nrows = numel(types);
    ncols = numel(fdata)/numel(types);
    fdata = transpose(reshape(fdata, nrows, ncols));
    
    
    %% Do GPS stuff
    % At this point, fdata is nReturns x mTypes
    % All time is in either tint or tgps_raw at this time
    
    % Logic:
    %   tgps_raw - time as reported by tGPS as sec since start of week
    %   tgps_raw_datenum - time as reported by tGPS as matlab datenum
    %   tgps_interp - interpolate using only tgps_raw
    %   tgps_interp_datenum - apply datenum to tgps_interp
    %   tgps - filter tgps_raw by filename and interpolate using tint
    %   tgps_datenum - apply datenum to tgps
    %
    %   tint - get raw internal time
    %   tfilename - use filename + tint to make tgps
    %   tfilename_datenum - apply datenum to tfilename
    
    % Check for any gps time data types and datenum conversions requested
    % by the user
    wantraw = any(~cellfun(@isempty, strfind(types, 'raw')));
    wantinterp = any(~cellfun(@isempty, strfind(types, 'interp')));
    wantgps = any(~cellfun(@isempty, strfind(types, 'tgps')));

    % If the user wants tgps_raw, the data is already good to go
    if wantraw  
        rawIndex = find(~cellfun(@isempty, strfind(types, 'raw')));
        for g = rawIndex
            if isdbg
                disp([types{g} ' data stored in column ' num2str(rawIndex)]);
            end
        end
        
    % If the user wants tgps_interp, the data will be interpolated using
    % raw only
    elseif wantinterp
        interpIndex = find(~cellfun(@isempty, strfind(types, 'interp')));
        
        for g = interpIndex
            % Count the number of bad gps points
            nanTime = isnan(fdata(:, g));
        
            % Some badGPS points
            if sum(nanTime) > 0
                % All pts bad -- cannot interpolate
                if sum(nanTime) >= size(fdata, 1)-1
                    if isdbg
                        disp([types{g} ' data stored in column ' num2str(g) ' is all NaNs']);
                    end
                else
                    if isdbg
                        disp('Some gps dropouts detected, interpolating with tgps_raw to fill gaps');
                    end
                    fdata(:, g) = interp1(find(~nanTime), fdata(~nanTime, g), ...
                        1:numel(nanTime), 'linear', 'extrap');
                end
            % All points are good    
            else
                if isdbg
                    disp(['All gps data in column ' num2str(g) ' is good']); 
                end

            end
        end
    
    % tgps checking
    elseif wantgps
        tgpsIndex = find(~cellfun(@isempty, strfind(types, 'tgps')));
        
        for g = tgpsIndex
            % Count the number of bad gps points and check time skip
            nanTime = isnan(fdata(:, g));
            maxTimeStepInSeconds = 1;
        
            % Some badGPS points
            
            if sum(nanTime) > 0 || any(abs(diff(fdata(:, g))) > maxTimeStepInSeconds)
                % All pts bad -- cannot interpolate
                if sum(nanTime) >= size(fdata, 1)-1
                    if isdbg
                        disp([types{g} ' data stored in column ' num2str(g) ' is all NaNs']);
                        if isempty(strfind(types{g}, 'datenum'))
                            types{g} = 'tfilename';
                        else
                            types{g} = 'tfilename_datenum';
                        end
                        disp(['Using ' types{g} ' values instead']);
                        fdata(:, g) = readLidar(files{i}, 'types', {'tint'}, 'exe', exe, 'returns', returns);
                        
                    end
                else
                    if isdbg
                        disp('Some gps dropouts or jumps detected, interpolating with good tgps_raw to fill gaps');
                        % Good determined as being within a scan length's
                        % time from the filename time (as gathered from
                        % tint and filename)
                    end
                    % Check to see if tint needed
                    maxTimeStepInSeconds = 1;
                    if any(abs(diff(fdata(:, g))) > maxTimeStepInSeconds)
                        % Read tint from file find the num of seconds long
                        % scan was according to the internal clock
                        %%%% Assume the internal clock never rests to 0 mid scan %%%%
                        tint = readLidar(files{i}, 'types', {'tint'}, 'exe', exe, 'returns', returns);
                        durationFromTint = range(tint); % as seconds
                        
                        % Determine which tgps points in fdata(:,g) are not
                        % at proper time by comparing to the start of week
                        fileDatenum = filename2Datenum(filepath);
                        weekStartFileDatenum = floor(fileDatenum) - (weekday(fileDatenum)-1); %datenum of current week of filename
                        fileDatenumAsSecOfWeek = (fileDatenum - weekStartFileDatenum)*24*60*60;
                        
                        goodGPSTime = abs(fdata(:, g) - fileDatenumAsSecOfWeek) < durationFromTint; %nans eval to false
                        
                        %also need to find the first contiguous set of good
                        %GPStime data
                        indJumps=find(diff(fdata(goodGPSTime, g))>1,1,'first');%assume the start of the file is good if there are jumps in the goodGPSTime section
                        
                        % Have some tgps times that are in the proper range
                        % and there are no jumps >1 sec in that good data
                        if sum(goodGPSTime) > 1 && isempty(indJumps)
                            goodtGPSMap = polyfit(find(goodGPSTime), fdata(goodGPSTime, g), 1);
                            fdata(:, g) = polyval(goodtGPSMap, 1:numel(goodGPSTime)); %in seconds
                        elseif sum(goodGPSTime) > 1 && ~isempty(indJumps)
                            %have good GPS data, but only using the first
                            %contiguous chunk for fit
                            goodtGPSMap = polyfit(find(goodGPSTime(1:indJumps)), fdata(goodGPSTime(1:indJumps), g), 1);
                            fdata(:, g) = polyval(goodtGPSMap, 1:numel(goodGPSTime)); %in seconds
                        else
                            % Have none, must use tint and filename
                            disp([types{g} ' data stored in column ' num2str(g) ' is all NaNs']);
                            if isempty(strfind(types{g}, 'datenum'))
                                types{g} = 'tfilename';
                            else
                                types{g} = 'tfilename_datenum';
                            end
                            disp(['Using ' types{g} ' values instead']);
                            fdata(:, g) = tint;
                        end
                    
                    % All tgps data is good, able to interpolate
                    else
                        fdata(:, g) = interp1(find(~nanTime), fdata(~nanTime, g), ...
                            1:numel(nanTime), 'linear', 'extrap');
                    end               
                    
                end
            % All points are good    
            else
                if isdbg
                    disp(['All gps data in column ' num2str(g) ' is good']); 
                end

            end
        end
    end %tgps control

    % User wants tint
    wanttint = any(~cellfun(@isempty, strfind(types, 'tint')));
    if wanttint
        tintIndex = find(~cellfun(@isempty, strfind(types, 'tint')));
        for g = tintIndex
            if isdbg
                disp([types{tintIndex} ' data stored in row ' num2str(tintIndex)]);
            end
        end
    end
    
    % The user either requested tfilename or tgps returned all nans
    wanttfilename = any(~cellfun(@isempty, strfind(types, 'filename')));
    if wanttfilename
        fnameIndex = find(~cellfun(@isempty, strfind(types, 'filename')));
        for g = fnameIndex
            % tint is in fdata(:,g), so we need to convert it to tgps sec since
            % week began
            fileDatenum = filename2Datenum(filepath);
            weekStartFileDatenum = floor(fileDatenum) - (weekday(fileDatenum)-1); %datenum of current week of filename
            fileDatenumAsSecOfWeek = (fileDatenum - weekStartFileDatenum)*24*60*60;
            fdata(:, g) = fileDatenumAsSecOfWeek + (fdata(:, g) - min(fdata(:, g)));
        end
    end %filename control
    
    % Convert gps to datenum if requested
    wantdatenum = any(~cellfun(@isempty, strfind(types, 'datenum')));
    if wantdatenum
        fileDatenum = filename2Datenum(filepath);
        datenumIndex = find(~cellfun(@isempty, strfind(types, 'datenum')));
        for g = datenumIndex
            if isdbg
                disp(['Converting ' types{g} ' data to datenum']);
            end
            fdata(:, g) = GPStime2datenum(fileDatenum, fdata(:, g));
        end
    end %datenum control
    
    
    %% Append the data to the full dataset
    data = [data; fdata];
    
end


%% Delete any temporary files/folders

for i=1:numel(zipdirs)
    rmdir(zipdirs{i}, 's');
end


%% Process data to build varargout

% Build scannum by looking at change in theta angles
%%%%%%%%%%%%%%%%%THIS IS A DUMB FIX. SCANNUM should be an input%%%%%%%%%%
%%%%%%%%%%%%%%%%%TRISTAN HELP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if nargout > 1 && ...
        p.Results.types{1}=='x' && p.Results.types{2}=='y' && p.Results.types{3}=='z'
    x_socs = data(:,1);
    y_socs = data(:,2);
    z_socs = data(:,3);
    
    theta = atan2(sqrt((x_socs).^2 + (y_socs).^2), z_socs) * 180/pi;
    scannum = nan(length(theta), 1);
    ind = [0; find(diff(theta) < -10); length(theta)];
    for i=1:length(ind)-1
        scannum(ind(i)+1:ind(i+1)) = i;
    end
    
    varargout{1} = scannum;
end


%% Make each point a row, then apply rotation matrix
if ~ismember('rmat', p.UsingDefaults)
    data = rotateXYZ(data, rmat);

end

end %fucntion close


function [good] = isValidPath(filenames)
% Checks for string or cell array of strings

    if iscell(filenames)
        
        good = true;

        % If a cell array was passed in, make sure each one is a string
        for i=1:numel(filenames)
            if ~ischar(filenames{i})
                good = false;
            end
        end

    else

        % Otherwise, just make sure a string was passed in
        good = ischar(filenames);

    end

end


function [good] = isValidReturns(returns)
% Checks for valid return string

    validreturns = {'allshots', 'allpoints', 'first', 'last'};

    good = any(validatestring(returns, validreturns));

end


function [good] = isValidTypes(data)
% Checks for valid type strings

    validTypes = {'a', 'b', 'd', 'R', 'f', 'h', 'p', 'e', 's', 'r', ...
                  'x', 'y', 'z', '0', '1', '2', '3', '4', '5', 'n', ...
                  'N', 'tgps_raw', 'tgps_raw_datenum', 'tgps_interp', ...
                  'tgps_interp_datenum', 'tgps', 'tgps_datenum', ...
                  'tint', 'tfilename', 'tfilename_datenum'};

    if iscell(data)
        
        good = true;
        
        % Make sure each data type is valid
        for i=1:numel(data)
            if ~any(validatestring(data{i}, validTypes))
                good = false;
                break;
            end
        end

    else

        good = false;

    end

end


function [good] = isValidExe(exe)
% Checks that executable exists

good = exist(exe, 'file');

end


function [good] = isValidRMat(rmat)
% Checks that the rotation matrix is 4x4

[m, n] = size(rmat);

good = ismatrix(rmat) && m==4 && n==4;

end


function [tgpsdatenum] = GPStime2datenum(filenamedatenum, tgps)
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
