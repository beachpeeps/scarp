% process_paros.m
% The following code processes the paros pressure sensors from psi to
% hydrostatic head (m), correcting for atmospheric pressure. 
%
% Original code is from Adam Young.
% Updated 2020 for SCaRP experiment, Julia Fiedler jfiedler@ucsd.edu

reefbreakFolder = '/volumes/group/';
parosFolder = [ reefbreakFolder 'paros/TORREY_PINES/'];
parosName = dir([parosFolder '/P*']);
atmosFolder = [ reefbreakFolder 'paros/ATM-123630/'];

outFolder = '~/Documents/SCARP/data/processed/';

% Define constants used throughout all processing
rho = 1024; %1024 kg/m^3 seawater density
g = 9.81; %gravitational constant m/s^2
psi2pascal = 6894.76; % convert psi to pascals (kg/m/s^2)

% the following info is from
% reefbreak/group/DeploymentNotes2019_2020.xls
sampleFrequency = 2; %Hz, as denoted in the excel spreadsheet file
clockdrift = [-6 -1 0 -1 -5 2 -4 1 3 0]; % negative: paros is fast, positive: paros is slow
% note that 0 clock drift here means that it is unknown!!

%%

% loop through all paros
for parosNum = 2:10
    
    % get list of all hourly files on the server
    fileList = dir([parosFolder '/' parosName(parosNum).name '/*00.dat']);
    nfile = length(fileList);
    
    %% rename outfile to have year in front
    % CAUTION: THIS ONLY APPLIES TO TORREY WINTER DATA, NOT A UNIVERSAL FIX!!!
    outName = strings(nfile,1);
    for indHour=1:nfile
        fileName = fileList(indHour).name(1:8);
        if fileName(1) == '0'
            outName(indHour) = ['2020' fileName];
        elseif fileName(1) == '1'
            outName(indHour) = ['2019' fileName];
        end
    end
    
    daten = datenum(outName,'yyyymmddHHMM'); % change to datenum
    [dateSORT,indDate] = sort(daten);
    
    %% calculate clock drift correction over all time (add/subtract drift to/off total)
    cdrift = clockdrift(parosNum); % negative: paros is fast, positive: paros is slow
    
    hourDriftCorrection = cdrift/(nfile-1);  %%% amount of seconds to change to each file
    hourDriftCorrection = hourDriftCorrection/(60*60*24);   %%% amount of a day
    
    f = waitbar(0,'0%','Name',['Processing ' parosName(parosNum).name]);
    %% loop through all hourly files, correct for hourly clock drift, convert to depth
    for nhour=1:nfile
        parosStartTime = dateSORT(nhour);
        waitbar(100*nhour/nfile,f,sprintf('%3.1f',100*nhour/nfile))
        % get start time of file
        startTime = parosStartTime + (hourDriftCorrection*(nhour-1));
        
        % load in data
        psi = importdata([parosFolder parosName(parosNum).name '/' datestr(dateSORT(nhour),'mmddHHMM') '.dat']);
        
        % break loop if there's a file but no data to import
        if isempty(psi)
            break
        end
        
        % make time column, correcting for clock drift
        lenFile = length (psi (:,1));
        sampleIncrement = 1./(60*60*24*sampleFrequency);
        timeLength = lenFile * sampleIncrement;
        endTime = (startTime + timeLength)-sampleIncrement;
        pSensorTime = [startTime : sampleIncrement : endTime]';
        
        % remove atmospheric signal from same time period
        % TODO: Assumes clockdrift is not extreme!! Does not account
        % for clockdrift in atm data. Is this important??? Highly doubt
        % it...
        psi_atmos = importdata([atmosFolder datestr(dateSORT(nhour),'yyyy') '/' datestr(dateSORT(nhour),'mmddHHMM') '.dat']);
        psi = psi(:,1)-psi_atmos(:,1);
        
        % convert to meters of hydrostatic head --- press = rho g h
        % so h= press/(rho *g);
        h = psi.*psi2pascal./(rho * g);
        
        % check if folder exists, if not write it
        if ~exist([outFolder parosName(parosNum).name '/'],'dir')
            mkdir([outFolder parosName(parosNum).name '/'])
        end
        
        % name of hourly file
        saveName = strcat(outFolder, parosName(parosNum).name, '/', outName(nhour), '.mat');
        
        % save depth and time variables to hourly file, write to command window
        % that you've done it (if you want).
        save(saveName,'h','pSensorTime');
%         disp(strcat(saveName, ' saved'))
        clear psi h psi_atmos 
    end
end