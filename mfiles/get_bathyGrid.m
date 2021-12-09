%make bathymetry grid for SCaRP data
%this is all very specific so I'm not sure it even belongs in this repo
clear
% the goal is to get a grid of Z values for time (continuous) and x,
% where Z(x,t)
MopID = 582;
dx = 1;
scarpdates = [datenum(2019,12,06,0,0,0) datenum(2020,03,04,0,0,0)];
t = scarpdates(1):1:scarpdates(2);
% 2021, Julia Fiedler

addpath /Volumes/group/MOPS/toolbox
addpath /Volumes/group/MOPS

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% get transect
n=MopNumber;
load(['M' num2str(n,'%05.0f') 'SM.mat'])

x = SM(1).X1D(1):dx:1200; %this is our crosshore distance, in meters
T =struct2table(SM);

% get SCARP profiles
d = T.Datenum;
d_ind = find(d>=scarpdates(1) & d<=scarpdates(2));
surveyType = categorical(T.Source);
surveyType = surveyType(d_ind);


%get exact dates and times OH MY GOD I HATE THIS
dateScarp = T.File(d_ind);

tScarp = nan(size(dateScarp));
for i= 1:length(dateScarp)
    a = char(dateScarp(i));
    if surveyType(i) == 'Trk'
        strdate =  a([59:66 80:83]);
        if sum(isletter(strdate)) == 0
            tScarp(i) = datenum([strdate, '00'],'yyyymmddHHMMSS');
        elseif sum(isletter(strdate)) == 4
            tScarp(i) = datenum([strdate(1:8), '000000'],'yyyymmddHHMMSS');
        end
    elseif surveyType(i) == 'Gps-Jumbo' || surveyType(i) == 'Gps'
        tScarp(i) = datenum([a(26:33), '000000'],'yyyymmddHHMMSS');
    end
end

%kluge to change time based on tides on truck survey
tScarp(8) = tScarp(8) + 15/24;

% check if offshore data
% also extend z out to ~15m or more
Z = T.Z1Dmean;
Z = Z(d_ind,:); %this is our bathy for SCaRP dates, 
% but it only extends to ~10m or so, and some data is missing

%load the survey files with the jetski tracks and rotate into X1D world
load(['M' num2str(n,'%05.0f') 'SA.mat'])
for i=1:length(d_ind)
    [Nmop,X]= UTM2MopxshoreX(SA(d_ind(i)).X,SA(d_ind(i)).Y);
    % now average the survey X,z values into unique 1m xshore X distance bins
    [ux, ~, xidx] = unique(X); % ux are unique 1m resolution offshore distances
    zavg=accumarray(xidx(:), SA(d_ind(i)).Z',[], @mean); % mean z for each ux
    X1D=min(ux):max(ux);
    Z1D=interp1(ux,zavg,X1D);
    Zmop(i,:) = interp1(X1D,Z1D,x);
end

ZZ = Zmop;
ZBill = [Z nan(length(d_ind),length(Zmop)-length(Z))]; %this extends the SM Z bathy to our current x-domain

k = ~isnan(ZBill); %find all times when Bill has made a fancy gridded product
ZZ(k) = ZBill(k); %now we have something gridded in x that includes offshore

%but we want something gridded in t, so:
[xx, ttScarp] = meshgrid(x,tScarp);
[xxnew,tt] = meshgrid(x,t);
ZZt = interp2(xx,ttScarp,ZZ,xxnew,tt); %interpolate it to a new grid (using prescribed time vector)

%make a new matrix that shows time elapsed since previous survey
ZZ_daysElapsed = ZZt;
ZZ_daysElapsed = double(~isnan(ZZ_daysElapsed));

for i = 1:length(x)
    daysElapsed(:,i) = get_daysElapsed(ZZ_daysElapsed(:,i));
end





%now we linearly interpolate in time to fill in the gaps
%this comes with lots of caveats but for now it is good enough
ZZfill = fillmissing(ZZt,'linear',1);
pcolor(x,t,ZZfill); shading flat

%%
for i=1:length(t)
plot(x,ZZfill(i,:))
xlim([0 300])
ylim([-10 5])
pause(0.1)
end

save('../mat/scarp_bathy.mat','ZZfill','x','t','daysElapsed','ZZt')


