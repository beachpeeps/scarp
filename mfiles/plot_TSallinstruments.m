%plot_TSallinstruments
% script to compare the lidar with the paros
% will spit out 3 figures: timeseries, spectra, and scatter plots
clear
%choose hover number and paros number and load data
hovern = 2;
parosn = 2;
datadir = '/Volumes/FiedlerBot8000/scarp/';
drone = load([datadir '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat'],'Processed');
truck = load([datadir '/mat/lidar/truck/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat'],'Processed');
load('../mat/paros.mat');
load('../mat/20200224/DroneStartStop_20200224.mat');
load(['../mat/20200224_H' num2str(hovern) '_wiggle.mat'])
%%
% make sure the truck and drone are sampling from the same time
tchunk = [tstart(hovern) tstop(hovern)];

tind = knnsearch(truck.Processed.t',tchunk');
truck.Processed.t = truck.Processed.t(tind(1):tind(2));
truck.Processed.Zinterp = truck.Processed.Zinterp(tind(1):tind(2),:);
truck.Processed.Zinterp2 = truck.Processed.Zinterp2(tind(1):tind(2),:);

% put things into datetime format
tstart = datetime(tstart,'ConvertFrom','datenum');
tstop = datetime(tstop,'ConvertFrom','datenum');

%% get paros data

% find which hour to pull from
tstarthour = dateshift(tstart,'start','hour');
tind = find(paros(parosn).t == tstarthour(hovern));

% check if hover spans 2 hours
tstophour = dateshift(tstop,'start','hour');
if tstophour(hovern)>tstarthour(hovern)
    tind = [tind tind+1];
end

% make a time vector for the paros hourly chunk
for i=1:length(tind)
tvec(:,i) = paros(parosn).t(tind(i)):seconds(0.5):paros(parosn).t(tind(i))+7166*seconds(0.5);
end
tvec = tvec(:);

%NOTE: THERE IS A 1 SECOND OFFSET in the paros data at P4!!! 
if parosn == 4
    tvec = tvec-seconds(1);
end
%NOTE: THERE IS A 0.5 SECOND OFFSET in the paros data at P4!!! 

if parosn == 2
    tvec = tvec-seconds(0.5);
end

% find location in gridded lidar data that matches paros location
xind = knnsearch(truck.Processed.x',paros(parosn).crossshore);

%% interpolate lidar data onto the same time grid
Hz_lidar = 10;
Hz_paros = 2;
tvecHover = tstart(hovern):seconds(1/Hz_lidar):tstop(hovern);


TSdrone = drone.Processed.Zinterp2(:,xind)-wiggle;
TSdrone = interp1(datetime(drone.Processed.t(~isnan(TSdrone)),'ConvertFrom','datenum'),TSdrone(~isnan(TSdrone)),tvecHover);

% note: check if this can be Zinterp2 or Zinterp!!!
% Things look a lot better visually if using Zinterp, but meh.
% using linear interpolation here to put things on to same time grid -
% CAUTION!
TStruck = truck.Processed.Zinterp2(:,xind);
%debird and dehuman
TStruck(TStruck>nanmean(TStruck)+3*nanstd(TStruck))=NaN;
% linearly interpolate to the same time grid as the the drone
TStruck = interp1(datetime(truck.Processed.t(~isnan(TStruck)),'ConvertFrom','datenum'),TStruck(~isnan(TStruck)),tvecHover);

%infill the nans in the TStruck timeseries if using Zinterp2 (CAUTION)!!!
TStruck = inpaint_nans(TStruck); 

% just make sure that the start/stop times for the paros are the same as
% lidar
tvecHover2Hz = tstart(hovern):seconds(1/Hz_paros):tstop(hovern);
if ~isempty(paros(parosn).offset)
pvec = paros(parosn).etaCorrected(:,tind(:))+paros(parosn).z-paros(parosn).offset;
else
    pvec = paros(parosn).etaCorrected(:,tind(:))+paros(parosn).z;
end
TSparos = interp1(tvec,pvec(:),tvecHover2Hz);


nfftL = 5*60*Hz_lidar; % 5 minute windows
nfftP = 5*60*Hz_paros; % 5 minute windows

% get spectra for each instrument (unfiltered)
[fmp, SppParos, ~, ~, nens, dof] = get_spectrum(TSparos, nfftP, Hz_paros, 0.05);
[fm, SppDrone, ~, ~, nens, dof] = get_spectrum(TSdrone, nfftL, Hz_lidar, 0.05);
[fm, SppTruck, ~, ~, nens, dof] = get_spectrum(TStruck, nfftL, Hz_lidar, 0.05);

%% ******prewhiten the timeseries to avoid ringing at beginning/end********
sizepad = 60*Hz_lidar; %60 second padding - should be enough??
cutoff = 0.25; %Hz cutoff to match paros cutoff

pad = ones(sizepad,1);
pre = nanmean(TStruck(1:5*Hz_lidar))*pad;
post = nanmean(TStruck(end-5*Hz_lidar:end))*pad;
ppad = [pre; TStruck'; post];

T = lowpass(ppad,cutoff,Hz_lidar); %using matlab's standard here...
T = T(sizepad+1:end-sizepad); %remove padding

pre = nanmean(TSdrone(1:5*Hz_lidar))*pad;
post = nanmean(TSdrone(end-5*Hz_lidar:end))*pad;
ppad = [pre; TSdrone'; post];

D = lowpass(ppad,cutoff,Hz_lidar); %using matlab's standard here...
D = D(sizepad+1:end-sizepad); %remove padding

% get spectra for filtered lidar timeseries
[fm, SppT, ~, ~, nens, dof] = get_spectrum(T, nfftL, Hz_lidar, 0.05);
[fm, SppD, ~, ~, nens, dof] = get_spectrum(D, nfftL, Hz_lidar, 0.05);


%% make timeseries figure
figure
plot(tvecHover2Hz,TSparos,'k')
hold on
plot(tvecHover,TSdrone,':')
plot(tvecHover,TStruck,':')
xlim([tstart(hovern) tstart(hovern)+minutes(8)])
plot(tvecHover,D,'b')
plot(tvecHover,T,'r')
legend('paros','drone','truck','droneFilt','truckFilt')
ylabel('z NAVD (m)')

%% make spectra figure
figure
loglog(fmp,SppParos)
hold on
loglog(fm,SppDrone)
loglog(fm,SppTruck)
loglog(fm,SppD)
loglog(fm,SppT)

legend('paros','drone','truck','droneFilt','truckFilt')
xlabel('Frequency (Hz)')
ylabel('Energy (m^2/Hz)')
title(['Lidar data filtered at ' num2str(cutoff) ' Hz'])

%% make scatter plot
% linearly interpolate onto same time grid for plotting porpoises
T_2Hz = interp1(tvecHover,T,tvecHover2Hz);
D_2Hz = interp1(tvecHover,D,tvecHover2Hz);

lidar_data = [D_2Hz; T_2Hz];
lidar_label = {'drone (m) @ 2Hz','truck (m) @ 2Hz'};

figure

for i=1:2
    ax(i) = subplot(1,2,i);
    hs(i) = scatter(TSparos,lidar_data(i,:),'o','filled');
    hs(i).MarkerFaceAlpha = 0.2;
    hold on
    
    bias(i) = mean(lidar_data(i,:)-TSparos);
    [rmse(i),Sk(i)] = get_Skill(lidar_data(i,:)-bias(i),TSparos);
    
    ax(i).XLabel.String = 'paros (m)';
    ax(i).XLim(1) = min([TSparos(:);D_2Hz(:);T_2Hz(:)])-0.1;
    ax(i).XLim(2) = max([TSparos(:);D_2Hz(:);T_2Hz(:)])+0.1;
    ax(i).YLim = ax(i).XLim;

    plot(ax(i).XLim,ax(i).XLim,':k','parent',ax(i))
    ax(i).YLabel.String = lidar_label(i);
    
    box on
    
    text(0.95,0.1,...
        {['bias = ' num2str(100*bias(i),'%2.2f') 'cm'],...
        ['rmse = ' num2str(100*rmse(i),'%2.2f') 'cm'],...
        ['Sk = ' num2str(Sk(i),'%2.2f')]},'units','normalized',...
        'horizontalalignment','right')
end

suptitle(['Hover ' num2str(hovern) ' : P' num2str(parosn)])
