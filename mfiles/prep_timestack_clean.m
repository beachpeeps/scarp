%prep timestacks

clear all
datadir = '/Volumes/FiedlerBot8000/scarp/';
savedir = [datadir '/mat/timestacks/'];
%% choose hover number and hover date and load data
hoverdate = '20191214';
% hoverdate = '20200224';
for hovern = 2:5




if str2double(hoverdate) == 20191214
    
    drone = load([datadir '/mat/lidar/drone/20191214_H' num2str(hovern) '_navd88_geoid12b_10cm.mat'],'Processed');
    truck = load([datadir '/mat/lidar/truck/20191214_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat'],'Processed');
    load('../mat/20191214/DroneStartStop_20191214.mat');
    load(['../mat/20191214_H' num2str(hovern) '_wiggle.mat']);
    xlimit = 1000;
    
elseif str2double(hoverdate) == 20200224
    
    drone = load([datadir '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat'],'Processed');
    truck = load([datadir '/mat/lidar/truck/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat'],'Processed');
    load('../mat/20200224/DroneStartStop_20200224.mat');
    load(['../mat/20200224_H' num2str(hovern) '_wiggle.mat']);
    xlimit = 500;
end
    

dewiggle = 1; %set to 0 if no de-wiggling, 1 if you want to de-wiggle
%% make sure the truck and drone are sampling from the same time
tchunk = [tstart(hovern) tstop(hovern)];

tind = knnsearch(truck.Processed.t',tchunk');
truck.Processed.t = truck.Processed.t(tind(1):tind(2));
truck.Processed.Zinterp = truck.Processed.Zinterp(tind(1):tind(2),:);
truck.Processed.Zinterp2 = truck.Processed.Zinterp2(tind(1):tind(2),:);

% put things into datetime format
tstart = datetime(tstart,'ConvertFrom','datenum');
tstop = datetime(tstop,'ConvertFrom','datenum');


%% interpolate lidar data onto the same time x grid
Hz_lidar = 10;
tvecHover = tstart(hovern):seconds(1/Hz_lidar):tstop(hovern);
tdrone = datetime(drone.Processed.t(2:end),'ConvertFrom','datenum');
[Tgrid,Xgrid] = meshgrid(tvecHover,drone.Processed.x);


% WIGGLE OR NO WIGGLE???
if dewiggle == true
    TSdrone = drone.Processed.Zinterp2-wiggle;
else
    TSdrone = drone.Processed.Zinterp2;
end

if sum(isnan(datenum(tdrone)))>0
    nanind = find(isnan(datenum(tdrone)));
    TSdrone(nanind,:) = [];
    tdrone(nanind) = [];
end

TSdrone = TSdrone(2:end,:);



TXdrone = interp2(datenum(tdrone),drone.Processed.x,TSdrone',datenum(Tgrid),Xgrid);
%manual clean up on Hover 1, Dec for marker pole
if hovern == 1 && str2double(hoverdate) == 20191214
rmInd = Xgrid>-78 & Xgrid<-76 & TXdrone>3.5;
TXdrone(rmInd) = nan;
end


TXdrone2 = inpaint_nans(TXdrone,2); %inpaint nans applies del^2 over whole matrix
%%
ttruck = truck.Processed.t;
ttruck = fillmissing(ttruck,'linear'); %there are weird missing times in the time vector??

TStruck = truck.Processed.Zinterp2;
TStruck(TStruck>nanmean(TStruck)+3*nanstd(TStruck))=NaN;
TStruck = TStruck';

%fixing errors in the time grid in Feb hovers - misinterpretation of
%"linescan" starts and stops due to noise leads to extra time grid spacing.
%Removing these blips makes the interpolation to the 10Hz time vector
%smoother.
if str2double(hoverdate) == 20200224
invalid = find(sum(isnan(TStruck))>=525);
TStruck(:,invalid) = [];
ttruck(invalid) = [];
end

TXtruck = interp2(ttruck,truck.Processed.x,TStruck,datenum(Tgrid),Xgrid,'nearest');
TXtruck2 = inpaint_nans(TXtruck,2);

%%
save([savedir hoverdate '_' num2str(hovern) '.mat'],'TXdrone','TXdrone2','TXtruck','TXtruck2','Hz_lidar','Tgrid','Xgrid','tvecHover')
disp(['saved data to ' savedir hoverdate '_' num2str(hovern) '.mat'])

end
