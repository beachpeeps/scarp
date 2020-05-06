clear
hFig = figure;
hoverN = 2;

addpath /Users/juliafiedler/Documents/SWASH_runup/mfiles/functions/
vizdir = '../viz/';

drone = load(['../mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hoverN,'%d') '.mat']);

truckdir = '../mat/lidar/truck/';
truckscans = dir([truckdir 'WaveScanDrone*.mat']);
% truckscans = dir([truckdir 'RiProcess*.mat']);

truck = load([truckdir truckscans(hoverN-1).name]);
% truck = load([truckdir truckscans(2).name]);

load('../mat/paros.mat','paros')

minTime = datenum(dateshift(datetime(nanmin(drone.Processed.t),'ConvertFrom','datenum'),'start','minute'));
maxTime = datenum(dateshift(datetime(nanmax(drone.Processed.t),'ConvertFrom','datenum'),'end','minute'));

[t_verified,t_predicted,verified,predicted,tideInfo] = getNOAAtide(minTime,maxTime,'9410230');

indTide = find(t_verified>minTime & t_verified<=maxTime);

tideMean = mean(verified(indTide));
%%
for pnum = 2:10
indHour = find(datenum(paros(pnum).t) == datenum(dateshift(datetime(minTime,'ConvertFrom','datenum'),'start','hour')));

tchunk = round((maxTime-minTime)*60*60*24*2);
tchunkStart = round((minTime-datenum(paros(pnum).t(indHour)))*60*60*24*2);
if tchunkStart==0
    tchunkStart = 1;
end

fitinHour = tchunkStart+tchunk-length(paros(pnum).eta);

if fitinHour>0
    pchunk(pnum,:) = paros(pnum).etaCorrected(tchunkStart:end,indHour);
elseif fitinHour<0
    pchunk(pnum,:) = paros(pnum).etaCorrected(tchunkStart:tchunkStart+tchunk,indHour);
end
end

%%
figure

for pnum = 2:10
    if ~isempty(paros(pnum).offset)
        zloc = mean(pchunk(pnum,:))+paros(pnum).z-paros(pnum).offset;
    else
        zloc = mean(pchunk(pnum,:))+paros(pnum).z;
    end
    
    scatter(paros(pnum).crossshore,zloc)
    hold on
    text(paros(pnum).crossshore,zloc,paros(pnum).name(1:3))
end
%%
meanDrone = nanmean(drone.Processed.Zinterp2);
perGood = sum(~isnan(drone.Processed.Zinterp2))./length(drone.Processed.Zinterp2);
meanDrone(perGood<0.5) = NaN;

minDrone = nanmin(drone.Processed.Zinterp2);
minDrone(perGood<0.5) = NaN;




meanTruck = nanmean(truck.Processed.Zinterp2);
perGood = sum(~isnan(truck.Processed.Zinterp2))./length(truck.Processed.Zinterp2);
meanTruck(perGood<0.5) = NaN;
minTruck = nanmin(truck.Processed.Zinterp2);
minTruck(perGood<0.5) = NaN;

 
droneh = plot(drone.Processed.x,meanDrone);
truckh = plot(truck.Processed.x,meanTruck);
dronehmin = plot(drone.Processed.x,minDrone);
truckhmin = plot(truck.Processed.x,minTruck);

%% 
  plot([-120 20],[1 1]*tideMean,':k')
    ylim([0.5 3])
    xlim([-120 25])
    title([datestr(minTime) ' Drone Hover ' num2str(hoverN,'%d')])
    xlabel('X-shore location (m)')
        ylabel('z NAVD (m)')

    box on

legend([droneh,truckh],'drone','truck')    
print(gcf, '-djpeg', [vizdir ' droneFlightCheck_1209_H' num2str(hoverN,'%d') '.jpeg'],'-r300');

%%

for pnum=1:10
    indHour(pnum) = find(datenum(paros(pnum).t) == datenum(dateshift(datetime(minTime,'ConvertFrom','datenum'),'start','hour')));
end

t = paros(pnum).t(indHour(pnum))+ (0:7166)*seconds(0.5)-seconds(1);
%% plot everything at pnum location
figure
ax(1) = subplot(2,1,1);
pnum = 4;
plot(t,paros(pnum).etaCorrected(:,indHour(pnum))-paros(pnum).offset+paros(pnum).z,'linewidth',2)
hold on

%find xlocation in processed drone/truck data

xloc = knnsearch(drone.Processed.x',paros(pnum).crossshore);

plot(datetime(drone.Processed.t,'ConvertFrom','datenum'),drone.Processed.Zinterp2(:,xloc))
plot(datetime(truck.Processed.t,'ConvertFrom','datenum'),truck.Processed.Zinterp2(:,xloc),'.')
plot(datetime(truck.Processed.t,'ConvertFrom','datenum'),truck.Processed.Zinterp(:,xloc))

xlim(datetime([minTime maxTime],'ConvertFrom','datenum'))
title(['Hover ' num2str(hoverN) ' : x= ' num2str(paros(pnum).crossshore,'%2.2f') 'm'])
legend(paros(pnum).name(1:3),'drone','truck')   
ylabel('zNAVD (m)')


ax(2) = subplot(2,1,2);
P4x = find(truck.xyzti.x>floor(paros(4).UTMx) & truck.xyzti.x<ceil(paros(4).UTMx));
plot(datetime(truck.xyzti.t(P4x),'ConvertFrom','datenum'),truck.xyzti.i(P4x),'.')


xlim(datetime([minTime maxTime],'ConvertFrom','datenum'))
% title(['Hover ' num2str(hoverN) ' : x= ' num2str(paros(pnum).crossshore,'%2.2f') 'm'])
% legend(paros(pnum).name(1:3),'drone','truck')   
ylabel('zNAVD (m)')
linkaxes(ax,'x')
% print(gcf, '-djpeg', [vizdir 'droneTSCheck_0224_H' num2str(hoverN,'%d') '_' paros(pnum).name(1:3) '.jpeg'],'-r300');


