%plot the truck and drone together
drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H3.mat')
truck = load('../mat/lidar/truck/RiProcess_NoRiPrecision - Scanner VZ-2000 Horizontal - 191209_174131_Scanner_VZ-2000_Horizontal_0 - originalpoints, (last).mat');
figureDir = '../viz/';
figureName = 'truckanddrone';

%% take a look at the length of time of each collect
figure
plot(drone.Processed.t,zeros(size(drone.Processed.t)))
hold on
plot(truck.Processed.t,ones(size(truck.Processed.t)))
ylim([-1 2])
%% 
% things are still in datenum here, watch out. Change:
secondsinday = 24*60*60;
t_truck = (truck.Processed.t-min(truck.Processed.t))*secondsinday;
t_drone = (drone.Processed.t-min(drone.Processed.t))*secondsinday;

% get sampling rate of each
dt_truck = nanmean(diff(t_truck));
fs_truck = 1/dt_truck;

dt_drone = nanmean(diff(t_drone));
fs_drone = 1/dt_drone;

% resample to a 5 Hz grid?
% y = resample(truck.Processed.

%%
clf
tTruck = datetime(truck.Processed.t,'ConvertFrom','datenum');
tDrone = datetime(drone.Processed.t,'ConvertFrom','datenum');
%%
figwidth = 6;
figheight = 3;
nrow = 1;
ncol = 1;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,'units','inches');
xloc = 335;
plot(tTruck+seconds(18),truck.Processed.Zinterp2(:,xloc),'linewidth',2)
hold on
plot(tDrone,drone.Processed.Zinterp2(:,xloc),'linewidth',2)
ylabel('z NAVD88, m')
title(['Water Level, x = ' num2str(drone.Processed.x(xloc)) ' m'])
legend('truck','drone')

ax(1).XLim = [datetime(2019,12,9,17,54,0) datetime(2019,12,9,17,58,0)];
ax(1).Position(2) = 0.15;

print(hFig, '-djpeg', [figureDir figureName '.jpeg'],'-r300');

