%plot the truck and drone together
drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H1.mat')
truck = load('../mat/lidar/truck/RiProcess_NoRiPrecision - Scanner VZ-2000 Horizontal - 191209_165234_Scanner_VZ-2000_Horizontal_0 - originalpoints, (last).mat');


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
plot(truck.Processed.t,truck.Processed.Zinterp2(:,350))
hold on
plot(drone.Processed.t,drone.Processed.Zinterp2(:,350))
