%plot the truck and drone together
drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H3.mat')
truck = load('../mat/lidar/truck/RiProcess_NoRiPrecision - Scanner VZ-2000 Horizontal - 191209_174131_Scanner_VZ-2000_Horizontal_0 - originalpoints, (last).mat');
figureDir = '../viz/';
figureName = 'truckanddrone';

%% take a look at the length of time of each collect
% figure
% plot(drone.Processed.t,zeros(size(drone.Processed.t)))
% hold on
% plot(truck.Processed.t,ones(size(truck.Processed.t)))
% ylim([-1 2])
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
% %%
% figwidth = 6;
% figheight = 3;
% nrow = 1;
% ncol = 1;
% 
% [hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,'units','inches');
% P3 = -37.2271;
% xloc = max(find(drone.Processed.x<-37.2271));
% 
% plot(tTruck+seconds(18),truck.Processed.Zinterp2(:,xloc),'linewidth',2)
% hold on
% plot(tDrone,drone.Processed.Zinterp2(:,xloc),'linewidth',2)
% ylabel('z NAVD88, m')
% title('Water Level at P3')
% legend('truck','drone')
% 
% ax(1).XLim = [datetime(2019,12,9,17,54,0) datetime(2019,12,9,17,58,0)];
% ax(1).Position(2) = 0.15;
% 
% hFig.Color = 'none';
% hFig.Renderer = 'painters';
% print(hFig, [figureDir figureName '.png'],'-r300');

%% interpolate on to 5 Hz grid

% get 5 Hz in datenum
dt = 1/(secondsinday*5);

tstart = datetime(2019,12,9,17,42,0);

tarray = datenum(tstart):dt:datenum(datetime(2019,12,9,17,46,0));
%%
droneZ = nan(length(tarray),441);

for xi = 1:441
    Z = drone.Processed.Zinterp2(:,xi);
    if length(Z(~isnan(Z)))>100
    droneZ(:,xi) = interp1(drone.Processed.t(~isnan(Z)),Z(~isnan(Z)),tarray,'nearest');
    end
end
%%
truckZ = nan(length(tarray),441);

for xi = 1:441
    Z = truck.Processed.Zinterp2(:,xi);
    if length(Z(~isnan(Z)))>100
    truckZ(:,xi) = interp1(truck.Processed.t(~isnan(Z)),Z(~isnan(Z)),tarray,'nearest');
    end
end



%%
x = drone.Processed.x;
cmap = lines(2);
% i = 500
for i=500:1200
    h(1) = plot(x,droneZ(i,:),'.','color',cmap(2,:),'linewidth',2)
    hold on
    h(2) = plot(x,truckZ(i-18*5,:),'.','color',cmap(1,:),'linewidth',2)
    pause(0.1)
    ylim([0.5 4])
    ylabel('z (m, navd88)')
    xlabel('x (m)')
    if i<1200
    delete(h(1))
    delete(h(2))
    end
    
end