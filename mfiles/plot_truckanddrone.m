%plot the truck and drone together
drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H3_05m.mat')
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
xD = drone.Processed.x;
lenX = length(xD);
droneZ = nan(length(tarray),lenX);

for xi = 1:length(drone.Processed.x)
    Z = drone.Processed.Zinterp2(:,xi);
    if length(Z(~isnan(Z)))>100
    droneZ(:,xi) = interp1(drone.Processed.t(~isnan(Z)),Z(~isnan(Z)),tarray,'next');
    end
end
%%
xT = truck.Processed.x;
lenX = length(xT);
truckZ = nan(length(tarray),lenX);

for xi = 1:lenX
    Z = truck.Processed.Zinterp2(:,xi);
    if length(Z(~isnan(Z)))>100
    truckZ(:,xi) = interp1(truck.Processed.t(~isnan(Z)),Z(~isnan(Z)),tarray,'next');
    end
end



%%
x = drone.Processed.x;
cmap = lines(2);

hFig = figure('visible','on');
set(0,'defaultaxesfontsize',12)
set(0,'defaultaxeslinewidth',2)
hFig.PaperUnits = 'centimeters';
hFig.PaperSize = [12 5];
hFig.PaperPosition = [0 0 12 5];
hFig.Position = hFig.PaperPosition*50;
hFig.Position(1:2) = [-1200 250];
hFig.PaperPositionMode = 'manual';
hFig.Color = 'w';

h(1) = scatter(xD,droneZ(500,:),10,cmap(2,:),'filled');
hold on
h(2) = scatter(xT,truckZ(500-18*5,:),10,cmap(1,:),'filled');

ylim([0.5 6])
ylabel('z (m, navd88)')
xlabel('x (m)')
box on
load('../mat/sensors.mat')


%
X = P.UTMEastings_Zone11_;
Y = P.UTMNorthings_Zone11_;
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(X,Y,THETA,XO,YO);
hold on
for i=1:6
    plot([XR(i) XR(i)],[0.5 6],'color',[0.5 0.5 0.5])
    hold on
    ht(i) = text(XR(i),4.6,P.Properties.RowNames(i),...
        'horizontalalignment','center','background','white',...
        'fontsize',12,'fontweight','bold');
end
ht(1).Position(2) = 5.3; ht(2).Position(2) = 4.8;
%

delete(h(1))
delete(h(2))

igif = 1;
for i=500:1200
    h(1) = scatter(xD,droneZ(i,:),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none');
    hold on
    h(2) = scatter(xT,truckZ(i-18*5,:),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none');
    pause(0.1)

    
    drawnow
    F = getframe(hFig);
    [im{igif},map] = frame2im(F);    
    igif = igif+1;
    
    if i<1200
    delete(h(1))
    delete(h(2))
    end
    
end

gifname = [figureDir figureName '.gif'];
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0);
    end
end