dataLocation = '/Volumes/FiedlerBot8000/scarp/';

hovern = 1;

% drone = load([dataLocation '/mat/lidar/drone/20191214_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
drone = load([dataLocation '/mat/lidar/drone/20191214_H' num2str(hovern) '_navd88_geoid12b_10cm.mat']);
% drone = load([dataLocation '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
% truck = load([dataLocation '/mat/lidar/truck/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);

truck = load([dataLocation '/mat/lidar/truck/20191214_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
load('../mat/sensors.mat','P','theta')
%%
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

%convert xyzti to lat/lon coords
% [lat,lon]=utm2ll(drone.xyzti.x,drone.xyzti.y,11*ones(size(drone.xyzti.y)));

% x = c.x;
% y = c.y;
% z = c.z;
% intensity = get_intensity(c);
xx = [-150 drone.Processed.x 50];
[XR, YR] = xyRotate(xx,zeros(size(xx)),-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[lat,lon]=utm2ll(XR,YR,11*ones(size(XR)));

%%
%get standard deviation for width
x = drone.Processed.x;
[XR2, YR2] = xyRotate(drone.xyzti.x,drone.xyzti.y,THETA,XO,YO);

stdDrone_alongshore = nan(size(x(1:end-10)));
meanDrone_alongshore = nan(size(x(1:end-10)));

for i = 1:length(x)-10
xind = find(XR2>=x(i) & XR2<=x(i+10));
stdDrone_alongshore(i) = nanstd(YR2(xind));
meanDrone_alongshore(i) = nanmean(YR2(xind));
end

[XR, YR] = xyRotate(x(1:end-10),meanDrone_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDrone,lonDrone]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-10),meanDrone_alongshore+2*stdDrone_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDroneUP,lonDroneUP]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-10),meanDrone_alongshore-2*stdDrone_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDroneDOWN,lonDroneDOWN]=utm2ll(XR,YR,11*ones(size(XR)));

%get standard deviation for width
x = truck.Processed.x;
[XR2, YR2] = xyRotate(truck.xyzti.x,truck.xyzti.y,THETA,XO,YO);

stdTruck_alongshore = nan(size(x(1:end-10)));
meanTruck_alongshore = nan(size(x(1:end-10)));

for i = 1:length(x)-10
xind = find(XR2>=x(i) & XR2<=x(i+10));
stdTruck_alongshore(i) = nanstd(YR2(xind));
meanTruck_alongshore(i) = nanmean(YR2(xind));
end

[XR, YR] = xyRotate(x(1:end-10),meanTruck_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latTruck,lonTruck]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-10),meanTruck_alongshore+2*stdTruck_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latTruckUP,lonTruckUP]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-10),meanTruck_alongshore-2*stdTruck_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latTruckDOWN,lonTruckDOWN]=utm2ll(XR,YR,11*ones(size(XR)));

%
dronePLot = load(['/Volumes/FiedlerBot8000/scarp/mat/lidar/drone/20191214_H' num2str(hovern) '_navd88_geoid12b_10cm_ParLot.mat']);
% dronePLot = load(['../mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm_ParLot.mat']);
%%
x = dronePLot.Processed.x;
[XR2, YR2] = xyRotate(dronePLot.xyzti.x,dronePLot.xyzti.y,THETA,XO,YO);

stdDronePLot_alongshore = zeros(size(x(1:end-25)));
meanDronePLot_alongshore = nan(size(x(1:end-25)));

for i = 1:length(x)-25
xind = find(XR2>=x(i) & XR2<=x(i+10));
stdDronePLot_alongshore(i) = nanstd(YR2(xind));
meanDronePLot_alongshore(i) = nanmean(YR2(xind));
end

[XR, YR] = xyRotate(x(1:end-25),meanDronePLot_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDronePLot,lonDronePLot]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-25),meanDronePLot_alongshore+2*stdDronePLot_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDronePLotUP,lonDronePLotUP]=utm2ll(XR,YR,11*ones(size(XR)));

[XR, YR] = xyRotate(x(1:end-25),meanDronePLot_alongshore-2*stdDronePLot_alongshore,-THETA,0,0);
XR = XR+XO;
YR = YR+YO;
[latDronePLotDOWN,lonDronePLotDOWN]=utm2ll(XR,YR,11*ones(size(XR)));

%remove nans for fill plotting
nanInd = find(isnan(latDronePLot));
latDronePLotUP(nanInd) = [];
lonDronePLotUP(nanInd) = [];
latDronePLotDOWN(nanInd) = [];
lonDronePLotDOWN(nanInd) = [];
latDronePLot(nanInd) = [];
lonDronePLot(nanInd) = [];


%remove nans for fill plotting
nanInd = find(isnan(latDrone));
latDroneUP(nanInd) = [];
lonDroneUP(nanInd) = [];
latDroneDOWN(nanInd) = [];
lonDroneDOWN(nanInd) = [];
latDrone(nanInd) = [];
lonDrone(nanInd) = [];



%%
truckcolor = [148 62 143]./255;
dronecolor = [45 114 178]./255;

truckcolor2 = [255 125 0]./255;
dronecolor2 = [204 194 0]./255;

close all; 

hFig = figure;
hFig.PaperUnits = 'centimeters';
hFig.PaperSize = [8 4];
hFig.PaperPosition = [0 0 8 4];
% hFig.Position = hFig.PaperPosition*100;

ax1 = gca;


hf = fill([lonDroneUP; flipud(lonDroneDOWN)],[latDroneUP; flipud(latDroneDOWN)],dronecolor);
hold on
hf.EdgeColor = 'none'; hf.FaceAlpha = 0.5;
plot(lonDrone,latDrone,'color',dronecolor)

hf2 = fill([lonTruckUP; flipud(lonTruckDOWN)],[latTruckUP; flipud(latTruckDOWN)],truckcolor);
hf2.EdgeColor = truckcolor;

hf3 = fill([lonDronePLotUP; flipud(lonDronePLotDOWN)],[latDronePLotUP; flipud(latDronePLotDOWN)],dronecolor);
hf3.EdgeColor = 'none'; hf3.FaceAlpha = 0.5;
plot(lonDronePLot,latDronePLot,'color',dronecolor)


F = load('../mat/20200224_H2_dronetruckstd.mat');
hf4 = fill([F.lonDroneUP; flipud(F.lonDroneDOWN)],[F.latDroneUP; flipud(F.latDroneDOWN)],dronecolor2);
hf4.EdgeColor = 'none'; hf4.FaceAlpha = 0.5;
plot(F.lonDrone,F.latDrone,'color',dronecolor2)

hf5 = fill([F.lonTruckUP; flipud(F.lonTruckDOWN)],[F.latTruckUP; flipud(F.latTruckDOWN)],truckcolor2);
hf5.EdgeColor = truckcolor2;

hf6 = fill([F.lonDronePLotUP; flipud(F.lonDronePLotDOWN)],[F.latDronePLotUP; flipud(F.latDronePLotDOWN)],dronecolor2);
hf6.EdgeColor = 'none'; hf6.FaceAlpha = 0.5;
plot(F.lonDronePLot,F.latDronePLot,'color',dronecolor2)


% xlim([-117.26125 -117.259])
% ylim([32.927 32.9276])

% text(,Yp(1)+1,'P1','color','r')
% text(Xp(2),Yp(2)+1,'P2','color','r')
% 
[lonVec, latVec, imag] = plot_google_map('Axis',ax1,'MapScale', 1,'MapType','satellite','ShowLabels',0,'Alpha',0.2,'Height',300);
lonVec = lonVec+1e-4;

h = image(lonVec,latVec,imag, 'Parent', ax1);
set(ax1,'YDir','Normal')
set(h,'tag','gmap')
set(h,'AlphaData',0.6)
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)

scatter(P.DeploymentLongitude(1:6),P.DeploymentLatitude(1:6),10,'k','filled')
hold on
scatter(P.DeploymentLongitude(7:10),P.DeploymentLatitude(7:10),5,'k','filled')
plot(lon,lat,':w')
text(P.DeploymentLongitude(3:6),P.DeploymentLatitude(3:6)-5e-5,P.Properties.RowNames(3:6),'FontSize',8)    
    
xlim([-117.26125 -117.259])
ylim([32.927 32.9278])
box on

hl = legend([hf hf2 hf4 hf5],{'Drone: Dec H1','Truck: Dec H1','Drone: Feb H2','Truck: Feb H2'});
hl.FontSize = 8;
hl.NumColumns = 2;
ax1.FontSize = 6;
ax1.Position(3) = 0.85;

ax1.DataAspectRatio = [1 1 1];
ax1.Position(2) = 0.03;
hl.Position(2) = 0.85;
hl.Position(1) = 0.43;


print(gcf, '-djpeg', ['../viz/paperFigures/aerialView.jpeg'],'-r300');
