load('../mat/lidar_bathy.mat');

load('../mat/sensors.mat','P','theta')
load('/Volumes/FiedlerBot8000/scarp/mat/lidar/drone/20200224_00582_TorreyRunup_H2_10cm.mat','Processed')
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);
%%

% 1) Latitude (decimal degrees)
% 2) Longitude (decimal degrees)
% 3) Northings (m UTM 11N)
% 4) Eastings (m UTM 11N)
% 5) Elevation (m NAVD88 Geoid12B)
% 6) Time (YYYYMMDDHH.HH)
% 7) Substrate (0= unknown/jetski, 1=sand, 2=cobble, 3=bedrock, 7=unknown/dolly)
jumbo = dlmread('../data/20200227_MOP582.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-100 & X<10 & Y>-10 & Y<10);

X0227 = X(tarea);
Y0227 = Y(tarea);
z0227 = z(tarea);


jumbo = dlmread('../data/20200221_MOP582.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-100 & X<10 & Y>-10 & Y<10);

X0221 = X(tarea);
Y0221 = Y(tarea);
z0221 = z(tarea);

%%
d_lidarinterp = interp1(x_lidar,d_lidar,Processed.x);
h_lidar = nanmean(Processed.Zinterp2)+d_lidarinterp;
%%
clf
scatter(X0221,z0221,'.')
hold on
scatter(X0227,z0227,'.')

plot(Processed.x,nanmean(Processed.Zinterp2),'linewidth',2)
plot(Processed.x,h_lidar,'linewidth',2)

ylabel('z navd88 (m)')
xlabel('cross-shore (m)')

hl = legend('jumbo 02/21','jumbo 02/27','Hover 2: 02/24 nanmean of Zinterp2','Hover 2: derived bathy');
hl.Location = 'southeast';
ylim([-1.5 3.5])

print(gcf, '-djpeg', ['../viz/lidarBathy.jpeg'],'-r300');
