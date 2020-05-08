clear
load('../mat/sensors.mat')

filename = 'sensors.kml';
iconDir = '/Users/juliafiedler/Documents/Illustrations/';

iconFilename = fullfile(iconDir,'psensor.png');
%%
kmlwritepoint(filename,P.DeploymentLatitude,P.DeploymentLongitude,...
    P.DeploymentNAVD88_m_Geoid12B_,...
    'Description',P.Properties.RowNames,'name',P.Properties.RowNames,...
    'Icon',iconFilename,'IconScale',0.4,'altitudemode','relativeToSeaLevel');

%%
cmd = 'open -a Google\ Earth\ Pro ';
fullfilename = fullfile(pwd,filename);
system([cmd fullfilename])

%%

% truck = load('../mat/lidar/truck/RiProcess_NoRiPrecision - Scanner VZ-2000 Horizontal - 191209_174131_Scanner_VZ-2000_Horizontal_0 - originalpoints, (last).mat');
truckdir = '../mat/lidar/truck/';
truckfile = dir([truckdir 'WaveScanDrone*']);
truck = load([truckdir truckfile(4).name]);
% WaveScanDrone_POINTCLOUDS_20200224_185757_Scanner_VZ-2000_Horizontal.mat');
% truck = load('../mat/lidar/truck/WaveScanDrone_POINTCLOUDS_20200224_174302_Scanner_VZ-2000_Horizontal.mat');
% k = boundary(truck.xyzti.x(1:100:end),truck.xyzti.y(1:100:end),0.8);
cmap = lines(3);
[lat,lon] = utm2ll(truck.xyzti.x(1:1000:end),truck.xyzti.y(1:1000:end),11);
lats = [max(lat) min(lat)];
lons = [max(lon) min(lon) ];
kmlwritepolygon('truck.kml',lats,lons,'EdgeColor',cmap(3,:),'LineWidth',3)

%%
system([cmd fullfile(pwd,'truck.kml')])

%%
% drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H3.mat');
drone = load('../mat/lidar/drone/20200224_00582_TorreyRunup_H5.mat');

cmap = lines(3);
% k = boundary(drone.xyzti.x(1:100:end),drone.xyzti.y(1:100:end));
[lat,lon] = utm2ll(drone.xyzti.x(1:1000:end),drone.xyzti.y(1:1000:end),11);
lats = [max(lat) min(lat)];
lons = [max(lon) min(lon) ];
kmlwritepolygon('drone.kml',lats,lons,'EdgeColor',cmap(2,:),'LineWidth',3)

%%
system([cmd fullfile(pwd,'drone.kml')])
system([cmd fullfile(pwd,'sensors.kml')])
