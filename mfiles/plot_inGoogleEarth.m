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

truck = load('../mat/lidar/truck/RiProcess_NoRiPrecision - Scanner VZ-2000 Horizontal - 191209_174131_Scanner_VZ-2000_Horizontal_0 - originalpoints, (last).mat');
k = boundary(truck.xyzti.x(1:100:end),truck.xyzti.y(1:100:end));


[lat,lon] = utm2ll(truck.xyzti.x(k),truck.xyzti.y(k),11);
kmlwritepolygon('truck.kml',lat,lon,'EdgeColor',lines(1),'LineWidth',3)

%%
system([cmd fullfile(pwd,'truck.kml')])

%%
drone = load('../mat/lidar/drone/20191209_00582_TorreyRunup_H3.mat');
k = boundary(drone.xyzti.x(1:100:end),drone.xyzti.y(1:100:end));

cmap = lines(2);
[lat,lon] = utm2ll(drone.xyzti.x(k),drone.xyzti.y(k),11);
kmlwritepolygon('drone.kml',lat,lon,'EdgeColor',cmap(2,:),'LineWidth',3)

%%
system([cmd fullfile(pwd,'drone.kml')])
