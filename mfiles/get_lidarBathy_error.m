clear


%

Xarray = -100:0.1:-0.1;

hoverdates = {'20191214','20200224'};
for hoveri = 1:2
clearvars -except hover* ylims ax hFig meanz err* Xarray
    hoverdate = hoverdates{hoveri};
% hoverdate = '20200224';

load('../mat/sensors.mat','P','theta')
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

if str2double(hoverdate) == 20191214
    
    jumbo1 = dlmread('../data/20191211_MOP_TP_DM.txt');
    jumbo2 = dlmread('../data/20191217_MOP_TP_DM.txt');
    
elseif str2double(hoverdate) == 20200224
    
    jumbo1 = dlmread('../data/20200227_MOP582.txt');
    jumbo2 = dlmread('../data/20200221_MOP582.txt');
    
end

yUTM = jumbo1(:,3);
xUTM = jumbo1(:,4);
z = jumbo1(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-100 & X<10 & Y>-15 & Y<15);

X1 = X(tarea);
Y1 = Y(tarea);
z1 = z(tarea);


yUTM = jumbo2(:,3);
xUTM = jumbo2(:,4);
z = jumbo2(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-100 & X<10 & Y>-15 & Y<15);

X2 = X(tarea);
Y2 = Y(tarea);
z2 = z(tarea);

% get lidar bathy data
hover = dir(['../mat/h_' hoverdate '*drone*']);
for i = 1:length(hover)
    hData = load(['../mat/' hover(i).name]);
    z_bore(i,:) = hData.h_est2;
    z_linear(i,:) = hData.h_est;
    meanEta(i,:) = hData.meanEta;
    z_linear_gradient(i,:) = hData.meanEta-hData.depth_preA';
    z_bore_gradient(i,:) = hData.meanEta-hData.depth2_preA';
end

hover = dir(['../mat/h_' hoverdate  '*truck*']);
for i = 1:length(hover)
    hDataTruck = load(['../mat/' hover(i).name]);
    z_boreTruck(i,:) = hDataTruck.h_est2;
    z_linearTruck(i,:) = hDataTruck.h_est;
    meanEtaTruck(i,:) = hDataTruck.meanEta;
    z_linear_gradientTruck(i,:) = hDataTruck.meanEta-hDataTruck.depth_preA';
    z_bore_gradientTruck(i,:) = hDataTruck.meanEta-hDataTruck.depth2_preA';
    
    
end
x = hData.x;

% neaten up the jumbos according to xgrid from lidar
z1 = bindata(X1,z1,min(x):0);
z2 = bindata(X2,z2,min(x):0);

meanzz = mean([z1 z2]');
% meanzz = z1;


meanz(hoveri,:) = interp1(min(x):0,meanzz,Xarray);


if hoveri== 1
    err_drone(hoveri,:) = nanmean(z_bore-meanz(hoveri,:));
    err_truck(hoveri,:) = nanmean(z_boreTruck-meanz(hoveri,:));
end
if hoveri == 2
    err_drone(hoveri,501:end) = nanmean(z_bore-meanz(hoveri,501:end));
    err_truck(hoveri,501:end) = nanmean(z_boreTruck-meanz(hoveri,501:end));
end
end
%%
close all
figwidth = 8;
figheight = 4;
nrow = 1;
ncol = 1;
options.units = 'centimeters';
options.axesfontsize = 9;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);

plot(Xarray,smoothdata(err_drone(1,:),'movmean',50),'r')
hold on
plot(Xarray,smoothdata(err_truck(1,:),'movmean',50),'b')

plot(Xarray(501:end),smoothdata(err_drone(2,501:end),'movmean',50),'-.r')
hold on
plot(Xarray(501:end),smoothdata(err_truck(2,501:end),'movmean',50),'-.b')
plot([-100 0],[0 0],':k')
hl = legend('Dec Drone','Dec Truck','Feb Drone','Feb Truck');
xlabel('cross shore (m)')
ylabel('Bias (m)')
ax.FontSize = 8;
ax.Position(1) = 0.15;
ax.Position(2) = 0.25;
ax.Position(3) = 0.8;
hl.Location = 'South';
hl.NumColumns = 2;
hl.Position(1) = hl.Position(1)-0.05;
hl.Position(2) = hl.Position(2)+0.05;
hl.Box = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'on';




print(hFig, '-dpdf', ['../viz/paperFigures/bathy_error.pdf']);
