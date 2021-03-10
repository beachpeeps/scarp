clear all
close all
dataLocation = '/Volumes/FiedlerBot8000/scarp/';

% c = lasdata('../data/las/20191208_00568_00636_2250_NoWaves_TorreyDelMar_000227.las');
% criprap = lasdata('../data/las/20191208_00568_00636_2250_NoWaves_TorreyDelMar_000228.las');

% c = lasdata('../data/las/20191214_00568_00636_0152_NoWaves_TorreyDelMar_000227.las');
% criprap = lasdata('../data/las/20191214_00568_00636_0152_NoWaves_TorreyDelMar_000228.las');

c = lasdata('/Volumes/FiedlerBot8000/scarp/data/las/20200224_00568_00590_0320_NoWaves_TorreyIntensive_000227.las');
criprap = lasdata('/Volumes/FiedlerBot8000/scarp/data/las/20200224_00568_00590_0320_NoWaves_TorreyIntensive_000228.las');
load('/Volumes/FiedlerBot8000/scarp/mat/lidar/drone/20200224_00582_TorreyRunup_H2_10cm.mat','Processed')
lidarFeb = Processed;
load('/Volumes/FiedlerBot8000/scarp/mat/lidar/drone/20191214_H1_navd88_geoid12b_10cm.mat','Processed')
lidarDec = Processed;


minTime = datenum(dateshift(datetime(nanmin(lidarDec.t),'ConvertFrom','datenum'),'start','minute'));
maxTime = datenum(dateshift(datetime(nanmax(lidarDec.t),'ConvertFrom','datenum'),'end','minute'));

[t_verified,t_predicted,verified,predicted,tideInfo] = getNOAAtide(minTime,maxTime,'9410230');

indTide = find(t_verified>minTime & t_verified<=maxTime);

tideDec = mean(verified(indTide));

minTime = datenum(dateshift(datetime(nanmin(lidarFeb.t),'ConvertFrom','datenum'),'start','minute'));
maxTime = datenum(dateshift(datetime(nanmax(lidarFeb.t),'ConvertFrom','datenum'),'end','minute'));

[t_verified,t_predicted,verified,predicted,tideInfo] = getNOAAtide(minTime,maxTime,'9410230');

indTide = find(t_verified>minTime & t_verified<=maxTime);

tideFeb = mean(verified(indTide));
% for hovern = 2:4
%     hovern = 3

% drone = load([dataLocation '/mat/lidar/drone/20191214_H' num2str(hovern) '_navd88_geoid12b_10cm.mat']);
% truck = load([dataLocation '/mat/lidar/truck/20191214_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
% drone = load([dataLocation '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
% truck = load([dataLocation '/mat/lidar/truck/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
%%
x = c.x;
y = c.y;
z = c.z;
intensity = get_intensity(c);

load('../mat/sensors.mat','P','theta')

THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(x,y,THETA,XO,YO);

tarea = find(XR<10 & YR>-10 & YR<10);

XR = XR(tarea);
YR = YR(tarea);
z = z(tarea);

[xq,yq] = meshgrid(-10:0.1:10);
z1 = griddata(XR,YR,z,xq,yq,'linear');


x = c.x;
y = c.y;
z = c.z;

[XR YR] = xyRotate(x,y,THETA,XO,YO);

tarea = find(XR<10 & YR>-10 & YR<10);

XR = XR(tarea);
YR = YR(tarea);
z = z(tarea);

z2 = griddata(XR,YR,z,xq,yq,'linear');

x = criprap.x;
y = criprap.y;
z = criprap.z;

[XR YR] = xyRotate(x,y,THETA,XO,YO);

tarea = find(XR>-10 & XR<10 & YR>-10 & YR<10);

XR = XR(tarea);
YR = YR(tarea);
z = z(tarea);

z1riprap = griddata(XR,YR,z,xq,yq,'linear');

%%

jumbo = dlmread('../data/20200227_MOP582.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-150 & X<10 & Y>-30 & Y<10);

X0227 = X(tarea);
Y0227 = Y(tarea);
z0227 = z(tarea);


%%
jumbo = dlmread('../data/20200221_MOP582.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-150 & X<10 & Y>-30 & Y<10);

X0221 = X(tarea);
Y0221 = Y(tarea);
z0221 = z(tarea);

%%
jumbo = dlmread('../data/20191217_MOP_TP_DM.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-150 & X<10 & Y>-30 & Y<10);

X1217 = X(tarea);
Y1217 = Y(tarea);
z1217 = z(tarea);
%%
jumbo = dlmread('../data/20191211_MOP_TP_DM.txt');
yUTM = jumbo(:,3);
xUTM = jumbo(:,4);
z = jumbo(:,5);
[X, Y] = xyRotate(xUTM,yUTM,THETA,XO,YO);

tarea = find(X>-150 & X<10 & Y>-30 & Y<10);

X1211 = X(tarea);
Y1211 = Y(tarea);
z1211 = z(tarea);

%% grid up December jumbos
XX = -120:0.5:5;
Xg1 = discretize(X1217,XX);
Xg2 = discretize(X1211,XX);

mDec17 = nan(size(XX));
mDec11 = nan(size(XX));
for i=1:length(XX)
    n1 = find(Xg1 == i);
    n2 = find(Xg2 == i);
    if ~isempty(n1)
        mDec17(i) = nanmean(z1217(n1));
    end
    
    if ~isempty(n2)
        mDec11(i) = nanmean(z1211(n2));
    end
end

% get rid of Nans for plotting
mDec11 = interp1(XX(~isnan(mDec11)),mDec11(~isnan(mDec11)),XX);
mDec17 = interp1(XX(~isnan(mDec17)),mDec17(~isnan(mDec17)),XX);

mD11 = smoothdata(mDec11,'movmean',5,'omitnan');
mD17 = smoothdata(mDec17,'movmean',5,'omitnan');

m = [mD11;mD17];
XXm = XX(~isnan(m(1,:)));
m = m(:,~isnan(m(1,:)));

%%
close all
figwidth = 8;
figheight = 4;
nrow = 1;
ncol = 1;
options.margin = 0.05;
options.units = 'centimeters';
options.axesfontsize = 8;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);
cmap = lines(2);
colorD = cmap(1,:);
colorF = cmap(2,:);

hJumbo(1) = plot(XX,mean([mD17;mD11]),'color',colorD);
hold on
hJumboFill(1) = fill([XXm fliplr(XXm)],[max(m) fliplr(min(m))],colorD);
hJumboFill(1).FaceAlpha = 0.15; hJumboFill(1).EdgeColor = 'none';
hold on

% grid up January jumbos

Xg1 = discretize(X0221,XX);
Xg2 = discretize(X0227,XX);

m1 = nan(size(XX));
m2 = nan(size(XX));
for i=1:length(XX)
    n1 = find(Xg1 == i);
    n2 = find(Xg2 == i);
    if ~isempty(n1)
        m1(i) = nanmean(z0221(n1));
    end
    
    if ~isempty(n2)
        m2(i) = nanmean(z0227(n2));
    end
end

% get rid of Nans for plotting
m1 = interp1(XX(~isnan(m1)),m1(~isnan(m1)),XX);
m2 = interp1(XX(~isnan(m2)),m2(~isnan(m2)),XX);

mJ1 = smoothdata(m1,'movmean',5,'omitnan');
mJ2 = smoothdata(m2,'movmean',5,'omitnan');

m = [mJ1;mJ2];
XXm = XX(~isnan(m(1,:)));
m = m(:,~isnan(m(1,:)));


hJumbo(2) = plot(XX,mean([mJ1;mJ2]),'color',colorF);
hold on
hJumboFill(2) = fill([XXm fliplr(XXm)],[max(m) fliplr(min(m))],colorF);
hJumboFill(2).FaceAlpha = 0.15; hJumboFill(2).EdgeColor = 'none';
hold on


%
[Xp, Yp] = xyRotate(P.UTMEastings_Zone11_,P.UTMNorthings_Zone11_,THETA,XO,YO);


% ax = subplot(3,1,1:2);
% plot(xq(:,:),z1(:,:),'k')
% hold on
plot(xq(101,:),z1riprap(101,:),'k')
% hs = scatter(X0227,z0227,'+k');
% hs.MarkerEdgeAlpha = 0.2;
% 
% hs = scatter(X0221,z0221,'+m');
% hs.MarkerEdgeAlpha = 0.2;
% 
% 
% hs = scatter(X1211,z1211,'+g');
% hs.MarkerEdgeAlpha = 0.2;
% 
% hs = scatter(X1217,z1217,'+c');
% hs.MarkerEdgeAlpha = 0.2;

scatter(Xp(1:6),P.DeploymentNAVD88_m_Geoid12B_(1:6),10,'k','filled')
text(Xp(1:6)+0.2,P.DeploymentNAVD88_m_Geoid12B_(1:6)-0.4,P.Properties.RowNames(1:6),'color','k','fontsize',6)
% text(Xp(2)+0.2,P.DeploymentNAVD88_m_Geoid12B_(2)-0.1,'P2','color','k')


% Mean Seal Level in NAVD88
NAVD88mllw = 0.058;
MSLmllw = 0.832;
MSLnavd88 = MSLmllw-NAVD88mllw;
MHWmllw = 1.402;
MHWnavd88 = MHWmllw-NAVD88mllw;


hlid(1) = plot(lidarDec.x,nanmean(lidarDec.Zinterp2),'-.','linewidth',1,'color',colorD);
hlid(2) = plot(lidarFeb.x,nanmean(lidarFeb.Zinterp2),'-.','linewidth',1,'color',colorF);



xlim([-120 10])
ylim([-2.5 5])

plot([ax.XLim(1) -27],[MSLnavd88 MSLnavd88],':k')
plot([ax.XLim(1) -10],tideFeb*[1 1],':','color',colorF)
plot([ax.XLim(1) -6],tideDec*[1 1],':','color',colorD)

% plot(ax.XLim,[MHWnavd88 MHWnavd88],':k')

text(-115,MSLnavd88,['MSL = ' num2str(MSLnavd88,'%2.2f') ' m'],'fontsize',6,'verticalalignment','middle','background','w')
text(-115,tideFeb-0.15,['Feb 24 '],'fontsize',6,'verticalalignment','middle','color',colorF)
text(-115,tideDec+0.15,['Dec 14 '],'fontsize',6,'verticalalignment','middle','color',colorD)

% text(-110,MHWnavd88,['MHW = ' num2str(MHWnavd88,'%2.2f') ' m'],'fontsize',12,'verticalalignment','middle','background','w')


% ylim([0 6])
% title('Average minimum foreshore')

% hl = legend(['Drone: Hover ' num2str(hovern)],['Truck: Hover ' num2str(hovern)],'mobile truck, 2/24 3am LST');
% hl = legend(['Drone: Hover ' num2str(hovern)],['Truck: Hover ' num2str(hovern)],'mobile truck, 12/14 2am LST');
hleg = legend([hlid(1),hlid(2)],'Dec Drone Lidar','Feb Drone Lidar')
hleg.Location = 'northwest'; hleg.Box = 'off';

ax.XLabel.String = 'Cross-shore (m)';
ax.YLabel.String = 'z navd88 (m)';
% ax(1).YLabel.String = 'alongshore (m)';
hc.Label.String = 'z navd88 (m)';


% 
% ax = subplot(3,1,3)
% % plot(xq(101,:),YR(101,:),'k')
% hold on
% % plot(xq(101,:),z1riprap(101,:),'k')
% hs = scatter(X0227,Y0227,'+k');
% hs.MarkerEdgeAlpha = 0.2;
% 
% hs = scatter(X0221,Y0221,'+m');
% hs.MarkerEdgeAlpha = 0.2;
% 
% 
% hs = scatter(X1211,Y1211,'+g');
% hs.MarkerEdgeAlpha = 0.2;
% 
% hs = scatter(X1217,Y1217,'+c');
% hs.MarkerEdgeAlpha = 0.2;
% 
% scatter(Xp(1:6),Yp(1:6),100,colorF,'filled')
% text(Xp(1)+0.2,Yp(1)-0.1,'P1','color',colorF)
% text(Xp(2)+0.2,Yp(2)-0.1,'P2','color',colorF)


% xlim([-10 10])
% ylim([0 6])
% title('Average minimum foreshore')

% hl = legend(['Drone: Hover ' num2str(hovern)],['Truck: Hover ' num2str(hovern)],'mobile truck, 2/24 3am LST');
% hl = legend(['Drone: Hover ' num2str(hovern)],['Truck: Hover ' num2str(hovern)],'mobile truck, 12/14 2am LST');

% hl.Location = 'northwest';
% 
% ax.XLabel.String = 'Cross-shore (m)';
% ax.YLabel.String = 'along-shore (m)';
% ax(1).YLabel.String = 'alongshore (m)';
ax.FontSize = 8;
ax.Position(2) = 0.25;
print(gcf, '-dpdf', ['../viz/paperFigures/transect.pdf'],'-r300','-painters');
% close all
% print(gcf, '-djpeg', ['../viz/views_20200224_H' num2str(hovern) '.jpeg'],'-r300');

% end