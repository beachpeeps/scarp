clear
clf
hoverdates = {'20191214','20200224'};
for hoveri=1:2
clearvars -except hover*
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
hover = dir(['../mat/h_' hoverdate '*']);
for i = 1:length(hover)

hData = load(['../mat/' hover(i).name]);
z_bore(i,:) = hData.h_est2;
z_linear(i,:) = hData.h_est;
meanEta(i,:) = hData.meanEta;
z_linear_gradient(i,:) = hData.meanEta-hData.depth_preA';
z_bore_gradient(i,:) = hData.meanEta-hData.depth2_preA';

end
x = hData.x;

% neaten up the jumbos according to xgrid from lidar
z1 = bindata(X1,z1,min(x):0);
z2 = bindata(X2,z2,min(x):0);




%%
ax(hoveri) = subplot(2,1,hoveri)
hj1 = fill([ min(x):0 0:-1:min(x)],[z1; flipud(z2)],[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on

kk = find(x>-20);
hdr1 = fill([x' fliplr(x')],[meanEta(1,:) fliplr(meanEta(3,:))],[0.5 0.5 0.9],'edgecolor',[0.5 0.5 0.9]);


kk = x<-10;
hc = plot(x(kk),z_bore(:,kk),':','LineWidth',1,'color',[0.4660    0.6740    0.1880 0.2]);
hb = plot(x(kk),movmean(mean(z_bore(:,kk)),80),'LineWidth',2,'color',[0.4660    0.6740    0.1880]);

hc = plot(x(kk),z_bore_gradient(:,kk),':','LineWidth',1,'color',[0.4660    0.6740    0.1880 0.2]);
hbg = plot(x(kk),movmean(mean(z_bore_gradient(:,kk)),80),'--','LineWidth',2,'color',[0.4660    0.6740    0.1880]);

%     
hl = plot(x(kk),z_linear(:,kk),':','LineWidth',1,'color',[0.6350    0.0780    0.1840 0.2]);
hl = plot(x(kk),movmean(mean(z_linear(:,kk)),80),'LineWidth',2,'color',[0.6350    0.0780    0.1840 ]);

hlg_dot = plot(x(kk),z_linear_gradient(:,kk),':','LineWidth',1,'color',[0.6350    0.0780    0.1840 0.2]);
hlg = plot(x(kk),movmean(mean(z_linear_gradient(:,kk)),80),'--','LineWidth',2,'color',[0.6350    0.0780    0.1840 ]);
    
    ylabel('elevation NAVD88 (m)')
    if hoveri==2 
        xlabel('cross-shore (m)')
    end
  ylim([min(min(z_bore(:,kk)))-0.5 max(z1)+0.5])
    xlim([min(x) 0])

    box on
    
    if hoveri ==1
    hleg = legend([hb hbg hl hlg],{'crest-tracking, w/ bore', 'gradient w/bore', 'crest-tracking, linear', 'gradient, linear'});
    hleg.Location = 'southeast';
    title('December 2019')
    end
    
     if hoveri ==2
    hleg = legend([hdr1 hj1],{'SWL range', '7-day bathy range'});
    hleg.Location = 'southeast';
    title('February 2020')
     end
    
    ax(hoveri).Title.Position(2) = 2.5;
end
%move the bottom axis around
ax(2).Position(3) = ax(2).Position(3)/2
ax(2).Position(1) = ax(2).Position(1)+ax(2).Position(3);
ax(2).Position(2) = 0.2;




print(gcf,'-dpng',['../viz/lidarBathy_means_dots.png'])
