clear
close all
figwidth = 11.4;
figheight = 16;
nrow = 4;
ncol = 2;
% options.bmargin = 0.15;
% options.visible = 'off';
options.units = 'centimeters';
options.axesfontsize = 8;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);
ax(1).Position(3) = 0.85;
% ax(2).Visible = 'off';
ax(2).Position = ax(4).Position;
ax(4).Visible = 'off';
ax(3).Position(3) = 0.325;
ax(2).Position(1) = 0.525;
ax(2).Position(3) = 0.425;
ax(5).Position(3) = 0.85;
ax(6).Visible = 'off';
ax(7).Position(3) = 0.85;
ax(8).Visible = 'off';
ax(7).Position(4) = 0.1;
ax(5).Position(4) = 0.1;

ax(5).Position(2) = ax(7).Position(2)+ax(7).Position(4)+0.01;

ax(3).Position(2) = ax(5).Position(2)+ax(5).Position(4)+0.04;
ax(2).Position(2) = ax(5).Position(2)+ax(5).Position(4)+0.04;

for i=1:4
ax(i).Position(4) = 0.27;
end

ax(1).Position(2) = ax(2).Position(2)+ax(2).Position(4)+0.04;

%move others down
% ax(3).Position(2) = 0.5;
% ax(2).Position(2) = ax(3).Position(2);

% ax(2).Position(4) = 0.27;
% ax(3).Position(4) = ax(2).Position(4);
% ax(1).Position(4) = ax(2).Position(4);
% ax(1).Position(2) = ax(3).Position(2)+ax(3).Position(3)+0.01;

for i=1:length(ax)
    text(0.1,0.9,num2str(i),'units','normalized','parent',ax(i))
end

print(hFig,'-dpdf','../viz/paperFigures/lidarBathy_means_withErrorstd.pdf')



%%
%get swash zone area

droneR = load(['../mat/20191214/Drone_Hover_01_L1_runupstats_10cm.mat']);
mn = mean(droneR.Tseries.Xrunup)-2*std(droneR.Tseries.Xrunup); 
mx = mean(droneR.Tseries.Xrunup)+2*std(droneR.Tseries.Xrunup); 
hfillD = fill([mn mn mx mx],[-2 3 3 -2],'b','parent',ax(1)) ;
hold on
droneR = load(['../mat/20200224/Drone_Hover_02_L1_runupstats_10cm.mat']);
mn = mean(droneR.Tseries.Xrunup)-2*std(droneR.Tseries.Xrunup); 
mx = mean(droneR.Tseries.Xrunup)+2*std(droneR.Tseries.Xrunup); 
hfillF = fill([mn mn mx mx],[-2 3 3 -2],'b','parent',ax(2)) ;
hold on
%%
Xarray = -100:0.1:-0.1;

hoverdates = {'20191214','20200224'};
for hoveri=1:2
clearvars -except hover* ylims ax hFig Xarray meanz err_drone err_truck std_truck std_drone hfill*
hoverdate = hoverdates{hoveri};

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

% meanz = mean([z1 z2]');




%%
% ax(hoveri) = subplot(2,1,hoveri)
axes(ax(hoveri))
hold on
hj1 = fill([ min(x):0 0:-1:min(x)],[z1; flipud(z2)],[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
hold on

kk = find(x>-20);
hdr1 = fill([x' fliplr(x')],[meanEta(1,:) fliplr(meanEta(3,:))],[0.5 0.5 0.9],'edgecolor',[0.5 0.5 0.9]);

smoothplot = 50;
kk = x<-10;
% kk(1:40) = 0;
hc = plot(x(kk),z_bore(:,kk),':','LineWidth',1,'color',[0.4660    0.6740    0.1880]);
hold on
hb = plot(x(kk),movmean(mean(z_bore(:,kk)),smoothplot,'omitnan'),'LineWidth',2,'color',[0.4660    0.6740    0.1880]);
% face
% kkT = x>-40 & x<-10
kkT = kk
% hbT = plot(x(kkT),movmean(mean(z_boreTruck(:,kkT)),smoothplot,'omitnan'),'LineWidth',2,'color',[0.4660    1    0.1880]);



% hc = plot(x(kk),z_bore_gradient(:,kk),':','LineWidth',0.5,'color',[0.4660    0.6740    0.1880 0.2]);
% hbg = plot(x(kk),movmean(mean(z_bore_gradient(:,kk)),smoothplot,'omitnan'),'--','LineWidth',2,'color',[0.4660    0.6740    0.1880]);

% hbgT = plot(x(kkT),movmean(mean(z_bore_gradientTruck(:,kkT)),smoothplot,'omitnan'),'--','LineWidth',2,'color',[0.4660    0.840    0.1880]);


%     
hl = plot(x(kk),z_linear(:,kk),':','LineWidth',1,'color',[0.6350    0.0780    0.1840]);
hl = plot(x(kk),movmean(mean(z_linear(:,kk)),smoothplot,'omitnan'),'LineWidth',2,'color',[0.6350    0.0780    0.1840 ]);

% hlT = plot(x(kkT),movmean(mean(z_linearTruck(:,kkT)),smoothplot,'omitnan'),'LineWidth',2,'color',[0.9    0.0780    0.9 ]);


% hlg_dot = plot(x(kk),z_linear_gradient(:,kk),':','LineWidth',0.5,'color',[0.6350    0.0780    0.1840 0.2]);
% hlg = plot(x(kk),movmean(mean(z_linear_gradient(:,kk)),smoothplot,'omitnan'),'--','LineWidth',2,'color',[0.6350    0.0780    0.1840 ]);
% hlgT = plot(x(kkT),movmean(mean(z_linear_gradientTruck(:,kkT)),smoothplot,'omitnan'),'--','LineWidth',2,'color',[0.2350    0.2780    0.4840 ]);

    
    ylabel('elevation NAVD88 (m)')
%     if hoveri==2 
%         xlabel('cross-shore (m)')
%     end     
  ylim([min(min(z_bore(:,kk)))-0.5 max(z1)+0.5])
    xlim([min(x) 0])

    box on
    
    if hoveri ==1
        %% add in the foreshore from the Bulk runupstats file


        load('../mat/20191214/Drone_Hover_04_L1_runupstats_10cm.mat', 'Bulk')
        plot(Bulk.foreshoreX,Bulk.foreshore,'parent',ax(1),'linewidth',2,'color','k')

%     hleg = legend([hb hbg hl hlg],{'crest-tracking, w/ bore', 'gradient w/bore', 'crest-tracking, linear', 'gradient, linear'});
    hleg = legend([hb hl],{'crest-tracking, w/ bore', 'crest-tracking, linear'});

    hleg.Location = 'southeast';
    title('December 2019')
    ylims = ax(1).YLim;
    end
    
     if hoveri ==2
         load('../mat/20200224/Drone_Hover_05_L1_runupstats_10cm.mat', 'Bulk')
         hforeshore = plot(Bulk.foreshoreX,Bulk.foreshore,'parent',ax(2),'linewidth',2,'color','k');

    hleg = legend([hdr1 hj1 hforeshore],{'SWL range', '7-day bathy range','LiDAR foreshore'});
    hleg.Location = 'southeast';
    title('February 2020')
     end
    
    ax(hoveri).Title.Position(2) = 2.25;
    meanzz = mean([z1 z2]');

% meanzz = z1;


meanz(hoveri,:) = interp1(min(x):0,meanzz,Xarray);


if hoveri== 1
    err_drone(hoveri,:) = nanmean(z_bore-meanz(hoveri,:));
    err_truck(hoveri,:) = nanmean(z_boreTruck-meanz(hoveri,:));
    std_drone(hoveri,:) = nanstd(z_bore);
    std_truck(hoveri,:) = nanstd(z_boreTruck);
%     h_insitu(hoveri,:) = hData.meanEta-meanz(hoveri,:);
%     hdiff_norm(hoveri,:,:) = ((hData.meanEta-z_bore)-h_insitu(hoveri,:))./h_insitu(hoveri,:);
end
if hoveri == 2
    err_drone(hoveri,501:1000) = nanmean(z_bore-meanz(hoveri,501:end));
    err_truck(hoveri,501:1000) = nanmean(z_boreTruck-meanz(hoveri,501:end));
    std_drone(hoveri,501:1000) = nanstd(z_bore);
    std_truck(hoveri,501:1000) = nanstd(z_boreTruck);
%     h_insitu(hoveri,501:end) = hData.meanEta-meanz(hoveri,501:end);
%     hdiff_norm(hoveri,:,:) = ((hData.meanEta-z_bore)-h_insitu(hoveri,:))./h_insitu(hoveri,:);
end
    
end
%%
hData = load('../mat/h_20191214_H3drone.mat');
axes(ax(3))
tstep = 10;
ti = 1;
tdiff = hData.wave(1).t(1)-hData.wave(2).t(1);
while ti<=length(hData.wave(1).t)
if mod(ti,tstep) == 0 
            h1 = plot(hData(1).x,hData.eta(1:1000,hData.wave(1).t(ti))+ti/tstep,'color',[0 0 0 0.3]); % wave
            hold on
            tdiff = hData.wave(1).t(ti)-hData.wave(2).t(ti);

%             plot(hData(1).x,hData.eta(1:1000,hData.wave(2).t(ti+tdiff))+ti/tstep,'color','r'); % wave

            kk = x>=hData.wave(1).xpeak(ti) & x<=hData.wave(1).xtrough(ti);
            h3(1) = plot(hData.wave(1).xpeak(ti),hData.wave(1).hpeak(ti)+ti/tstep,'or','markersize',3); %roller face
%             h1(2) = plot(wave.x(xi:xtroughInd(ti)),eta(xi:xtroughInd(ti),ti)+ti/tstep,'c','linewidth',2);
            h3(2) = plot(hData.wave(1).xtrough(ti),hData.wave(1).htrough(ti)+ti/tstep,'sqk','markersize',3); %roller face
            
            if (ti+tdiff)<=numel(hData.wave(2).xpeak)
            plot(hData.wave(2).xpeak(ti+tdiff),hData.wave(2).hpeak(ti++tdiff)+ti/tstep,'or','markersize',3); %roller face
%             h1(2) = plot(wave.x(xi:xtroughInd(ti)),eta(xi:xtroughInd(ti),ti)+ti/tstep,'c','linewidth',2);
            plot(hData.wave(2).xtrough(ti++tdiff),hData.wave(2).htrough(ti++tdiff)+ti/tstep,'sqk','markersize',3); %roller face
            end
%             h3(1) = plot(x(kk),eta(kk,wave.t(ti))+ti/tstep,'k','linewidth',3); % wave
%             xplot = find(~isnan(ggt(:,wave.t(ti))));
%             h3(3) = scatter(x(xplot),eta(xplot,wave.t(ti))+ti/tstep,50,cphase(xplot,wave.t(ti)),'filled');
            
end
        ti = ti+1;
end
%
% caxis([0 5])

xlim([-100 -10])
% hc = colorbar;
ax(3).YAxis.TickLabels =num2cell(ax(3).YAxis.TickValues*tstep/10);
% ylabel(hc,'C_p (m/s)')
% xlabel('cross-shore (m)')

ax(3).YTickMode = 'manual';
ylim([0 18])
ylabel('time (s)')
% ylabel(hc,'C_p (m/s)')
%%

ss = {'a','c','b'};
for i=1:3
    text(0.05,0.95,ss{i},'units','normalized','fontsize',14,'fontweight','bold','parent',ax(i))
end
%move the bottom axis around


% ax(2).Title.Position(2) = 2.25;
% ax(2).Position(3) = ax(2).Position(3)/2
% ax(2).Position(1) = ax(2).Position(1)+ax(2).Position(3);
% ax(2).Position(2) = 0.2;

for i=1:2
    ax(i).YLim = [-1.25 2.75];
end
% ax(2).YLim = [-0.25 2.5];

%%
truckcolor2 = [255 125 0]./255;
axes(ax(5))
kk = Xarray<=-10;
hh = plot(Xarray(kk),smoothdata(err_drone(1,(kk)),'movmean',50),'color',[0.4660    0.6740    0.1880]);
hold on
hs = plot(Xarray(kk),smoothdata(err_truck(1,(kk)),'movmean',50),'color',truckcolor2);



kk2 = Xarray>=-50 & Xarray<=-10;
hh2 = plot(Xarray(kk2),smoothdata(err_drone(2,kk2),'movmean',50),'-.','color',[0.4660    0.6740    0.1880]);
hold on
hs2 = plot(Xarray(kk2),smoothdata(err_truck(2,kk2),'movmean',50),'-.','color',truckcolor2);
plot([-100 0],[0 0],':k')
[hl,icons,plots,txt] = legend([hh hs hh hh2],'hovering','stationary','Dec 2019','Feb 2020');
hl.NumColumns = 2;
hl.Location = 'Southwest';
hl.FontSize = 8;

ylabel('Bias (m)')

axes(ax(7))
kk = Xarray<=-10;
hh = plot(Xarray(kk),smoothdata(std_drone(1,(kk)),'movmean',50),'color',[0.4660    0.6740    0.1880]);
hold on
hs = plot(Xarray(kk),smoothdata(std_truck(1,(kk)),'movmean',50),'color',truckcolor2);



kk2 = Xarray>=-50 & Xarray<=-10;
hh2 = plot(Xarray(kk2),smoothdata(std_drone(2,kk2),'movmean',50),'-.','color',[0.4660    0.6740    0.1880]);
hold on
hs2 = plot(Xarray(kk2),smoothdata(std_truck(2,kk2),'movmean',50),'-.','color',truckcolor2);
plot([-100 0],[0 0],':k')
% [hl,icons,plots,txt] = legend([hh hs hh hh2],'hovering','stationary','Dec 2019','Feb 2020');
% hl.NumColumns = 2;
% hl.Location = 'Southwest';
% hl.FontSize = 8;
xlabel('cross shore (m)')
ylabel('std (m)')

% ax(7) = ax(5);
% axes(ax(7))
% hl2 = legend([hh hh2],'Dec 2019','Feb 2020');
% hl2.NumColumns = 2;
% hl2.Location = 'NorthEast';
% hl2.FontSize = 8;

for i=1:8
    ax(i).FontSize = 8;
end
ax(5).XTickLabel = [];
ax(5).YLim = [-0.4 0.4];
    text(0.05,0.9,'d','units','normalized','fontsize',14,'fontweight','bold','parent',ax(5));
    text(0.05,0.9,'e','units','normalized','fontsize',14,'fontweight','bold','parent',ax(7));
    
ax(7).YLim = [0 0.4];
ax(7).YTick = [0:0.1:0.4];
ax(7).YTickLabel(end) = {' '};
ax(7).YTickMode = 'manual';
hl.FontSize = 8;

hfillF.EdgeColor = 'none';
hfillF.FaceAlpha = 0.1;
uistack(hfillF,'bottom')

hfillD.EdgeColor = 'none';
hfillD.FaceAlpha = 0.1;
uistack(hfillD,'bottom')

% print(gcf,'-dpng',['../viz/lidarBathy_means_dots.png'])
print(hFig,'-dpdf','../viz/paperFigures/lidarBathy_means_withErrorstd.pdf')

