%gif of waves
clear
% script to compare the drone with the truck lidar
% will spit out 2 figures: pcolor of coherence between truck and drone, and
% spectra @ a paros location
%choose hover number and paros number and load data
hovern = 1;
parosn = 2;
datadir = '/Volumes/FiedlerBot8000/scarp/';
figureDir = '../viz/';
figureName = 'truckanddrone';
% hoverdate = '20200224';
hoverdate = '20191214';

hoverdatefolder = ['../mat/' hoverdate '/'];
% drone = load([datadir '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
drone = load([datadir '/mat/lidar/drone/20191214_H' num2str(hovern) '_navd88_geoid12b_10cm.mat']);

truck = load([datadir '/mat/lidar/truck/' hoverdate '_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
droneR = load([hoverdatefolder 'Drone_Hover_' num2str(hovern,'%02.0f') '_runupstats_10cm.mat']);
truckR = load([hoverdatefolder 'Truck_Hover_' num2str(hovern,'%02.0f') '_runupstats_10cm.mat']);
load('../mat/paros.mat');
load([hoverdatefolder 'DroneStartStop_' hoverdate '.mat']);
load(['../mat/' hoverdate '_H' num2str(hovern) '_wiggle.mat'])
dewiggle = 1; %set to 0 if no de-wiggling, 1 if you want to de-wiggle
%%
% make sure the truck and drone are sampling from the same time
tchunk = [tstart(hovern) tstop(hovern)];

tind = knnsearch(truck.Processed.t',tchunk');
truck.Processed.t = truck.Processed.t(tind(1):tind(2));
truck.Processed.Zinterp = truck.Processed.Zinterp(tind(1):tind(2),:);
truck.Processed.Zinterp2 = truck.Processed.Zinterp2(tind(1):tind(2),:);

% put things into datetime format
tstart = datetime(tstart,'ConvertFrom','datenum');
tstop = datetime(tstop,'ConvertFrom','datenum');


%% interpolate lidar data onto the same time x grid
Hz_lidar = 10;
tvecHover = tstart(hovern):seconds(1/Hz_lidar):tstop(hovern);
tdrone = datetime(drone.Processed.t(2:end),'ConvertFrom','datenum');
[Tgrid,Xgrid] = meshgrid(tvecHover,drone.Processed.x);


% WIGGLE OR NO WIGGLE???

% get info for standard deviation plot
% dewiggled stuff
TSdrone = drone.Processed.Zinterp2(2:end,:)-wiggle(2:end);
TXdrone_noW = interp2(datenum(tdrone),drone.Processed.x,TSdrone',datenum(Tgrid),Xgrid);
stdTXdrone_noW = nanstd(TXdrone_noW,0,2);

TXdrone2_noW = inpaint_nans(TXdrone_noW,2); %inpaint nans applies del^2 over whole matrix
TXdrone2 = TXdrone2_noW;
stdTXdrone2_noW = nanstd(TXdrone2,0,2);

%with the wiggles
TSdrone = drone.Processed.Zinterp2(2:end,:);
TXdrone_W = interp2(datenum(tdrone),drone.Processed.x,TSdrone',datenum(Tgrid),Xgrid);
stdTXdrone_W = nanstd(TXdrone_W,0,2);

TXdrone2_W = inpaint_nans(TXdrone_W,2); %inpaint nans applies del^2 over whole matrix
stdTXdrone2_W = nanstd(TXdrone2_W,0,2);


if dewiggle == true
    TXdrone = TXdrone_noW;
    TXdrone2 = TXdrone2_noW;
else
    TXdrone = TXdrone_W ;
    TXdrone2 = TXdrone2_W; 
end

ttruck = datetime(truck.Processed.t,'ConvertFrom','datenum');

TStruck = truck.Processed.Zinterp2;
% TStruck(TStruck>nanmean(TStruck)+3*nanstd(TStruck))=NaN;

TStruck(isnan(datenum(ttruck)),:) = [];
ttruck(isnan(datenum(ttruck))) = [];

TXtruck = interp2(datenum(ttruck),truck.Processed.x,TStruck',datenum(Tgrid),Xgrid);

TXtruck2 = inpaint_nans(TXtruck,2);



dt = 1/Hz_lidar;

IGlength = 25; %seconds
IGfilt = 25/dt;

% xinterp = [-30:options.dx:10];

% for i=1:nt
%     Zinterp2(i,:) = interp1(xx,Processed.Zinterp2(i,:),xinterp);
%     Zinterp(i,:) = interp1(xx,Processed.Zinterp(i,:),xinterp);
% end

Md = mean(movmin(TXdrone2,IGfilt,1,'omitnan'),2); % get moving minimum for foreshore on IG freq
Mt = mean(movmin(TXtruck2,IGfilt,1,'omitnan'),2); % get moving minimum for foreshore on IG freq

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

h(1) = scatter(x,TXdrone(:,500),10,cmap(2,:),'filled');
hold on
h(2) = scatter(x,TXtruck(:,500),10,cmap(1,:),'filled');

ylim([0.5 6])
ylabel('z (m, navd88)')
xlabel('x (m)')
box on

hold on
for i=1:6
    plot(paros(i).crossshore*[1 1],[0.5 6],'color',[0.5 0.5 0.5])
    hold on
    ht(i) = text(paros(i).crossshore,4.6,paros(i).name(1:3),...
        'horizontalalignment','center','background','white',...
        'fontsize',12,'fontweight','bold');
end
ht(1).Position(2) = 5.3; ht(2).Position(2) = 4.8;
%

delete(h(1))
delete(h(2))

igif = 1;
for i=500:1200
    h(1) = scatter(x,TXdrone(:,i),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','none');
    hold on
    h(2) = scatter(x,TXtruck(:,i),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','none');
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
%%
gifname = [figureDir figureName '.gif'];
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0);
    end
end