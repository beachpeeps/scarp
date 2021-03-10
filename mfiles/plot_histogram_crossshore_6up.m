%histogram of the runup
addpath /Users/juliafiedler/Documents/MATLAB/Inpaint_nans/Inpaint_nans/
addpath /Users/juliafiedler/Documents/SWASH_runup/mfiles/functions/
% script to compare the drone with the truck lidar
%choose hover number and load data

hovern = 2;
datadir = '/Volumes/FiedlerBot8000/scarp/';
hoverdate = '20200224';

hoverdatefolder = ['../mat/' hoverdate '/'];
tstack  = load([datadir '/mat/timestacks/' hoverdate '_' num2str(hovern) '.mat']);

drone = load([datadir '/mat/lidar/drone/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);
truck = load([datadir '/mat/lidar/truck/20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm.mat']);


droneR = load([hoverdatefolder 'Drone_Hover_' num2str(hovern,'%02.0f') '_runupstats_10cm.mat']);
truckR = load([hoverdatefolder 'Truck_Hover_' num2str(hovern,'%02.0f') '_runupstats_10cm.mat']);
load([hoverdatefolder 'DroneStartStop_' hoverdate '.mat']);
load(['../mat/' hoverdate '_H' num2str(hovern) '_wiggle.mat'])
%% interpolate lidar data onto the same time x grid

% get info for standard deviation plot
% dewiggled stuff
stdTXdrone_noW = nanstd(tstack.TXdrone,0,2);
stdTXdrone2_noW = nanstd(tstack.TXdrone2,0,2);

%with the wiggles
k = ~isnan(drone.Processed.t);
TSdrone = drone.Processed.Zinterp2(k,:);
TXdrone_W = interp2(drone.Processed.t(k),drone.Processed.x,TSdrone',datenum(tstack.Tgrid),tstack.Xgrid);
stdTXdrone_W = nanstd(TXdrone_W,0,2);

TXdrone2_W = inpaint_nans(TXdrone_W,2); %inpaint nans applies del^2 over whole matrix
stdTXdrone2_W = nanstd(TXdrone2_W,0,2);



TStruck = truck.Processed.Zinterp2;
ttruck = truck.Processed.t;

% if str2double(hoverdate) == 20200224
% % invalid = find(sum(isnan(TStruck))>=525);
invalid = find(sum(isnan(TStruck),2)>=450);
% 
TStruck(invalid,:) = [];
ttruck(invalid) = [];
% end
TXtruck = interp2(ttruck,truck.Processed.x,TStruck',datenum(tstack.Tgrid),tstack.Xgrid);
TXtruck2 = inpaint_nans(TXtruck,2); %inpaint nans applies del^2 over whole matrix


dt = 1/tstack.Hz_lidar;
IGlength = 25; %seconds
IGfilt = 25/dt;

Md = mean(movmin(tstack.TXdrone2,IGfilt,1,'omitnan')); % get moving minimum for foreshore on IG freq
Mt = mean(movmin(tstack.TXtruck2,IGfilt,1,'omitnan')); % get moving minimum for foreshore on IG freq

%%
addpath ~/Documents/MATLAB/cmocean_v1.4/cmocean/

close all 
clear ax hFig
figwidth = 14;
figheight = 12;
ncol = 2;
nrow = 3;
[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,'Units','centimeters','visible','off','PaperUnits','centimeters');


hFig.PaperUnits = 'centimeters';
hFig.PaperSize = [figwidth figheight];

% hFig.Position = [3000 445 20*figwidth 20*figheight]
% hFig.PaperPositionMode = 'manual';

%

%
clear edges N

for colN = 1:2
    %each column data to plot. (left) un-interpolated (right) de-wiggled
    %and interpolated
    if colN ==1
        truckLidar = TXtruck;
        droneLidar = TXdrone_W;
    elseif colN==2
        truckLidar = TXtruck2;
        droneLidar = tstack.TXdrone2;
    end
    
    

axes(ax(colN))
xrunup = -15:0.1:1;
edgeS = [1.2:0.001:2.8];
clear N* R2e*
for i=1:length(xrunup)
xind = knnsearch(tstack.Xgrid(:,1),xrunup(i));

zd = truckLidar(xind,:)';
% upperlimit = 3*nanstd(zd) + nanmean(zd);
% zd(zd>upperlimit) = NaN;
% zd = zd(~isnan(zd));

[N(i,:),edges] = histcounts(zd,edgeS);
N(i,:) = 100*N(i,:)./sum(N(i,:));
[f,x] = ecdf(zd);
R2e(i) = x(find(f>0.98,1));
meanSS(i) = mean(zd);
minSS(i) =x(find(f>0.001,1));
end
binCenter = edges(1:end-1)+diff(edges);
%
pcolor(xrunup,binCenter,N'); shading flat

hold on



axes(ax(colN+2))
for i=1:length(xrunup)
xind = knnsearch(tstack.Xgrid(:,1),xrunup(i));

zd = droneLidar(xind,:)';
% upperlimit = 3*nanstd(zd) + nanmean(zd);
% zd(zd>upperlimit) = NaN;
% zd = zd(~isnan(zd));

[Nd(i,:),edges] = histcounts(zd,edgeS);
Nd(i,:) = 100*Nd(i,:)./sum(Nd(i,:));
[f,x] = ecdf(zd);
R2ed(i) = x(find(f>0.98,1));
meanSSd(i) = mean(zd);
minSSd(i) = x(find(f>0.001,1));
end
binCenter = edges(1:end-1)+diff(edges);
%
pcolor(xrunup,binCenter,Nd'); shading flat

hold on

axes(ax(colN+4))

pcolor(xrunup,binCenter,Nd'-N'); shading flat
hold on

end
orange = [0 0.8 1]
ss = {'a','b','c','d','e','f'};
for i=1:6
plot(xrunup,R2ed,'color',orange,'linewidth',1,'parent',ax(i))
plot(xrunup,R2e,':','color',orange,'linewidth',1,'parent',ax(i))
% plot(xrunup,meanSSd,':r','linewidth',1,'parent',ax(i))
set(ax(i),'layer','top')
text(0.1,0.9,ss{i},'units','normalized','parent',ax(i),'horizontalalignment','right','fontsize',16,'fontweight','bold')
ax(i).FontSize = 8;
end


for i=1:4
    ax(i).CLim = [0 0.5];
end


cmap1 = cmocean('thermal'); 

% cmap1 = colormap('parula'); 
% cmap1 = flipud(cmap1);


cmap1(1,:) = [1 1 1];

for i=1:6
    ax(i).XMinorTick = 'on';
    if i<5
ax(i).XTickLabel = [];
ax(i).Colormap = cmap1;
    end
end







ss = {'TERRESTRIAL','HOVER','HOVER-TERRESTRIAL'};
leftcol = [1 3 5];
rightcol = leftcol+1;
for i=1:3
    text(0.95,0.1,ss{i},'units','normalized','parent',ax(leftcol(i)),'horizontalalignment','right','fontsize',8)
    text(0.95,0.1,ss{i},'units','normalized','parent',ax(rightcol(i)),'horizontalalignment','right','fontsize',8)
    ax(leftcol(i)).YLabel.String = ' z (m)';
    ax(rightcol(i)).YLabel = [];
    ax(rightcol(i)).YTickLabel = [];

end

for i = leftcol
    ax(i).Position(3) = 0.39;
    ax(i).Position(1) = 0.07;
end

for i= rightcol
    ax(i).Position(3) = 0.43;
end

for i=[2 4 6]
    axes(ax(i))
    cb = colorbar;
    hyy = ylabel(cb,'%')
    ax(i).Position(1) = ax(i).Position(1)-0.06;
end
hyy.Position(1) = 2.06;
%
cmap2 = cmocean('balance','pivot',0); 

for i=5:6
% ax(i).Colormap = cmap2;
ax(i).XLabel.String = 'Cross shore x (m)';
end

ax(5).CLim = [-0.5 0.5];
ax(6).CLim = ax(5).CLim;

ax(5).Colormap = cmocean('balance','pivot',0); 
ax(6).Colormap = ax(5).Colormap;


print(hFig, '-djpeg', ['../viz/paperFigures/histCompare_10cm_H' num2str(hovern,'%1.0f') '_' hoverdate '_6up_v2.jpeg'],'-r300');
% close 

% %%
% %calculate standard deviation
% figure
% plot(drone.Processed.x,stdTXdrone_W)
% hold on
% plot(drone.Processed.x,stdTXdrone2_W)
% plot(drone.Processed.x,stdTXdrone_noW)
% plot(drone.Processed.x,stdTXdrone2_noW)
% plot(drone.Processed.x,nanstd(TXtruck,0,2))
% plot(drone.Processed.x,nanstd(TXtruck2,0,2))
% 
% legend('D wiggle no interp','D no wiggle no interp','D no wiggle interp','T no interp','T interp')
% 
% indX = find(drone.Processed.x>=0 & drone.Processed.x<=1); % find x locations on foreshore area (-2m to 2m)
% 
% stdTXdW = mean(stdTXdrone_W(indX));
% stdTXdnoW = mean(stdTXdrone_noW(indX));
% stdTXt = mean(nanstd(TXtruck2(indX,:),0,2));
