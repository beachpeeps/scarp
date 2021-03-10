%plot_timestackCompare
%choose hover number  and load data
hovern = 3;
datadir = '/Volumes/FiedlerBot8000/scarp/';

hoverdate = '20191214';
% hoverdate = '20200224';

% switch the limits for analysis, depending on hoverdate, to stop at x=0
if str2double(hoverdate) == 20191214
    xlimit = 1000;
    hovers = [1 3 4];
elseif str2double(hoverdate) == 20200224
    xlimit = 500;
    hovers = 2:5;
end

   
tstack  = load([datadir '/mat/timestacks_LARGE/' hoverdate '_' num2str(hovern) '.mat']);

droneR = load(['../mat/' hoverdate '/Drone_Hover_' num2str(hovern,'%02.0f') '_L1_runupstats_10cm.mat']);
truckR = load(['../mat/' hoverdate '/Truck_Hover_' num2str(hovern,'%02.0f') '_L1_runupstats_10cm.mat']);


hoverdate = '20200224';
hovern = 3;

tstack2  = load([datadir '/mat/timestacks_LARGE/' hoverdate '_' num2str(hovern) '.mat']);
droneR2 = load(['../mat/' hoverdate '/Drone_Hover_' num2str(hovern,'%02.0f') '_L1_runupstats_10cm.mat']);
truckR2 = load(['../mat/' hoverdate '/Truck_Hover_' num2str(hovern,'%02.0f') '_L1_runupstats_10cm.mat']);

load('../mat/paros.mat')    
%%
close all
figwidth = 8;
figheight = 15;
nrow = 2;
ncol = 2;
options.margin = 0.06;
options.visible = 'off';
options.units = 'centimeters';
options.axesfontsize = 6;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);

hFig.PaperSize = [figwidth figheight];
% hFig.PaperPosition = [0 0 figwidth figheight];
indplot = 5500:8500;
cmap = lines(2);
truckcolor = cmap(2,:);
dronecolor = cmap(1,:);
dronecolor = 'b';

drone = tstack.TXdrone(:,indplot);
truck = tstack.TXtruck(:,indplot);
t = datenum(tstack.tvecHover(indplot));
x = tstack.Xgrid(:,1);

axes(ax(1))
pcolor(x,t,drone')
shading flat
% caxis([1.2 3.5])


axes(ax(2))
pcolor(x,t,truck')
shading flat

indplot2 = length(t)-3000:length(t);

drone = tstack2.TXdrone(:,indplot2);
truck = tstack2.TXtruck(:,indplot2);
t = datenum(tstack2.tvecHover(indplot2));
x = tstack2.Xgrid(:,1);

axes(ax(3))
pcolor(x,t,drone')
shading flat
% caxis([1.2 3.5])


axes(ax(4))
pcolor(x,t,truck')
shading flat


ax(1).Title.String = 'Hovering';
ax(2).Title.String = 'Terrestrial';

for i=1:2
    axes(ax(i))
hold on
plot(droneR.Tseries.Xrunup,droneR.Tseries.T,'color',dronecolor,'parent',ax(i),'linewidth',1.5)
plot(truckR.Tseries.Xrunup,truckR.Tseries.T,'color',truckcolor,'parent',ax(i),'linewidth',1.5)
% datetick('y','MM:SS')
startTime = datenum(dateshift(datetime(droneR.Tseries.T(indplot(1000)),'ConvertFrom','datenum'),'start','minute'));
ylim([startTime datenum(datetime(startTime,'ConvertFrom','datenum')+minutes(3))])
ax(i).YTick = datenum(startTime):30/60/60/24:ax(i).YLim(2);
for n=2:5
    plot([paros(n).crossshore paros(n).crossshore],ax(1).YLim,'--','color',[0 0 0 0.3]);
    ht(n) = text(paros(n).crossshore,ax(1).YLim(2)+2/(60*60*24),['P' num2str(n,'%1.0f')],'rotation',90,'fontsize',8);
end
% ax(i).CLim = [1 2.5];
ax(i).Title.Position(2) = ax(i).Title.Position(2)+datenum(seconds(15));
end

for i=3:4
    axes(ax(i))
hold on
plot(droneR2.Tseries.Xrunup,droneR2.Tseries.T,'color',dronecolor,'parent',ax(i),'linewidth',1.5)
plot(truckR2.Tseries.Xrunup,truckR2.Tseries.T,'color',truckcolor,'parent',ax(i),'linewidth',1.5)
% datetick('y','MM:SS')
startTime = datenum(dateshift(datetime(droneR2.Tseries.T(indplot2(1000)),'ConvertFrom','datenum'),'start','minute'));
ylim([startTime datenum(datetime(startTime,'ConvertFrom','datenum')+minutes(3))])
xlabel('X (m)')
ax(i).YTick = datenum(startTime):30/60/60/24:ax(i).YLim(2);

for n=2:5
    plot([paros(n).crossshore paros(n).crossshore],ax(3).YLim,'--','color',[0 0 0 0.3]);
end

% ax(i).Title.Position(2) = ax(i).Title.Position(2)+datenum(seconds(15));
end

ax(1).YTickLabel = datestr(ax(1).YTick,'MM:SS');
ax(3).YTickLabel = datestr(ax(3).YTick,'MM:SS');

ax(2).YLabel.String = [];
ax(2).YTickLabel = [];
ax(4).YLabel.String = [];
ax(4).YTickLabel = [];



ax(1).YLabel.String = 'MM:SS';
ax(3).YLabel.String = 'MM:SS';



for i=1:4
ax(i).Position(1) = ax(i).Position(1)+0.032;
% ax(i).Position(4) = 0.65;
% ax(i).XLabel.Position(2) = ax(1).YLim(1) - 15/(60*60*24);
ax(i).XLabel.FontSize = 8;
ax(i).YLabel.FontSize = 8;
ax(i).XLim = [-120 30];
if i>2
ax(i).XTick = [-120:30:30];
elseif i<3
    ax(i).XTickLabel = [];
end
% ax(i).Position(1) = ax(i).Position(1)+0.04;
ax(i).FontSize = 8;
ax(i).CLim = [1 3];

end


axes(ax(4))
hc = colorbar;
hc.Label.String = 'z navd (m)';
hc.Location = 'southoutside';

hc.Position(1) = ax(3).Position(1);
hc.Position(3) = 0.775;
hc.Position(4) = 0.03;
hc.Position(2) = 0.08;

for i=3:4
    ax(i).Position(2) = 0.18;
end
for i=1:2
    ax(i).Position(2) = 0.55;
end

ss = {'a','b','c','d'};
for i=1:4
        text(0.05,0.95,ss{i},'units','normalized','fontsize',12,'fontweight','bold','parent',ax(i))
end

ax(4).Position(4) = ax(3).Position(4);
ax(3).XTickLabel(end) = {' '};


ax(1).YLabel.Position(1) =  -160;
ax(3).YLabel.Position(1) =  -160;


ha1 = annotation(hFig,'arrow',[0.20 0.26],[0.75 0.79]);
ha2 = annotation(hFig,'arrow',[0.6325 0.6925],[0.75 0.79]);
ha1.HeadWidth = 8;
ha2.HeadWidth = 8;
ha1.HeadLength = 8;
ha2.HeadLength = 8;

ha1.Color = [0.3 0.3 0.3];
ha2.Color = [0.3 0.3 0.3];

print(hFig, '-djpeg', '../viz/paperFigures/timestackCompare_LARGE.jpeg','-r300');
