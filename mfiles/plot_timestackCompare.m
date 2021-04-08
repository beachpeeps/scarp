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
load('../mat/paros.mat')    
%%
close all
figwidth = 8;
figheight = 15;
nrow = 1;
ncol = 2;
options.margin = 0.06;
options.visible = 'off';
options.units = 'centimeters';
options.axesfontsize = 6;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);

hFig.PaperSize = [figwidth figheight];
% hFig.PaperPosition = [0 0 figwidth figheight];
indplot = 4050:8050;


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


ax(1).Title.String = 'Drone';
ax(2).Title.String = 'Truck';

for i=1:2
    axes(ax(i))
hold on
plot(droneR.Tseries.Xrunup,droneR.Tseries.T,'color','b','parent',ax(i),'linewidth',1.5)
plot(truckR.Tseries.Xrunup,truckR.Tseries.T,'color','r','parent',ax(i),'linewidth',1.5)
datetick('y','MM:SS')
startTime = datenum(dateshift(datetime(droneR.Tseries.T(indplot(1000)),'ConvertFrom','datenum'),'start','minute'));
ylim([startTime datenum(datetime(startTime,'ConvertFrom','datenum')+minutes(4))])
xlabel('X (m)')

for n=2:5
    plot([paros(n).crossshore paros(n).crossshore],ax(1).YLim,':w');
    ht(n) = text(paros(n).crossshore,ax(1).YLim(2)+2/(60*60*24),['P' num2str(n,'%1.0f')],'rotation',90,'fontsize',8);
end
ax(i).CLim = [1 2.5];
ax(i).Title.Position(2) = ax(i).Title.Position(2)+datenum(seconds(15));
end

ax(2).YLabel.String = [];
ax(2).YTickLabel = [];

axes(ax(1))


ax(1).YLabel.String = 'MM:SS';

axes(ax(2))
hc = colorbar;
hc.Label.String = 'z navd (m)';
hc.Location = 'southoutside';


for i=1:2
ax(i).Position(2) = 0.22;
ax(i).Position(4) = 0.65;
ax(i).XLabel.Position(2) = ax(1).YLim(1) - 15/(60*60*24);
ax(i).XLabel.FontSize = 8;
ax(i).YLabel.FontSize = 8;
ax(i).XTick = [-100:25:0];
ax(i).Position(1) = ax(i).Position(1)+0.04;
ax(i).FontSize = 8;
end

hc.Position(1) = ax(1).Position(1);
hc.Position(3) = 0.775;
hc.Position(4) = 0.03;
hc.Position(2) = 0.1;
ax(2).Position(3) = ax(1).Position(3);

print(hFig, '-djpeg', '../viz/paperFigures/timestackCompare_LARGE.jpeg','-r300');
