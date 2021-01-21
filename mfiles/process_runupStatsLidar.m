clear
% rootFolder = '/Volumes/FiedlerBot8000/scarp/';
% folderDir = [rootFolder 'mat/lidar/drone/'];
% folderDir = [rootFolder 'mat/lidar/truck/'];

% filename = dir([folderDir 'Final*.mat']);
% filename = dir([folderDir '20191214*_10cm.mat']);
% filename = dir([folderDir '20200224*_10cm.mat']);

% hoverdate = '20200224';
hoverdate = '20191214';


datadir = '/Volumes/FiedlerBot8000/scarp/';
options.drone = 1;


% filename = dir([folderDir 'Wave*.mat']);
% load('../mat/20200224/DroneStartStop_20200224.mat','tstart','tstop')
% load('../mat/20191214/DroneStartStop_20191214.mat','tstart','tstop')


for i=[1 3 4]
% for i=[2 4]


tstack = [datadir '/mat/timestacks/' hoverdate '_' num2str(i) '.mat'];
% [Spec,Info,Bulk,Tseries] = get_runupStatsLidar_L1([folderDir filename(i).name],'tchunk',[tstart(i) tstop(i)]);
[Spec,Info,Bulk,Tseries] = get_runupStatsLidar_L1(tstack,options);

if options.drone == 1
processedFilename = ['../mat/' hoverdate '/Drone_Hover_' num2str(i,'%02.0f') '_L1'];
elseif options.drone == 0
processedFilename = ['../mat/' hoverdate '/Truck_Hover_' num2str(i,'%02.0f') '_L1'];
end


save([processedFilename '_runupstats_10cm.mat'],'Spec','Info','Bulk','Tseries')
    
end


