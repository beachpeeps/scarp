clear
% filename = '~/Desktop/Copy of DeploymentNotes2019-2020.xls';
% P = readtable(filename,'Range','A18:O27','ReadVariableNames',false,'ReadRowNames',true);
% VarNames = readtable(filename,'Range','A1:O2','ReadRowNames',true);
% P.Properties.VariableNames = VarNames.Properties.VariableNames;
load('../mat/sensors.mat')
figureDir = '../viz/';
figureName = 'transect3DwithLidar_2lines';
%%
figwidth = 8;
figheight = 6;
nrow = 1;
ncol = 1;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,'units','inches');

X = P.UTMEastings_Zone11_;
Y = P.UTMNorthings_Zone11_;
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(X,Y,THETA,XO,YO);

scatter3(XR,YR,P.DeploymentNAVD88_m_Geoid12B_,50,'k');
hold on
text(XR+1,YR-1,P.DeploymentNAVD88_m_Geoid12B_,P.Properties.RowNames)

%%

% Get the truck lidar
% clear
filedir = '../data/las/truck/lastPoints/';
filename = dir([filedir '*.las']);
cmap = lines(8);
% nfiles = [1 2 3 4];
hold on
% for nfile = nfiles
nfile = 5;
    cvar = cmap(1,:);
    hlidar(nfile) = plot_lidarData(nfile,cvar,filedir, filename,THETA,XO,YO);
% end


%%
xlabel('Cross-shore (m), local coordinate system, P1 = Origin, 262 deg')
ylabel('Along-shore (m)')

hold on

% n = 1;
% for i=nfiles
%     hlidar(i).CData = cmap(n,:);
% n = n+1;
% end

filedir = '../data/las/drone/';
filename = dir([filedir '*.las']);

hlidar(2) = plot_lidarData(3,cmap(2,:),filedir, filename,THETA,XO,YO);


legend([hlidar(5) hlidar(2)],{'truck', 'drone'})
%%
view(-25,55)
print(hFig, '-djpeg', [figureDir figureName '.jpeg'],'-r300');
