clear
filename = '~/Desktop/Copy of DeploymentNotes2019-2020.xls';
P = readtable(filename,'Range','A18:O27','ReadVariableNames',false,'ReadRowNames',true);
VarNames = readtable(filename,'Range','A1:O2','ReadRowNames',true);
P.Properties.VariableNames = VarNames.Properties.VariableNames;
%%
clf
thetadeg = 352;
X = P.UTMEastings_Zone11_;
Y = P.UTMNorthings_Zone11_;
THETA = deg2rad(thetadeg);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(X,Y,THETA,XO,YO);

scatter3(XR,YR,P.DeploymentNAVD88_m_Geoid12B_,50,'k');
hold on
text(XR+1,YR-1,P.DeploymentNAVD88_m_Geoid12B_,P.Properties.RowNames)

%%
% Get the truck lidar
% clear
filedir = '~/Desktop/';
filename = dir([filedir '*.las']);

nfiles = [1 2 3 4];

for nfile = nfiles
    cvar = 'r';
    hlidar(nfile) = plot_lidarData(nfile,cvar,filedir, filename,THETA,XO,YO);
end
view(0,90)

%%
xlabel('Cross-shore (m), local coordinate system, P1 = Origin, 262 deg')
ylabel('Along-shore (m)')


cmap = lines(8);
n = 1;
for i=nfiles
    hlidar(i).CData = cmap(n,:);
n = n+1;
end

hlidar(28) = plot_lidarData(28,'m',filedir, filename,THETA,XO,YO);


legend(hlidar([ nfiles 28]),filename([nfiles 28]).name,'interpreter','none')
