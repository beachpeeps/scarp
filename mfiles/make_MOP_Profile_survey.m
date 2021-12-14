% function bathy = make_MOP_Profile_survey(bathyfile,MopID,dx,datevec,outfileDir,bathyname)
% This function uses Bill's CPGMOPS codes/data on the reefbreak server to
% generate average cross-shore profiles for a given mop number.
%
% Input:
% bathydir = directory where the gridded bathymetry resides
% MopID = number of MOP line
% dx = spacing for bathymetry
% datevec = vector of time period, e.g. [dateStart dateEnd] for when to
% produce bathy
%
% Output:
% bathy structure with x, z (navd88), dx, date created, and months used
%
% 2021, Julia Fiedler
Mopnum = 582; % desired mop number
county = 'D';
bathyname = 'mopcpg_scarp';
baseDir = '~/Documents/';


dataFolder = [baseDir 'IPA_' county num2str(Mopnum,'%04.0f') '/data/']; 
outfileDir = dataFolder;


addpath /Volumes/group/MOPS/toolbox
addpath /Volumes/group/MOPS

load('MopTableUTM.mat','Mop');  % Load "Mop" table array
% 
% if nargin<6
%     minmax = [];
% end
bathyfile = '~/Documents/Repositories/scarp/mat/scarp_bathy.mat';
load(bathyfile)
dx = 1;
xx = repmat(x,90,1);
% bathy.name = strcat('mean profile for months: ', month(datetime(1, mon(1), 1), 's'), '-', month(datetime(1, mon(end), 1), 's'));

%bring in riprap data from SCaRP experiment
riprapfile = '/Volumes/jfiedler/scarp/mat/timestacks_LARGE/20191214_1.mat';
load('/Users/juliafiedler/Documents/Repositories/scarp/mat/sensors.mat','P')

load(riprapfile,'Xgrid','TXdrone')
riprap = nanmin(TXdrone');
xriprap = Xgrid(:,1);
ind = find(xriprap>-10 & xriprap<10);
riprap = riprap(ind);
xriprap = -xriprap(ind);

%move this origin to the backbeach origin of the MOP line instead of P1
%from SCaRP
P0 = [P.UTMEastings_Zone11_(1)];
Mop0 = [Mop.BackXutm(Mopnum)];
xriprap = xriprap-(P0-Mop0)-5; %something is off in the alignment, but P1 was at 582.3 so maybe that's the issue???

bb = interp1(xriprap,riprap,x);

kk = x>=-3 & x<=3;
ZZfill(:,kk) = repmat(bb(kk),90,1);
zmid = smoothdata(ZZfill,2,'gaussian',15/dx);
zup = smoothdata(ZZfill,2,'gaussian',3/dx);

kk =  xx<=25;
zmid(kk) = zup(kk);
z = zmid;

clf
plot(x,bb)
hold on
plot(x,z(6,:))
% xlim([-10 300])
% 
% %check
% % for i=1:length(t)
% % plot(x,z(i,:))
% % xlim([0 300])
% % ylim([-10 5])
% % pause(0.1)
% % end
% 
% 
% extend to z=10
for i=1:90
z(x<-3) = nan;
ind = find(~isnan(z(i,:)'),1,'first');
mdl = fitlm(x(ind:ind+5),z(i,ind:ind+5));
mb = mdl.Coefficients.Estimate;
% 
newslope = mb(1)+x*mb(2);
xind = find(newslope(:)<10,1,'first');
% 
z(i,xind:ind) = newslope(xind:ind);
end

x = -x;

xx1 = [min(x):2:-202 -200:1:-50 -49:0.5:-29.5 -29:0.25:x(xind)];

for i=1:90
zz(i,:) = interp1(x,z(i,:),xx1);
end

for i=1:90
% 
bathy(i).x = xx1;
bathy(i).z = zz(i,:);
bathy(i).dx = 'variable from 2 to 0.25m';

bathy(i).dateCreated = datetime;
bathy(i).t = t(i);
bathy(i).Mop = Mopnum;
end
% 
save([outfileDir bathyname '.mat'],'bathy')
disp(['bathymetry created and saved to ' outfileDir bathyname '.mat'])

rmpath /Volumes/group/MOPS/toolbox
rmpath /Volumes/group/MOPS