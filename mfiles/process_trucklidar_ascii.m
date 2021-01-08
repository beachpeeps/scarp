% function process_trucklidar(filedir, filename, savedir)
% filedir = '/Volumes/group/LiDAR/Waves/20191122_TorreyPines_Runup_Test/RiProcessMethod1/';




clear
% filename = '~/Desktop/Copy of DeploymentNotes2019-2020.xls';
% P = readtable(filename,'Range','A18:O27','ReadVariableNames',false,'ReadRowNames',true);
% VarNames = readtable(filename,'Range','A1:O2','ReadRowNames',true);
% P.Properties.VariableNames = VarNames.Properties.VariableNames;

%%
filedir = '~/Documents/Repositories/scarp/data/las/truck/';
filedir = '/Volumes/FiedlerBot8000/scarp/data/las/truck/';
filedir = '../data/las/';

filename = dir([filedir '2019*.txt']);
% filename = dir([filedir 'WaveScan_*.txt']);
%%
for nfile = 1
% nfile = 1;
%%
header = {'PID','x','y','z','range','amp','time','reflectance','theta','targetIndex','targetCount'};
A = dlmread([filedir filename(nfile).name],',',1,0);
%%
% phi,x,y,z,range,amp,time,reflectance
tGPS = A(:,7);

% c = lasdata([filedir filename(nfile).name]);
% tGPS = get_gps_time(c);
[tUTC,tLST] = GPStoUTC_midnightStart(tGPS,filename(nfile).name);
% data is in 100 Hz
%%
%sort on T
[T, sortedInd] = sort(tUTC);
A = A(sortedInd,:);

PID = A(:,1);
x = A(:,2);
y = A(:,3);
z = A(:,4);
range =A(:,5);
amp = A(:,6);
reflectance =A(:,8);
theta = A(:,9);
targetIndex = A(:,10);
targetCount = A(:,11);

% quick cleanup for visualization porpoises
notspray = find(reflectance>-27);

theta = theta(notspray);



% make structure
xyzti.x = x(notspray);
xyzti.y = y(notspray);
xyzti.z = z(notspray);
xyzti.t = T(notspray);
xyzti.i = reflectance(notspray);


%% rotate around P1, and transect angle

%% put into local coordinate system, rotate around P1

Paros = load('../mat/sensors.mat','P','theta');

THETA = deg2rad(Paros.theta);
XO = Paros.P.UTMEastings_Zone11_(1);
YO = Paros.P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(xyzti.x,xyzti.y,THETA,XO,YO);
%%
clf
ind = 1:1000:length(xyzti.z);
hold on
scatter(xyzti.x(ind),xyzti.z(ind),10,xyzti.t(ind))
axis equal
%%
% clf
% ind = 1:1000:length(z);

% scatter(x(ind),z(ind),10,tGPS(ind))
% caxis([-40 15])
% view(40,35)
%%
clf
plot(XR,'.')
%%
% for ifile = 3:5
% xgrid = [20:0.1:120];

% for ifile = 9
% xgrid = [0:0.1:150];
% 
% ygrid = [780:0.1:820];

%% loop on t
% scanN = find(abs(diff(x)>60));
%% loop on scangle
% scangle = get_scan_angle(c);



%%
% scanN = find(abs(diff(scangle))>200);
% 
% for i=1:length(scanN)
% tt(i) = T(scanN(i));
% end
%%
scanN = find(abs(diff(theta))>5);
scanN = find(abs(diff(theta))>20); % dec 14 truck is THETA, not PHI (polar vs azimuthal angle)


indstart = nan(length(scanN)-1,1);
%find when to start the scan
for i=1:length(scanN)-1
    sctemp = theta(scanN(i)+1:scanN(i+1));
%     ltemp = min(sctemp);
    minloc = find(sctemp == min(sctemp),1);
    indstart(i) = scanN(i)+minloc;
end

% add first scan
sctemp = theta(1:scanN(1));
minloc = find(sctemp == min(sctemp),1);
indstart = [minloc; indstart];


%%
clf
plot(theta,'.')
hold on
plot(scanN,theta(scanN),'or')
plot(indstart,theta(indstart),'og')
% amp = get_intensity(c);
%%

tot_scans = numel(scanN);
scanDIFF =scanN-indstart;
% scanDIFF = diff(scanN);
ind_scanstart = indstart;
numPTSscan = scanDIFF;
I = tot_scans -2;
J = max(scanDIFF);

Tmat = NaN(I,J);
Rmat = NaN(I,J);
Zmat = NaN(I,J);
Amat = NaN(I,J);

for i=1:I
    Tmat(i,1:scanDIFF(i)) = xyzti.t(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Rmat(i,1:scanDIFF(i)) = XR(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Zmat(i,1:scanDIFF(i)) = xyzti.z(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Amat(i,1:scanDIFF(i)) = xyzti.i(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));

    
    ind_good = find(~isnan(Tmat(i,:)));
    Tmat(i,1:length(ind_good)) = fliplr(Tmat(i,ind_good));
    Tmat(i,1+length(ind_good):J) = NaN;
    Amat(i,1:length(ind_good)) = fliplr(Amat(i,ind_good));
    Amat(i,1+length(ind_good):J) = NaN;
    Rmat(i,1:length(ind_good)) = fliplr(Rmat(i,ind_good));
    Rmat(i,1+length(ind_good):J) = NaN;
    Zmat(i,1:length(ind_good)) = fliplr(Zmat(i,ind_good));
    Zmat(i,1+length(ind_good):J) = NaN;
end

%%
[I,J] = size(Rmat);
xi_interp = [-100:0.1:10]; % HARD CODE for region we care about

Zinterp = nan(I,length(xi_interp));
for m=3:I
    X = reshape(Rmat(m,:),1,J);
    Z = reshape(Zmat(m,:),1,J);
    [C,ia,ic] = unique(X,'stable');
    X = X(ia);
    Z = Z(ia);
    if sum(~isnan(X))>2
    Zinterp(m,:) = interp1(X(~isnan(X)),Z(~isnan(X)),xi_interp);
    end
end
clf
pcolor(xi_interp,Tmat(:,1),Zinterp); shading flat
%%
[M,~]=size(Rmat);
xi = xi_interp;
dxi=mean(diff(xi))/2;
Zinterp2=nan(M,length(xi));
Ainterp=nan(M,length(xi));
%%
for a=1:M
    r=Rmat(a,:);
    z=Zmat(a,:);
    amp=Amat(a,:);
    for i=1:length(xi)
        ind=r>=xi(i)-dxi & r<=xi(i)+dxi;
        if ~isempty(z(ind))
            Zinterp2(a,i)=median(z(ind),'omitnan');
            Ainterp(a,i)=median(amp(ind),'omitnan');
        end
    end
end

% xi_interp = -xi_interp;

%denoise
%% save data
Processed.x = xi_interp;
Processed.Zinterp = Zinterp;
Processed.Zinterp2 = Zinterp2;
Processed.Ainterp = Ainterp;
Processed.t = Tmat(:,1)';


savedir = '../mat/lidar/truck/';
% save([savedir filename(nfile).name(1:end-4)],'Processed','xyzti');
% save([savedir '20200224_00582_TorreyRunup_H' num2str(nfile+1,'%1.0f') '_10cm'],'Processed','xyzti','-v7.3');
save([savedir '20191214_00582_TorreyRunup_H' num2str(nfile,'%1.0f') '_10cm'],'Processed','xyzti','-v7.3');


%% 
hFig = figure;
pcolor(xi_interp,Processed.t(1:1:1500),Processed.Zinterp2(1:1:1500,:));
shading flat
xlabel('Cross-shore (m)')
datetick('y','MM:SS')
ylabel('Time (MM:SS)')
title(filename(nfile).name, 'Interpreter', 'none');
hc = colorbar;
caxis([0.25 2.5])
hc.Label.String = 'Z navd88 (m) ';
print(hFig, '-djpeg', ['../viz/' filename(nfile).name(1:end-4) '_timestack.jpg'],'-r300');


% %%
% hFig = figure;
% plot(Processed.t,Zinterp2(:,6)-nanmean(Zinterp2(:,6)))
% hold on
% plot(Processed.t,Zinterp2(:,11)-nanmean(Zinterp2(:,11)))
% plot(Processed.t,Zinterp2(:,16)-nanmean(Zinterp2(:,16)))
% ylim([-0.6 0.6])
% xi_interp([6 11 16])
% legend('5m','10m','15m')
% ylabel('Drift: Z-mean(Z) (m)')
% datetick('x','MM:SS')
% xlabel('Time (MM:SS')
% title(filename(nfile).name, 'Interpreter', 'none');

% print(hFig, '-djpeg', [savedir filename(nfile).name '_drift.jpg'],'-r300');


%%

% print(hFig, '-djpeg', [savedir filename '.jpg'],'-r300');

% caxis([0.5 3])
% %% 
% hold on
% load sensorlocs.mat
% 
% plot([Paros(1).x Paros(1).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(2).x Paros(2).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(3).x Paros(3).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(4).x Paros(4).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(5).x Paros(5).x],[Tmat(2,1) Tmat(end,1)],'r')
% 
% 
%
% %%
% hFig = figure
% hFig.PaperSize = [7 3];
% hFig.PaperPosition = [0 0 7 3];
% %%
% n = 1;
% wavecolor = [51 102 153]/256;
% Zinterp0 = Zinterp;
% Zinterp0(isnan(Zinterp0)) = 0;
% % for i = 1:1:500
%     i= 5;
%     
%     h1 = plot(XR(indstart(i):scanN(i)),xyzti.z(indstart(i):scanN(i)),'.b');
%     hold on
% %     h2 = plot(xi_interp,Zinterp(i,:));
%     hw = fill([xi_interp fliplr(xi_interp)],[Zinterp0(i,:)'; zeros(401,1)],wavecolor);
%     xlim([-120 0])
%     ylim([0 4])
%     pause(0.1)
%     xlabel('m from lidar')
%     ylabel('Z (m, NAVD88)')
%     drawnow
%     F = getframe(gcf);
%     [im{n},map] = frame2im(F);
% %     delete(h1)
% %     delete(h2)
% %     n= n+1
% % end
% 
% %% 
% % % %%
% currdir = pwd;
% cd /Users/juliafiedler/Documents/Torrey/viz/gifs/
% gifname = [filename(nfile).name '.gif'];
% for idx = 1:length(im)
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.01);
%     else
%         imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.01);
%     end
% end
% cd(currdir)
% 

% for i=1:length(tt)-1
% xtemp = x(scanN(i)+1:scanN(i+1));
% xu = find(unique(xtemp));
% ztemp = z(scanN(i)+1:scanN(i+1));
% xtemp = xtemp(xu);
% ztemp = ztemp(xu);
% 
% zinterp(i,:) = interp1(xtemp,ztemp,xgrid,'nearest');
end