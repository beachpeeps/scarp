function [Processed,xyzti] = process_truckanddronelidar(filedir, filename, savedir, drone, dxgrid)
% PROCESS_TRUCKANDDRONELIDAR 
% This function does the initial gridding of the lidar data, puts things
% into the local coordinate system, and saves to a mat file for future QC
% and processing. 
%
% INPUTS:
% filedir: location of the las files to look at
% filename: name of the las file in filedir
% savedir: where to save the output
% drone: 0 or 1, to say if we're processing drone or truck data
% 
% OUTPUTS:
% Processed = 
%   struct with fields:
% 
%            x: [1×441 double], x, local coordinate system
%      Zinterp: [15781×441 double], (t,x)
%     Zinterp2: [15781×441 double], (t,x)
%      Ainterp: [15781×441 double], (t,x)
%            t: [1×15781 double], time in UTC
% xyzti = 
%       struct with fields x,y,z,t,i, sorted by time.
%       coordinates are in UTM eastings and northings and z NAVD88
%       t is UTC, i is amplitude
%
% copyright Julia Fiedler 2020 jfiedler@ucsd.edu
% a fair bit of this has been gleaned from USACE's lidar processing codes
% circa 2014

%% get las object and re-organize
c = lasdata([filedir filename]);
% get time variable from las object
tGPS = get_gps_time(c);
%convert GPS time to UTC time
[tUTC,tLST] = GPStoUTC(tGPS,filename);
disp(['Processing lidar las file, LOCAL start time: ' datestr(min(tLST)) ', stop time: ' datestr(max(tLST))])
%% get xyzit and sort on time variable
[T, x, y, z, amp, sortedInd] = sortLASobject(tUTC,c);
% make structure
xyzti.x = x;
xyzti.y = y;
xyzti.z = z;
xyzti.t = T;
xyzti.i = amp;


%% rotate around P1, and transect angle
load('../mat/sensors.mat','P')
thetadeg = 352;
THETA = deg2rad(thetadeg);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(x,y,THETA,XO,YO);
%% initial look at cross-shore transect
% take subset of points for easy plotting
ind = 1:1000:length(z);
hold on
scatter(XR(ind),z(ind),10,T(ind))
xlabel('X (m)')
ylabel('Z (m, navd88)')
colorbar
title(['x-shore point cloud: ' filename])
%% initial look at scans in cross-shore
clf
plot(T,XR,'.')
xlabel('time (UTC)')
ylabel('X (m)')
title(['x over time: ' filename])

%% Start the gridding
scangle = get_scan_angle(c);
scangle = scangle(sortedInd);
scangle = double(scangle);
if drone ==1 
    % need to 
    scangle = scangle*0.006;
end

scanN = find(abs(diff(scangle))>30); % not sure how to define this
scanN = find(abs(diff(scangle))>20); % seems to work better for dec 9 truck


indstart = nan(length(scanN)-1,1);
%find when to start the scan
for i=1:length(scanN)-1
    sctemp = scangle(scanN(i)+1:scanN(i+1));
%     ltemp = min(sctemp);
if drone ==1
    maxloc = find(sctemp == min(sctemp),1);
elseif drone == 0
    maxloc = find(sctemp == max(sctemp),1);
end
    indstart(i) = scanN(i)+maxloc;
end

% add first scan
sctemp = scangle(1:scanN(1));
minloc = find(sctemp == max(sctemp),1);
indstart = [scanN(1)+maxloc-1; indstart];


clf
plot(scangle,'.')
hold on
plot(scanN,scangle(scanN),'or')
plot(indstart,scangle(indstart),'og')
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
    Tmat(i,1:scanDIFF(i)) = T(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Rmat(i,1:scanDIFF(i)) = XR(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Zmat(i,1:scanDIFF(i)) = z(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Amat(i,1:scanDIFF(i)) = amp(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));

    
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
xi_interp = [-200:dxgrid:20]; % HARD CODE for region we care about

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

pcolor(Zinterp); shading flat
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
            Zinterp2(a,i)=median(z(ind));
            Ainterp(a,i)=median(amp(ind));
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


save(['../mat/lidar/drone/' filename(1:end-4)],'Processed','xyzti');
% save(['../mat/lidar/truck/' filename(1:end-4)],'Processed','xyzti');

%% 
hFig = figure;
pcolor(xi_interp,Processed.t(1:1:3000),Processed.Zinterp2(1:1:3000,:));
shading flat
xlabel('m from lidar')
datetick('y','MM:SS')
ylabel('Time (MM:SS)')
title(filename, 'Interpreter', 'none');
hc = colorbar;
caxis([0.5 2])
hc.Label.String = 'Z ';
print(hFig, '-djpeg', ['../viz/' filename '_timestack.jpg'],'-r300');


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


% for i=1:length(tt)-1
% xtemp = x(scanN(i)+1:scanN(i+1));
% xu = find(unique(xtemp));
% ztemp = z(scanN(i)+1:scanN(i+1));
% xtemp = xtemp(xu);
% ztemp = ztemp(xu);
% 
% zinterp(i,:) = interp1(xtemp,ztemp,xgrid,'nearest');
% end