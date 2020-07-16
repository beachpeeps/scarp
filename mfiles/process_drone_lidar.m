% process_drone_lidar
% needs to be cleaned up, but uploading to github as is, in most
% embarassing form.
% loosely based on ERDC FRF's lidar processing code from 2014
% April 21,2020 Julia Fiedler jfiedler@ucsd.edu
filedir = '/Volumes/FiedlerBot8000/scarp/data/las/drone/';
filename = dir([filedir '20200224*.las']);
%%
for nfile = 1:5
c = lasdata([filedir filename(nfile).name]);
t = get_gps_time(c);

% data is in 100 Hz
%%
%convert t to datenum
secondsinoneday = 60*60*24;
sunstart = datenum(2019,12,8,0,0,0);
t_in_days = t./secondsinoneday;
T = sunstart+t_in_days;

[T,I] = sort(T);
%% get values
x = c.x;
y = c.y;
z = c.z;
intensity = get_intensity(c);

x = x(I);
y = y(I);
z = z(I);
zz = z;
amp = intensity(I);

% cutoff bad data
% TO DO: CAN WE find this automatically? This is not good programming.
% newInd = 1.3e6:length(x);

% x = x(newInd);
% y= y(newInd);
% z = z(newInd);
% zz = z;
% amp = amp(newInd);
% T = T(newInd);



% make structure
xyzti.x = x;
xyzti.y = y;
xyzti.z = z;
xyzti.t = T;
xyzti.i = amp;

%% plot data in 3D
% view(40,35)

ind = 1300000:1000:length(x);
scatter3(x(ind),y(ind),zz(ind),30,T(ind))



%% put into local coordinate system

load('../mat/sensors.mat')

THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[xr YR] = xyRotate(x,y,THETA,XO,YO);
%%
clf
plot(y,'.')

%% loop on r
scanN = find(diff(y)<-4);
%% loop on scangle
scangle = get_scan_angle(c);



%%
% scanN = find(abs(diff(scangle))>200);
% 
% for i=1:length(scanN)
% tt(i) = T(scanN(i));
% end
%%
% scangle = get_scan_angle(c);
scangle = scangle(I);
scangle = scangle(newInd);
%%
% scangle = double(scangle);
% scanN = find(abs(diff(scangle))>5);

indstart = nan(length(scanN)-1,1);
%find when to start the scan
for i=1:length(scanN)-1
    sctemp = y(scanN(i):scanN(i+1));
%     ltemp = min(sctemp);
    minloc = find(sctemp == min(sctemp),1);
    indstart(i) = scanN(i)+minloc-1;
end

% add first scan
sctemp = y(1:scanN(1));
minloc = find(sctemp == min(sctemp),1);
indstart = [scanN(1)+minloc-1; indstart];

%%
clf
plot(y,'.')
hold on
plot(scanN,y(scanN),'or')
plot(indstart,y(indstart),'og')
% amp = get_intensity(c);
%%

tot_scans = numel(scanN);
% scanDIFF =scanN-indstart;
scanDIFF = diff(scanN);
ind_scanstart = indstart;
numPTSscan = scanDIFF;
I = tot_scans -2;
J = max(scanDIFF);

Tmat = NaN(I,J);
Rmat = NaN(I,J);
Zmat = NaN(I,J);
Amat = NaN(I,J);

for i=40:I
    Tmat(i,1:scanDIFF(i)) = T(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Rmat(i,1:scanDIFF(i)) = xr(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
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
xi_interp = [-200:0.25:50]; % HARD CODE for region we care about

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

%%
[M,~]=size(Rmat);
xi = xi_interp;
dxi=mean(diff(xi))/2;
Zinterp2=nan(M,length(xi));
Ainterp=nan(M,length(xi));
%%
f = waitbar(0,'Getting Zinterp2');

for a=1:M
    waitbar(a/M,f)
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

filedir = '../data/las/drone/';
filename = dir([filedir '*.las']);


savedir = '../mat/lidar/drone/';
save([savedir filename(nfile).name(1:end-4) '.mat'],'Processed','xyzti');
%% 
hFig = figure;
pcolor(xi_interp,Processed.t(1:1:end),Processed.Zinterp2(1:1:end,:));
shading flat
xlabel('m from lidar')
datetick('y','MM:SS')
ylabel('Time (MM:SS)')
title(filename(nfile).name, 'Interpreter', 'none');
hc = colorbar;
caxis([2 6])
hc.Label.String = 'Z ';
% print(hFig, '-djpeg', [savedir filename(nfile).name(1:23) '_timestack.jpg'],'-r300');


%%
hFig = figure;
plot(Processed.t,Zinterp2(:,200)-nanmean(Zinterp2(:,200)))
hold on
plot(Processed.t,Zinterp2(:,150)-nanmean(Zinterp2(:,150)))
plot(Processed.t,Zinterp2(:,300)-nanmean(Zinterp2(:,300)))
ylim([-0.06 0.06])
xi_interp([200 150 300])
legend('20m','15m','30m')
ylabel('Drift: Z-mean(Z) (m)')
datetick('x','MM:SS')
xlabel('Time (MM:SS')
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
% for i = 1001:1:1500
%     
%     h1 = plot(x(indstart(i):scanN(i)),zz(indstart(i):scanN(i)),'.k');
%     xlim([-450 -300])
%     ylim([-34 -31])
%     pause(0.1)
%     xlabel('m from lidar')
%     ylabel('Z (m, NAVD88)')
%     drawnow
%     F = getframe(gcf);
%     [im{n},map] = frame2im(F);
%     delete(h1)
%     n= n+1
% end
% 
% % 
% %%
% currdir = pwd;
% % cd /Users/juliafiedler/Documents/LidarAdam/gifs/
% gifname = [filename.name '.gif'];
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
% end