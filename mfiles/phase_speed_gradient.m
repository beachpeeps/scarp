% load('20200224_00582_TorreyRunup_H2.mat')

x = x(1:xlimit,1);
z =eta(1:xlimit,:); % interpolated data
t = datenum(tstack.tvecHover);


%make mask of non-interpolated values
zmask = tstack.TXdrone(1:xlimit,:); % non-interpolated data
zmask(~isnan(zmask)) = 1;

%time in seconds, t = 0 is mean of record
ts = (t-nanmean(t))*24*60*60;

%array should be x: rows, t: columns
nt = length(t);
nx = length(x);
[m,n] = size(z);
if (m~=nx)
    z = z';
    zmask = zmask';
end
%sample spacing and period
dx = mean(diff(x));
dt = nanmean(diff(ts));

%compute the gradient of the surface elevation
[gt,gx] = gradient(z);
%gradient amplitude
ga = sqrt(gt.^2+gx.^2);

%find gradient above threshold, directed normal to crests
k = ga > 0.05 & gx < 0 & gt > 0;
% k = gx < 0 & gt > 0;

gg = NaN(size(ga));
gg(k) = ga(k);
ggx = NaN(size(ga));
ggt = ggx.*zmask;
ggx(k) = gx(k);
ggt(k) = gt(k);
%%
%%%%%%%%%%%%%%%
%compute depth: time-average all variables before computing depth
%%%%%%%%%%%%%%%

% ggxa = nanmedian(ggx')';
% ggta = nanmedian(ggt')';

%angle of gradient
% ang_avg = atan2(ggxa,ggta); % average angle
ang = atan2(ggx,ggt); %get angle at each point

%angle parallel to x-t crest lines
ang = ang+pi/2; 
anga = nanmedian(ang')'; % time average all angles

%phase speed and bottom depth
cphase = dx/dt*tan(ang); %<-- here tan applied again because of rotation
cphase_preaverage = dx/dt*tan(anga);

%%

depth = (cphase.^2)./9.8;
depth_postA = nanmedian(depth,2); % average after all calcs are done

depth_preA = (cphase_preaverage.^2)./9.8; % pre averaged depth

twindow = floor(1./fp/2);

%depth with bore height correction from Martins et al. (2017)
%get bore height
H = NaN(size(z));
for j = 1:nt % for each timestep
    kj = find(~isnan(gg(:,j))); % find 'wave' points
    nkj = length(kj); % get number of wave points per timestep
    if(nkj > 0)
        zz = z(:,j); %surface at time step j
        for j2 = 1:nkj % for each wave point
%             H(kj(j2),j) = range(zz(max(1,kj(j2)-10):min(kj(j2)+10,nx)));

if j>window_length && j<nt-window_length && ~isnan(twindow(j2))
            H(kj(j2),j) = range(z(kj(j2),j-twindow(j2):j+twindow(j2)));
end

%                 %find local minimum
%                 if kj(j2)+2*window_length<xlimit
%               H(kj(j2),j) = range(zz(kj(j2):kj(j2)+2*window_length)) ;
%                 else
%                    H(kj(j2),j) = range(zz(kj(j2):kj(end))) ;
%                 end
        end
    end
end
depth2_preA = (cphase_preaverage(1:xlimit).^2)./9.8-nanmedian(H')'/2;
depth2 = (cphase.^2)./9.8-H./2;
depth2_postA = nanmedian(depth2')';



%%
figure
clf
xchunk = 1:xlimit;
tchunk = wave(1).t(2)-50:wave(1).t(end)+50;
pcolor(z(xchunk,tchunk))
shading flat
ylabel('cross-shore index (dx = .1m)')
xlabel('time index (dt = 0.1 s)')
set(gca,'FontSize',14)
g = colorbar;
ylabel(g,'elevation (m)')
% axis('equal')
xlim([0 length(tchunk)])
hold on
% print -dpng fig1a.png
quiver(ggt(xchunk,tchunk),ggx(xchunk,tchunk),10,'r')
plot(wave(1).t-tchunk(1),10*wave(1).xpeak+1000,'ok')
% axis('equal')
% print -dpng fig1b.png

%%
clf
tstep = 5;
ti = 1;
while ti<=length(wave(1).t)
if mod(ti,tstep) == 0 
            h1 = plot(x(1:xlimit),eta(1:xlimit,wave(1).t(ti))+ti/tstep,'color',[0 0 0 0.5]); % wave
            hold on
            
            kk = x>=wave(1).xpeak(ti) & x<=wave(1).xtrough(ti);
            h3(1) = plot(wave(1).xpeak(ti),wave(1).hpeak(ti)+ti/tstep,'ok','markersize',10); %roller face
%             h1(2) = plot(wave.x(xi:xtroughInd(ti)),eta(xi:xtroughInd(ti),ti)+ti/tstep,'c','linewidth',2);
            h3(2) = plot(wave(1).xtrough(ti),wave(1).htrough(ti)+ti/tstep,'sqk','markersize',10); %roller face
            
%             h3(1) = plot(x(kk),eta(kk,wave.t(ti))+ti/tstep,'k','linewidth',3); % wave
            xplot = find(~isnan(ggt(:,wave(1).t(ti))));
            h3(3) = scatter(x(xplot),eta(xplot,wave(1).t(ti))+ti/tstep,50,cphase(xplot,wave(1).t(ti)),'filled');
            
end
        ti = ti+1;
end

caxis([0 5])
ylim([9 31])
xlim([-80 -30])
hc = colorbar;
ylabel('timestep (0.5s)')
ylabel(hc,'C_p (m/s)')
xlabel('x (m)')
title([hoverdate ' H' num2str(hovern) ' crestTrack_truck.png'])
print(gcf,'-dpdf',['../viz/' hoverdate '_H' num2str(hovern) '_crestTrack.pdf'],'-bestfit')

%%
% xlim([100 130])
% ylim([-60 -5])

figure 
clf
kk = find(x<-5);
plot(x(kk),-depth_preA(kk),'LineWidth',1)
xlabel('cross-shore (m)')
ylabel('depth (m)')
ylim([-6 0])
set(gca,'FontSize',14)
title('Lidar Bathy, hover 2')
ylim([-6 0])
print -dpng fig2a.png

hold on
plot(x(kk),-depth2_preA(kk),'LineWidth',1)
legend('raw','w/ bore height correction')
print -dpng fig2b.png
% xlim([-65 -10])
ylim([-2.5 0])
print -dpng fig3.png


figure
zmin = nanmin(z')';
plot(x,zmin,x,-depth2_preA+1.7,'LineWidth',2)
set(gca,'FontSize',14)
legend('Lidar bathy','Lidar min elevation')
% xlim([-65 10])
ylim([0 5])
ylabel('elevation (m)')
xlabel('cross-shore (m)')
print -dpng fig4.png

%%
figure
r = nan(1,xlimit);
for i=1:xlimit
    if sum(~isnan(cphase(:,i)))>10
        val = find(~isnan(cphase(:,i)) & cphase(:,i)<6);
        
        mdl = corrcoef(H(val,i),cphase(val,i));
        if ~isnan(mdl(1)) && numel(mdl)>1
        r(i) = mdl(2);
        end
    end
end

hgrad = plot(x(1:xlimit),r,'x');
hold on
plot(x, zeros(size(x)),':k')
xlabel('x (m)')
ylabel('Correlation r')
title('Correlation of C_p and H')
print -dpng fig5.png

