% track wave crests
clear all

for hovern = [1 3:4]
    clearvars -except hovern
hoverdate = '20191214';
% hoverdate = '20200224';

% switch the limits for analysis, depending on hoverdate, to stop at x=0
if str2double(hoverdate) == 20191214
    xlimit = 1000;
elseif str2double(hoverdate) == 20200224
    xlimit = 500;
end

datadir = '/Volumes/FiedlerBot8000/scarp/';

tstack  = load([datadir '/mat/timestacks/' hoverdate '_' num2str(hovern) '.mat']);

%% decide if you want to look at truck or drone data
eta =tstack.TXdrone2-nanmean(tstack.TXdrone2,2); % interpolated
eta2 = tstack.TXdrone-nanmean(tstack.TXdrone2,2); % not interpolated

% eta = tstack.TXtruck2-nanmean(tstack.TXtruck2,2);
% eta2 = tstack.TXtruck-nanmean(tstack.TXtruck2,2);

[nx,nt] = size(eta);
x = tstack.Xgrid(:,1);
dx =x(2)-x(1);

%%
% CC =eta(1:xlimit,:);
% XX = Xgrid(1:xlimit,:);
% TT = datenum(Tgrid(1:xlimit,:));
%
% [crestmag,crestdir] = imgradient(CC,'sobel');
% threshold = 0.03;
% [cmask,threshold] = edge(CC,'sobel','nothinning',threshold);
%
% % imshowpair(crestmag(:,1:1000),crestdir(:,1:1000),'montage')
%
% % make masks
% % only want the incoming waves, so make sure direction is positive
% cmask(crestdir<5) = 0;
% % remove outliers
% cmask(crestdir>60) = 0;
% cmask(abs(crestdir)>(2*nanstd(crestdir(:))+nanmedian(crestdir(:)))) = 0;
%
% cmask = double(cmask);
% cmask(cmask==0) = nan;
%
% %make sure that we aren't using interpolated data for this
% etamask = TXtruck;
% etamask(~isnan(etamask)) = 1;
%
% % crest directions
% % mask out everything but crests
% cdir = crestdir.*cmask.*etamask(1:xlimit,:);
%
% % use median instead of mean - outliers will otherwise skew answers
% cdirmean = nanmedian(cdir,2);
%
%
% % get bore heights
%
% % this will make the range filter emphasize the x axis (we want H as
% % measured in x for every timestep t.)
% nhood = ones(31,1);
% H= rangefilt(CC,nhood);
% Hh = H(1:xlimit,:).*cmask.*etamask(1:xlimit,:);
% etam = eta2(1:xlimit,:).*cmask.*etamask(1:xlimit,:);
%
% % get phase speeds of wave crests
% C = 1./tand(cdir); % per (x,t)
% Cmean = 1./tand(cdirmean); % per (x)
%
% %remove outliers
% C(abs(C)>(nanstd(C(:))+nanmedian(C(:)))) = NaN;
%
% % these will be rewritten in the crest tracking code, if you want to use
% % them just comment out the "crest tracking" section
% h = Cmean.^2./9.81;
% h2 = Cmean.^2./9.81 -nanmedian(Hh,2)./2;

%% crest tracking
% get fp
addpath ~/Documents/MATLAB/jlab/
pp = 12;
dt = 1/tstack.Hz_lidar;
[PSI,~] = sleptap(length(tstack.tvecHover),pp);


fp = nan(nx,1);
[f,Snn(1,:)]=mspec(dt,detrend(eta(1,:))',PSI);
SS = find(f>0.04 & f<0.25);
for i=2:nx
    [f,Snn(i,:)]=mspec(dt,detrend(eta(i,:))',PSI);
    fp(i) = f(SS(Snn(i,SS) == max(Snn(i,SS))));
end
rmpath ~/Documents/MATLAB/jlab/

%%
% locate the wave crests at the edge of the domain
xloc_startcrest = round(2/dx); % 1m from edge of domain
% xloc_startcrest = 1; % edge of domain

window_t = floor(tstack.Hz_lidar/fp(xloc_startcrest)/2); % get timechunk for finding peaks at xloc_startcrest

% get elevation data at edge of domain
z = eta(xloc_startcrest,:);

% re-organize elevation data
length_z = nt - mod(nt,window_t);
zz = reshape(z(1:length_z),window_t,length_z/window_t);

% find max points for each timechunk
indPeak = find(zz == max(zz));

%  quick plot to see if method works
plot(z)
hold on
plot(indPeak,z(indPeak),'ok')
title('check peak-finding at edge of domain')

% nudge the indPeak back a bit in time to account for xloc_startcrest
indPeak(z(indPeak)<0) = [];
%% crest-tracking.
% For each wave crest identified at outer edge of domain, track the crest
% as it comes onshore, by time-stepping and looking for the peak location
% in a sliding cross-shore window around the previous crest location.


figure

window_length = round(5/dx); % window size in meters, used by Martins et al.
max_time = 500; % maximum time value to track wave (arbitrary high number)
nwave = length(indPeak);

for nw = 1:nwave-1
    
    %pre-allocate
    xx = nan(1,max_time); % all x-values of the crest per each timestep
    H = nan(1,max_time); % wave height (crest-trough) per each timestep
    
    
    % plug in starting values to get algorithm chugging
    tti = 1; %starting time index (relative to wave start)
    xi = 1; %starting x-index
    window_x = xi:xi+window_length/2; % starting x window indices
    
    
    % for each timestep, begining at start of wave
    for ti=indPeak(nw):indPeak(nw)+max_time 
        
        
        % find the index of maximum elevation within the cross-shore window
        xi = window_x(eta(window_x,ti) == max(eta(window_x,ti)));
        xi = xi(end); % in case there are multiple max points   
        
%         if tti>1 && x(xi)<xx(tti-1)
%             x(xi) = nan;
%             continue
%         end
        
        
        
        % if there's no forward progress of the search window in x
        if tti>15 && mean(xx(tti-15:tti-1))>x(xi)
            
            % account for a bore, get most shoreward point of wave crest
            xi = window_x(find(eta(window_x,ti)>0.8*max(eta(window_x,ti)),1,'last'));
            
            % stop crest tracking if the surface is flat
            % or if there is no forward progress
            % or if there is too much forward progress
            if isempty(xi) || tti>15 && mean(xx(tti-15:tti-1))>x(xi) || x(xi)>(xx(tti-1)+dx.*window_length/3)
                break
            else
            end
            
        end
        
        %stop crest tracking at the shoreline if somehow it's still going
        if xi> xlimit
            disp('pau')
            break
        end
        
        %get the actual x-location for each timestep 
        xx(tti) = x(xi);
        
        %get the height of the wave crest for each timestep 
        hpeak(tti) = eta(xi,ti);
        
        %find the location and height of the wave trough
        xsearch = xi:xi+window_length/2; %only search in front of the wave crest
        htrough(tti) = min(eta(xsearch,ti)); %minimum value is the trough height
        
        % get the most shoreward location of trough value (good for
        % plotting only, not used in any calculation.)
        xtrough(tti) = max(x(xsearch(eta(xsearch,ti) == min(eta(xsearch,ti))))); 
        
        
        H(tti) = hpeak(tti)-htrough(tti); %wave height is crest-trough
        
        % plot a timestack of the waves and the roller faces (only use every 3rd timestep)
        tstep = 1;
        if mod(tti,tstep) == 0 && nw == 13 %use the 40th wave (just to show how it works!)
            h1 = plot(x(1:xlimit),eta(1:xlimit,ti)+tti/tstep,'b'); % wave
            hold on
            h3(1) = plot(x(xi),hpeak(tti)+tti/tstep,'ok'); %roller face
            h1(2) = plot([x(xi) xtrough(tti)],[hpeak(tti) htrough(tti)]+tti/tstep,'linewidth',2); %roller face
        end
        
        
        % move to next timestep, establish new cross-shore window
        tti = tti+1;
        if xi<26
            window_x = xi:xi+window_length/2;
        else
            window_x = xi-window_length/5:xi+window_length/2;
        end
        
    end
    
    % record each wave's x-loc, wave-height, and time.
    wave_x(nw,:) = xx;
    wave_h(nw,:) = H;
    wave_t(nw,:) = indPeak(nw):indPeak(nw)+max_time-1;
    
end

xlabel('x (m)')
ylabel('timestep (0.1s)')
title([ hoverdate ' Hover ' num2str(hovern) ': Sample timestack for crest tracking, wave 13'])
clear H xx

%% disregard waves without enough data points (COULD BE FIXED)
minSec = 10/dt; % minimum number of seconds we need for crest-tracking
lwave = sum(~isnan(wave_x),2);
wave_x(lwave<minSec,:) = [];
wave_h(lwave<minSec,:) = [];
wave_t(lwave<minSec,:) = [];
[nwave,ntime] = size(wave_x);
%%
figure
tchunk = 500:2500;
pcolor(x(1:xlimit),dt*(tchunk),eta(1:xlimit,tchunk)')
shading flat
hold on
%  plot(wave_x(1:50,:),dt.*wave_t(1:50,:),'r')
for i=1:nwave
    plot(wave_x(i,:),dt.*wave_t(i,:),'r','linewidth',2)
    hold on
end
xlabel('x (m)')
ylabel('time (s)')
caxis([-1 1])
hc = colorbar;
ylabel(hc,'elevation (m)')
title([ hoverdate ' Hover ' num2str(hovern) ': wave trajectories, n = ' num2str(nwave)])

%% get moving slope
Cp = nan(nwave,xlimit);
Cp_H = nan(nwave,xlimit);
%


for i =1:nwave
    
    wavex = wave_x(i,~isnan(wave_x(i,:)));
    waveh = wave_h(i,~isnan(wave_x(i,:)));
        %
    wavet = 1:length(wavex);
    [a,b] = unique(wavex);
    %
    wt = interp1(a,wavet(b),x);
    wh = interp1(a,waveh(b),x);
   
    nwaveplot = 50;
    n = 1;
    while n < find(x == max(wavex))
        nn = n:n+window_length-1; %<-- this is a 5m window, and doesn't change. Issue?
        
        % get subset of x,t points in moving window
        xloc = x(nn);
        tloc = wt(nn);
        
        %analyze only chunk of wave with data
        xloc = xloc(~isnan(tloc));
        tloc = tloc(~isnan(tloc));
        
        %make sure there are enough points to regress to
        if numel(xloc)<20
            break
        end
        
        % simple linear regression to get slope
        E = [ones(size(tloc)) tloc];
        c = E\xloc;
        
        Cp(i,n) = tstack.Hz_lidar*c(2); % wave speed is the slope
        Cp_H(i,n) = nanmedian(wh(nn)); % get median wave height of subwindow
%         if i==nwaveplot
%             
% %             pause(0.01)
%             delete(h1)
%             h1 = plot(xloc,tloc,'.r','parent',ax1);
%             
%             plot(x(n),Cp(i,n),'.r','parent',ax2)
%             
%         end
        n = n+1;
        
        
        
        
    end
end

% get wave speed from moving slope
Cp(Cp<=0) = nan;
Cp_H(isnan(Cp)) = nan;
Cpm = nanmedian(Cp); %average in x
h_linear = Cpm.^2./9.81;
h_bore = nanmedian(Cp.^2./9.81-Cp_H./2);


%%
figure
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

hold(ax1,'on')
    hold(ax2,'on')
    hold(ax3,'on')


ax1.XLim = [x(1) x(end)];
ax2.XLim = ax1.XLim;
ax3.XLim = ax1.XLim;
plot(wave_x,1:max_time,'color',[0.8 0.8 0.8 0.4],'parent',ax1);

hwave = plot(wave_x(nwaveplot,:),1:max_time,'r','parent',ax1);
axes(ax1)
legend(hwave,['wave traj.' num2str(nwaveplot)])



hcp = plot(x(1:xlimit),Cp,'.','parent',ax2,'color',[0.8 0.8 0.8 0.4]);
h1 = plot(x(1:xlimit),Cp(nwaveplot,:),'r','parent',ax2);

hm = plot(x(1:xlimit),Cpm,'m','parent',ax2);

axes(ax2)
legend([hcp(1) h1 hm],{'C_p',['C_p wave ' num2str(nwaveplot)],'median C_p'})

plot(x(1:xlimit),-h_linear,'parent',ax3);

plot(x(1:xlimit),-h_bore,'parent',ax3);
axes(ax3)
legend('linear','with bore correction')

ax1.Title.String = ['wave trajectories and slope analysis region for wave ' num2str(nwaveplot)];
ax2.YLabel.String = 'phase speed (m/s)';
ax3.XLabel.String = 'x (m)';
ax3.YLabel.String = 'h (m)';
ax3.YLim = [-3 0];
ax3.Legend.Location = 'southeast';
ax1.XLim = [x(1) x(end)];
        ax2.XLim = ax1.XLim;
        ax2.YLim = [0 10];


clear nx nt

%%

%
%
%
% %%
% %make an empty matrix to store wave information
% wave_h = nan(length(indPeak)-10,300);
% wave_x = wave_h;
% wave_t = wave_h;
% %
%
% for waven=1:length(indPeak)-10
%     timen=indPeak(waven);
%     L = 1:100;
%
%
%     for n=1:300
%         ind = find(eta(L,timen) == max(eta(L,timen))); %<-- switch
%
%         ind = L(ind);
%         L = [ind-20:ind+20];
%         if L(end)>=xlimit
%             L = L(1):xlimit;
%         end
%         if L(1)<=1
%             L=1:L(end);
%         end
%
%         ind2 = find(eta(L,timen) == max(eta(L,timen)));
%         %         if ~isempty(ind2)
%         %             break
%         %         end
%         ind2 = L(ind2);
%
%         if n>15 && nanmean(wave_x(waven,n-15:n-1))>Xgrid(ind2)
%             break
%         end
%
%         wave_x(waven,n) = Xgrid(ind2);
%
%         if n>15 && wave_x(waven,n)-nanmean(wave_x(waven,n-15:n-1))<=0
%             break
%         end
%
%         if ind2>2 && ind2+2<xlimit
%             htemp = TXdrone2(ind2-2:ind2+2,timen-2:timen+2);
%         elseif ind2<2
%             htemp = TXdrone2(1:ind2+2,timen-2:timen+2);
%         elseif ind2+2>xlimit
%             htemp = TXdrone2(ind2-2:xlimit,timen-2:timen+2);
%         end
%
%         wave_h(waven,n) = nanmedian(htemp(:));
%         wave_t(waven,n) = timen;
%
%         L = ind2-5:ind2+5;
%         if L(end)>=xlimit
%             L = L(1):xlimit;
%         end
%         if L(1)<=1
%             L=1:L(end);
%         end
%
%         timen = timen+1;
%     end
%
%     % get new wave
%
% end
% %
% %disregard waves that don't have enough data points (COULD BE FIXED)
% lwave = sum(~isnan(wave_x'));
% wave_x(lwave<80,:) = [];
% wave_h(lwave<80,:) = [];
% wave_t(lwave<80,:) = [];
% [nwave,ntime] = size(wave_x);
% %
% %interp to same Xgrid
% n = 1;
% for i=1:nwave
%     xtemp = wave_x(i,:);
%     ttemp = wave_t(i,:);
%     xtemp = xtemp(~isnan(xtemp));
%     ttemp = ttemp(~isnan(xtemp));
%
%     [xtemp,tind] = unique(xtemp);
%     ttemp = ttemp(tind);
%
%     if length(ttemp)>10
%         wave_tinterp(n,:) = interp1(xtemp(~isnan(xtemp)),ttemp(~isnan(xtemp)),Xgrid(1:xlimit),'nearest');
%
%         n = n+1;
%     end
% end
%
%
%
% %%
% % align to max starting x location
% maxn = max(wave_x(:,1));
%
% for i=1:nwave
%     indAlign(i) = knnsearch(wave_x(i,:)',maxn);
% end
%
% %adjust for new times added
% wave_xx = [nan(nwave,max(indAlign)) nan(size(wave_x))];
% wave_hh = wave_xx;
%
% for n=1:nwave
%     wave_xx(n,max(indAlign)-indAlign(n)+1:max(indAlign)-indAlign(n)+ntime) = wave_x(n,:);
%     wave_hh(n,max(indAlign)-indAlign(n)+1:max(indAlign)-indAlign(n)+ntime) = wave_h(n,:);
%
% end
%
% timeAdj = ones(nwave,1)*(max(indAlign):-1:1);
% wave_tt = [wave_t(:,1)-timeAdj wave_t];
% wave_tt = wave_tt-wave_tt(:,1);
%
%
% %% make plot of wave crests
% figure
% tchunk = 1000:3000;
%
% pcolor(Xgrid(:,1),tchunk,TXdrone2(:,tchunk)'); shading flat
% hold on
% % scatter(wave_x(:),wave_t(:),'or','filled')
% for i=1:nwave
%     scatter(wave_x(i,:),wave_t(i,:),'o','filled')
%     hold on
%     %
% end
% ylim([tchunk(1) tchunk(end)])
%
% xlabel('x(m)')
% ylabel('t(0.1s)')
% hc = colorbar;
% ylabel(hc,'z (m)')
% caxis([1.5 4])
% % print(gcf,'-dpng',['../viz/' hoverdate '_H' num2str(hovern) '_lidarBathy_cresttrack.png'])
%
%
% %% get moving slope
% Cp = nan(nwave,xlimit);
% Cp_H = nan(nwave,xlimit);
%
% for i =1:nwave
%     tt = 1:346;
%     wavex = wave_xx(i,~isnan(wave_xx(i,:)));
%     waveh = wave_hh(i,~isnan(wave_xx(i,:)));
%
%     wavet = tt(~isnan(wave_xx(i,:)));
%     [a,b] = unique(wavex);
%
%     wx = interp1(a,wavet(b),Xgrid(1:xlimit));
%     wh = interp1(a,waveh(b),Xgrid(1:xlimit));
%     plot(Xgrid(1:xlimit),wx,'o')
%     %
%
%     chunk = 100;
%     n = 1;
%     while n <xlimit-chunk
%         nn = n:n+chunk-1;
%
%
%         xloc = Xgrid(nn);
%         tloc = wx(nn);
%         xloc = xloc(~isnan(tloc))';
%         tloc = tloc(~isnan(tloc))';
%
%         E = [ones(size(tloc)) tloc];
%         c = E\xloc;
%         Cp(i,n) = Hz_lidar*c(2); % wave speed is the slope
%         Cp_H(i,n) = range(wh(nn)); % get max wave height over range of subwindow (COULD BE BETTER)
%         n = n+1;
%     end
% end
% %% get wave speed from moving slope
% Cp(Cp<=0) = nan;
% Cp_H(isnan(Cp)) = nan;
% Cpm = nanmedian(Cp);
%
% h = Cpm.^2./9.81;
% h2 = nanmedian(Cp.^2./9.81-Cp_H./2);
%
%
%% run the other two codes to get a plot with all the things
phase_speed_gradient
compare_jumbo_h_est
end
%
% %% check for wave height dependence on speed
% Cp_x = Xgrid(1:xlimit,1:nwave)';
% Cp_x(isnan(Cp)) = nan;
% figure
% hs = scatter(Cp_x(:),Cp(:),50',Cp_H(:),'+');
% hs.MarkerEdgeAlpha = 0.1;
% hc = colorbar;
% ylabel('Wave speed C (m/s)')
% ylabel(hc,'Wave Height H (m)')
% xlabel('X (m)')
% caxis([0 0.5])
% % xlim([0 10])
% ylim([ 0 10])
% % title('Matlab representation of toddler art')
% % print(gcf,'-dpng',['../viz/' hoverdate '_H' num2str(hovern) '_lidarBathy_Hdependence.png'])
%
%
% %% JUNK PILE BELOW
% %% OLD STUFF, maybe can be revived
% % track the peak of the waves over 5m x 10 sec
% % dx = 5*10; %50
% % dt = 10*10;
% %
% % nt = floor(length(tvecHover)./dt);
% % tstart = 1;
% % for tvec = 1:10000
% %
% %     tchunk = tstart:tstart+dt;
% % n = TXdrone2(100:200,tchunk);
% % crest = edge(n,'Sobel');
% % % imagesc(crest)
% %  theta = 0:90;
% %
% % [R,xp] = radon(crest,theta);
% % maxTheta = theta(var(R) == max(var(R)));
% % slope = tand(maxTheta);
% % c_est(tvec) = 1/slope;
% % tstart = tstart+1;
% % end
%
% % imshow(edges(:,1:1400),[])
% %%
% % % do radon transformation
% % n = TXdrone2(1:200,:);
% % [R,xp] = radon(TXdrone2(1:200,:),theta);
% % clear R
% % for i=1:10
% % [R(i).r,xp] = radon(edges(1+50*(i-1):50*i,:),theta);
% % plot(theta,var(R(i).r))
% % hold on
% % end
%
%
% %% make a movie
% % % etamt = etam;
% % close all
% % figure
% % scatter(X1211,z1211,'.')
% % hold on
% % scatter(X1217,z1217,'.')
% % ylim([-1 3])
% % xlim([-100 10])
% % kk = Xgrid(:,1)>-20;
% % plot(Xgrid(kk),nanmin(TXdrone2(kk,:)'),'linewidth',2)
% % kk = Xgrid(:,1)>-150;
% % plot(Xgrid(kk),nanmedian(TXdrone2(kk,:)'),'linewidth',2)
% % plot(Xgrid(kk),nanmean(TXdrone2(kk,:)'),'linewidth',2)
% %
% %
% % for i=1000:2100
% %     hl1 = plot(Xgrid(1:xlimit),TXdrone2(1:xlimit,i),'color',[0.5 0.5 0.5 0.3]);
% %     hold on
% %     %     hl2 = plot(Xgrid(1:xlimit),TXdrone(1:xlimit,i),'.k');
% %     hl4 = plot(Xgrid(1:xlimit),TXtruck2(1:xlimit,i),'c');
% %     hl2 = plot(Xgrid(1:xlimit),etam(1:xlimit,i)+nanmean(TXdrone2(1:xlimit,:),2),'linewidth',2);
% %
% %     hl3 = plot(Xgrid(1:xlimit),etamt(1:xlimit,i)+nanmean(TXdrone2(1:xlimit,:),2),'r');
% %
% %
% %         pause(0.1)
% %         delete(hl1)
% %         delete(hl2)
% %         delete(hl3)
% %         delete(hl4)
% % end
