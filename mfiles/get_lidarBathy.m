% track wave crests
clearvars -except fracNaN_drone
hoverdate = '20191214';
hoverdate = '20200224';

% switch the limits for analysis, depending on hoverdate, to stop at x=0
if str2double(hoverdate) == 20191214
    xlimit = 1000;
    hovers = [1 3 4];
elseif str2double(hoverdate) == 20200224
    xlimit = 500;
    hovers = 2:5;
end
% for hovern = hovers
for hovern = 3;
    clearvars -except hovern fracNaN hovers hoverdate xlimit
    
    
    datadir = '/Volumes/FiedlerBot8000/scarp/';
    
    tstack  = load([datadir '/mat/timestacks/' hoverdate '_' num2str(hovern) '.mat']);
    
    %% decide if you want to look at truck or drone data
    eta =tstack.TXdrone2-nanmean(tstack.TXdrone2,2); % interpolated
    eta2 = tstack.TXdrone-nanmean(tstack.TXdrone2,2); % not interpolated
    
%     eta = tstack.TXtruck2-nanmean(tstack.TXtruck2,2);
%     eta2 = tstack.TXtruck-nanmean(tstack.TXtruck2,2);
    
    [nx,nt] = size(eta);
    x = tstack.Xgrid(:,1);
    dx =x(2)-x(1);
    
    % get fraction of NaNs in un-interpolated gridded data
    fracNaN(hovern,:) = sum(isnan(eta2'))./nt;
    
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
    xloc_startcrest = round(1/dx); % 1m from edge of domain
    
    % get timechunk for finding peaks at xloc_startcrest
    window_t = floor(tstack.Hz_lidar/fp(xloc_startcrest)/2);
    
    % get elevation data at edge of domain
    z = eta(xloc_startcrest,:);
    
    % re-organize elevation data
    length_z = nt - mod(nt,window_t);
    zz = reshape(z(1:length_z),window_t,length_z/window_t);
    
    % find max points for each timechunk
    indPeak = find(zz == max(zz));
    
    % nudge the indPeak back a bit in time to account for xloc_startcrest
    indPeak(indPeak<10) = [];
    
    
    % make sure the peak is actually a peak
    indPeak(z(indPeak)<0) = [];
    isPeak = (z(indPeak) - z(indPeak-2)).*(z(indPeak) - z(indPeak+2));
    indPeak(isPeak<0) = [];
    
    %  quick plot to see if method works
    plot(z,'.')
    plot(z)
    hold on
    plot(indPeak,z(indPeak),'or')
    title('check peak-finding at edge of domain')
    
    
    %% crest-tracking.
    % For each wave crest identified at outer edge of domain, track the crest
    % as it comes onshore, by time-stepping and looking for the peak location
    % in a sliding cross-shore window around the previous crest location.
    
    
    figure
    
    window_length = round(5/dx); % window size in meters, used by Martins et al.
    max_time = 500; % maximum time value to track wave (arbitrary high number)
    nwave = length(indPeak);
    
    for nw =1:nwave-2
        
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
            
            % if there's no forward progress of the search window in x
            if tti>5 && mean(xx(tti-5:tti-1))>=x(xi)
                
                % account for a bore, get most shoreward point of wave crest
                xi = window_x(find(eta(window_x,ti)>0.75*max(eta(window_x,ti)),1,'last'));
                
                % stop crest tracking if the surface is flat
                % or if there is no forward progress
                % or if there is too much forward progress
                if isempty(xi) || tti>5 && mean(xx(tti-5:tti-1))>=x(xi) || x(xi)>-10 && (x(xi)-xx(tti-1))>1.5
                    disp('hi')
                    break
                else
                end
                
            end
            
            % stop if too much forward progress
            if x(xi)>-10 && (x(xi)-xx(tti-1))>1.5
                break
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
            xsearch = xi:xi+window_length; %only search in front of the wave crest
            htrough(tti) = min(eta(xsearch,ti)); %minimum value is the trough height
            
            % get the most shoreward location of trough value (good for
            % plotting only, not used in any calculation.)
            xtrough(tti) = max(x(xsearch(eta(xsearch,ti) == min(eta(xsearch,ti)))));
            xtroughInd(tti) = find(x==xtrough(tti));
            
            
            H(tti) = hpeak(tti)-htrough(tti); %wave height is crest-trough
            
            % disregard spurious waves less than 10cm
            if tti>3 && mean(H(tti-3:tti))<0.1
                break
            end
            
            
            % plot a sample timestack of the waves and the roller faces (only use every 1st timestep)
            tstep = 1;
            if mod(tti,tstep) == 0 && nw == 11 %use the 40th wave (just to show how it works!)
                h1 = plot(x(1:xlimit),eta(1:xlimit,ti)+tti/tstep,'b'); % wave
                hold on
                h3(1) = plot(x(xi),hpeak(tti)+tti/tstep,'ok'); %roller face
                h1(2) = plot(x(xi:xtroughInd(tti)),eta(xi:xtroughInd(tti),ti)+tti/tstep,'c','linewidth',2);
                h3(2) = plot(xtrough(tti),htrough(tti)+tti/tstep,'og'); %roller face
                
                wave(1).xtrough(tti) = xtrough(tti);
                wave(1).hpeak(tti) = hpeak(tti);
                wave(1).htrough(tti) = htrough(tti);
                wave(1).xpeak(tti) = x(xi);
                wave(1).t(tti) = indPeak(nw)+tti;
                wave(1).nwave = nw;
                
            end
            
             if mod(tti,tstep) == 0 && nw == 10 %use the preceeding wave for plotting
                h1 = plot(x(1:xlimit),eta(1:xlimit,ti)+tti/tstep,'b'); % wave
                hold on
                h3(1) = plot(x(xi),hpeak(tti)+tti/tstep,'ok'); %roller face
                h1(2) = plot(x(xi:xtroughInd(tti)),eta(xi:xtroughInd(tti),ti)+tti/tstep,'c','linewidth',2);
                h3(2) = plot(xtrough(tti),htrough(tti)+tti/tstep,'og'); %roller face
                
                wave(2).xtrough(tti) = xtrough(tti);
                wave(2).hpeak(tti) = hpeak(tti);
                wave(2).htrough(tti) = htrough(tti);
                wave(2).xpeak(tti) = x(xi);
                wave(2).t(tti) = indPeak(nw)+tti;
                wave(2).nwave = nw;
                
            end
            
            % move to next timestep, establish new cross-shore window
            tti = tti+1;
            if xi<26
                window_x = xi:xi+window_length/2;
            else
                window_x = xi-window_length/2:xi+window_length/2;
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
    
    %% disregard waves without enough data points (COULD BE FIXED)
    minSec = 5/dt; % minimum number of seconds we need for crest-tracking
    lwave = sum(~isnan(wave_x),2);
    wave_x(lwave<minSec,:) = [];
    wave_h(lwave<minSec,:) = [];
    wave_t(lwave<minSec,:) = [];
    [nwave,ntime] = size(wave_x);
    %% plot wave trajectories on pcolor plot
    figure
    tchunk = 500:2500;
    pcolor(x(1:xlimit),dt*(tchunk),eta2(1:xlimit,tchunk)')
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
    window_length = 100;
    window_length_adapt = -x*1.5;
% figure
% ax1 = subplot(2,1,1);
% ax2 = subplot(2,1,2);
    
    for i =1:nwave
        
        wavex = wave_x(i,~isnan(wave_x(i,:)));
        waveh = wave_h(i,~isnan(wave_x(i,:)));

        wavet = 1:length(wavex);
        [a,b] = unique(wavex);
        %
        wt = dt.*interp1(a,wavet(b),x(1:xlimit));
        wh = interp1(a,waveh(b),x(1:xlimit));
        
%         p = polyfit(wavet(b)*dt,a,3);
%         polyFit = polyval(p,wt);
%         dp = p.*[ 3 2 1 0];
% %         Cp(i,:) = dp(1)*wt+dp(2);
%         Cp(i,:) = dp(1)*wt.^2+dp(2)*wt+dp(3);

        Cp_H(i,:) = wh;
        %

%         
        nwaveplot = 20; % pick which wave to make a plot of
        n = round(window_length_adapt(1)/2)+1;
        while n < min(find(x == max(wavex)),xlimit-window_length_adapt(1)/2-1)
            
            nn = n-window_length/2:n+window_length/2; %<-- this is a XXm window, and doesn't change. Issue?
            nn = round(n-round(window_length_adapt(n)/2):n+round(window_length_adapt(n)/2)); %<--changing XXm window with x locaiton

            
            % get subset of x,t points in moving window
            xloc = x(nn);
            tloc = wt(nn);
            
            %analyze only chunk of wave with data
            xloc = xloc(~isnan(tloc));
            tloc = tloc(~isnan(tloc));
            
            %make sure there are enough points to regress to
            if numel(xloc)<20 %this is arbitrary
                break
            end
            
            % simple linear regression to get slope
%             E = [ones(size(tloc)) tloc];
%             c = E\xloc;

%                     p = polyfit(tloc,xloc,1);
                    b = robustfit(tloc,xloc);
                    warning('off','stats:statrobustfit:IterationLimit');
% [p2,stats] = polyfit(tloc,xloc,2);
% dp2 = p2.*[2 1 0];
% cpFit = dp(1)*tloc+dp(2);
%         polyFit = polyval(p,wt);
%         dp = p.*[ 1 0];
%         Cp(i,n) = dp(1);
        Cp(i,n) = b(2);
%         Cp(i,n) = dp(1)*nanmean(wt(nn))+dp(2);
            
            
%             Cp(i,n) = tstack.Hz_lidar*c(2); % wave speed is the slope
            Cp_H(i,n) = nanmedian(wh(nn)); % get median wave height of subwindow
            
            %mid loop diagnostic plot for troubleshooting porpoises
%          
    % nudge the indPeak back a bit in time to account for xloc_startcrest
    indPeak(indPeak<10) = [];
    
    
    % make sure the peak is actually a peak
    indPeak(z(indPeak)<0) = [];
    isPeak = (z(indPeak) - z(indPeak-2)).*(z(indPeak) - z(indPeak+2));
    indPeak(isPeak<0) = [];
    
    %  quick plot to see if method works
    plot(z,'.')
    plot(z)
    hold on
    plot(indPeak,z(indPeak),'or')
    title('check peak-finding at edge of domain')
    
    
    %% crest-tracking.
    % For each wave crest identified at outer edge of domain, track the crest
    % as it comes onshore, by time-stepping and looking for the peak location
    % in a sliding cross-shore window around the previous crest location.
    
    
    figure
    
    window_length = round(5/dx); % window size in meters, used by Martins et al.
    max_time = 500; % maximum time value to track wave (arbitrary high number)
    nwave = length(indPeak);
    
    for nw =1:nwave-2
        
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
            
            % if there's no forward progress of the search window in x
            if tti>5 && mean(xx(tti-5:tti-1))>=x(xi)
                
                % account for a bore, get most shoreward point of wave crest
                xi = window_x(find(eta(window_x,ti)>0.75*max(eta(window_x,ti)),1,'last'));
                
                % stop crest tracking if the surface is flat
                % or if there is no forward progress
                % or if there is too much forward progress
                if isempty(xi) || tti>5 && mean(xx(tti-5:tti-1))>=x(xi) || x(xi)>-10 && (x(xi)-xx(tti-1))>1.5
                    disp('hi')
                    break
                else
                end
                
            end
            
            % stop if too much forward progress
            if x(xi)>-10 && (x(xi)-xx(tti-1))>1.5
                break
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
            xsearch = xi:xi+window_length; %only search in front of the wave crest
            htrough(tti) = min(eta(xsearch,ti)); %minimum value is the trough height
            
            % get the most shoreward location of trough value (good for
            % plotting only, not used in any calculation.)
            xtrough(tti) = max(x(xsearch(eta(xsearch,ti) == min(eta(xsearch,ti)))));
            xtroughInd(tti) = find(x==xtrough(tti));
            
            
            H(tti) = hpeak(tti)-htrough(tti); %wave height is crest-trough
            
            % disregard spurious waves less than 10cm
            if tti>3 && mean(H(tti-3:tti))<0.1
                break
            end
            
            
            % plot a sample timestack of the waves and the roller faces (only use every 1st timestep)
            tstep = 1;
            if mod(tti,tstep) == 0 && nw == 11 %use the 40th wave (just to show how it works!)
                h1 = plot(x(1:xlimit),eta(1:xlimit,ti)+tti/tstep,'b'); % wave
                hold on
                h3(1) = plot(x(xi),hpeak(tti)+tti/tstep,'ok'); %roller face
                h1(2) = plot(x(xi:xtroughInd(tti)),eta(xi:xtroughInd(tti),ti)+tti/tstep,'c','linewidth',2);
                h3(2) = plot(xtrough(tti),htrough(tti)+tti/tstep,'og'); %roller face
                
                wave(1).xtrough(tti) = xtrough(tti);
                wave(1).hpeak(tti) = hpeak(tti);
                wave(1).htrough(tti) = htrough(tti);
                wave(1).xpeak(tti) = x(xi);
                wave(1).t(tti) = indPeak(nw)+tti;
                wave(1).nwave = nw;
                
            end
            
             if mod(tti,tstep) == 0 && nw == 10 %use the preceeding wave for plotting
                h1 = plot(x(1:xlimit),eta(1:xlimit,ti)+tti/tstep,'b'); % wave
                hold on
                h3(1) = plot(x(xi),hpeak(tti)+tti/tstep,'ok'); %roller face
                h1(2) = plot(x(xi:xtroughInd(tti)),eta(xi:xtroughInd(tti),ti)+tti/tstep,'c','linewidth',2);
                h3(2) = plot(xtrough(tti),htrough(tti)+tti/tstep,'og'); %roller face
                
                wave(2).xtrough(tti) = xtrough(tti);
                wave(2).hpeak(tti) = hpeak(tti);
                wave(2).htrough(tti) = htrough(tti);
                wave(2).xpeak(tti) = x(xi);
                wave(2).t(tti) = indPeak(nw)+tti;
                wave(2).nwave = nw;
                
            end
            
            % move to next timestep, establish new cross-shore window
            tti = tti+1;
            if xi<26
                window_x = xi:xi+window_length/2;
            else
                window_x = xi-window_length/2:xi+window_length/2;
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
    
    %% disregard waves without enough data points (COULD BE FIXED)
    minSec = 5/dt; % minimum number of seconds we need for crest-tracking
    lwave = sum(~isnan(wave_x),2);
    wave_x(lwave<minSec,:) = [];
    wave_h(lwave<minSec,:) = [];
    wave_t(lwave<minSec,:) = [];
    [nwave,ntime] = size(wave_x);
    %% plot wave trajectories on pcolor plot
    figure
    tchunk = 500:2500;
    pcolor(x(1:xlimit),dt*(tchunk),eta2(1:xlimit,tchunk)')
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
    window_length = 100;
    window_length_adapt = -x*1.5;
% figure
% ax1 = subplot(2,1,1);
% ax2 = subplot(2,1,2);
    
    for i =1:nwave
        
        wavex = wave_x(i,~isnan(wave_x(i,:)));
        waveh = wave_h(i,~isnan(wave_x(i,:)));

        wavet = 1:length(wavex);
        [a,b] = unique(wavex);
        %
        wt = dt.*interp1(a,wavet(b),x(1:xlimit));
        wh = interp1(a,waveh(b),x(1:xlimit));
        
%         p = polyfit(wavet(b)*dt,a,3);
%         polyFit = polyval(p,wt);
%         dp = p.*[ 3 2 1 0];
% %         Cp(i,:) = dp(1)*wt+dp(2);
%         Cp(i,:) = dp(1)*wt.^2+dp(2)*wt+dp(3);

        Cp_H(i,:) = wh;
        %

%         
        nwaveplot = 20; % pick which wave to make a plot of
        n = round(window_length_adapt(1)/2)+1;
        while n < min(find(x == max(wavex)),xlimit-window_length_adapt(1)/2-1)
            
            nn = n-window_length/2:n+window_length/2; %<-- this is a XXm window, and doesn't change. Issue?
            nn = round(n-round(window_length_adapt(n)/2):n+round(window_length_adapt(n)/2)); %<--changing XXm window with x locaiton

            
            % get subset of x,t points in moving window
            xloc = x(nn);
            tloc = wt(nn);
            
            %analyze only chunk of wave with data
            xloc = xloc(~isnan(tloc));
            tloc = tloc(~isnan(tloc));
            
            %make sure there are enough points to regress to
            if numel(xloc)<20 %this is arbitrary
                break
            end
            
            % simple linear regression to get slope
%             E = [ones(size(tloc)) tloc];
%             c = E\xloc;

%                     p = polyfit(tloc,xloc,1);
                    b = robustfit(tloc,xloc);
                    warning('off','stats:statrobustfit:IterationLimit');
% [p2,stats] = polyfit(tloc,xloc,2);
% dp2 = p2.*[2 1 0];
% cpFit = dp(1)*tloc+dp(2);
%         polyFit = polyval(p,wt);
%         dp = p.*[ 1 0];
%         Cp(i,n) = dp(1);
        Cp(i,n) = b(2);
%         Cp(i,n) = dp(1)*nanmean(wt(nn))+dp(2);
            
            
%             Cp(i,n) = tstack.Hz_lidar*c(2); % wave speed is the slope
            Cp_H(i,n) = nanmedian(wh(nn)); % get median wave height of subwindow
            
            %mid loop diagnostic plot for troubleshooting porpoises
%                     if i==nwaveplot
%             
%             %             pause(0.01)
%                         delete(h1)
%                         h1 = plot(xloc,tloc,'.r','parent',ax1);
%                         plot(x(n),Cp(i,n),'.r','parent',ax2)
%                         hold(ax2,'on')
%                         plot(x(n),Cp2(i,n),'.k','parent',ax2);
%                 pause
%                     end
            n = n+1; 
            
        end
    end
    
    % get wave speed from moving slope
    Cp(Cp<=0) = nan;
    Cp_H(isnan(Cp)) = nan;
    Cpm = nanmedian(Cp); %average in x
    h_linear = Cpm.^2./9.81;
    h_bore = nanmedian(Cp.^2./9.81-Cp_H./2);
    
    
    %% diagnostic plot of stacked wave trajectories
    plot_diagnostic_get_lidarBathy
   
    %% run the other two codes to get a plot with all the things
    phase_speed_gradient
    compare_jumbo_h_est
end
           if i==nwaveplot
%             
%             %             pause(0.01)
%                         delete(h1)
%                         h1 = plot(xloc,tloc,'.r','parent',ax1);
%                         plot(x(n),Cp(i,n),'.r','parent',ax2)
%                         hold(ax2,'on')
%                         plot(x(n),Cp2(i,n),'.k','parent',ax2);
%                 pause
%                     end
            n = n+1; 
            
        end
    end
    
    % get wave speed from moving slope
    Cp(Cp<=0) = nan;
    Cp_H(isnan(Cp)) = nan;
    Cpm = nanmedian(Cp); %average in x
    h_linear = Cpm.^2./9.81;
    h_bore = nanmedian(Cp.^2./9.81-Cp_H./2);
    
    
    %% diagnostic plot of stacked wave trajectories
    plot_diagnostic_get_lidarBathy
   
    %% run the other two codes to get a plot with all the things
    phase_speed_gradient
    compare_jumbo_h_est
end
