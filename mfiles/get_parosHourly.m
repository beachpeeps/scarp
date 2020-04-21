%plot_parosHourly
% NEEDS TO BE CLEANED UP, NOT FOR PRIMETIME
% April 21,2020 Julia Fiedler jfiedler@ucsd.edu

clear
load('/Users/juliafiedler/Documents/Repositories/scarp/mat/sensors.mat')
parosFolder = '~/Documents/SCARP/data/processed/';
parosName = dir([parosFolder 'P*']);
addpath /Users/juliafiedler/Documents/SWASH_runup/mfiles/functions/


for pnum = 1:10
paros(pnum).name = parosName(pnum).name;

%Make all on the same hourly time
fileList = dir([parosFolder parosName(pnum).name '/2*']);

minTime = datetime(fileList(1).name(1:12),'InputFormat','yyyyMMddHHmm');
maxTime = datetime(fileList(end).name(1:12),'InputFormat','yyyyMMddHHmm');

% hourVec = minTime:hours(1):maxTime;


for nfile = 1:length(fileList)
    load([parosFolder parosName(pnum).name '/' fileList(nfile).name],'h','pSensorTime')
    
    ptemp(nfile) = mean(h);
    Htemp(nfile) = 4*nanstd(h);
    ttemp(nfile) = dateshift(datetime(pSensorTime(20),'ConvertFrom','datenum'),'start','hour');
    etatemp(:,nfile) = h;
end

[tsort,indsort] = sort(ttemp);

paros(pnum).h = ptemp(indsort);
paros(pnum).t = tsort;
paros(pnum).Hs = Htemp(indsort);
paros(pnum).eta = etatemp(:,indsort);

paros(pnum).latitude = P.DeploymentLatitude(pnum);
paros(pnum).longitude = P.DeploymentLongitude(pnum);
paros(pnum).z = P.DeploymentNAVD88_m_Geoid12B_(pnum);
paros(pnum).UTMx = P.UTMEastings_Zone11_(pnum);
paros(pnum).UTMy = P.UTMNorthings_Zone11_(pnum);

clear hourVec h ptemp ttemp t Htemp etatemp 
end
%% get tide
[t_verified,t_predicted,verified,predicted,tideInfo] = getNOAAtide(minTime,maxTime,'9410230');
% [t,verified,tideInfo] = getNOAAtide_hourly(datestr(minTime,'yyyymmdd'),datestr(maxTime,'yyyymmdd'),'9410230');
t = datetime(t_verified,'ConvertFrom','datenum');
%% plot timeseries

savedir = '../viz/';
savename = 'hourlyparos';
figure
ax(1) = subplot(2,1,1)
for pnum=1:10
    plot(paros(pnum).t,paros(pnum).h+paros(pnum).z)
    hold on
end
plot(t,verified,'k')
xlim([datetime(2019,11,20) datetime(2020,3,15)])
ylabel('z NAVD (m)')
title('Pressure Sensors 1:10')
legend(P.Properties.RowNames)
ax(2) = subplot(2,1,2)
plot(paros(6).t,paros(6).Hs)
xlim([datetime(2019,11,20) datetime(2020,3,15)])
ylim([0 1.5])
linkaxes(ax,'x')
print(gcf, '-djpeg', [savedir savename]);
%% plot xshore array

% get all on the same time vector
minTime = datetime(2019,11,20);
maxTime = datetime(2020,3,15);
dateVec = minTime+minutes(10):minutes(20):maxTime-minutes(10);
Pnum = nan(10,length(dateVec));

%TODO: fix paros 1, what's going on at index at new year?? ind = 1108/1107
for pnum=2:10
    Pmat(pnum,:) = interp1(paros(pnum).t,paros(pnum).h,dateVec);
end

% tide = interp1(t,verified,dateVec);
%%

X = P.UTMEastings_Zone11_;
Y = P.UTMNorthings_Zone11_;
THETA = deg2rad(theta);
XO = P.UTMEastings_Zone11_(1);
YO = P.UTMNorthings_Zone11_(1);

[XR YR] = xyRotate(X,Y,THETA,XO,YO);
%%
hFig = figure

tideHour = interp1(t,verified,dateVec);
igif = 1;
for indHour = 100:500
    clf
    for pnum = 2:10
        scatter(XR(pnum),Pmat(pnum,indHour)+paros(pnum).z)
        hold on
        text(XR(pnum),Pmat(pnum,indHour)+paros(pnum).z,P.Properties.RowNames(pnum))
        plot([-130 20],[1 1]*tideHour(indHour),':k')
    end
    ylim([-0.5 2.5])
    title(datestr(dateVec(indHour)))
    xlabel('X-shore location (m)')
    
    print(['~/Desktop/temp/Frame' num2str(igif,'%02.0f')], '-dpng','-r150');
    drawnow
    F = getframe(hFig);
    [im{igif},map] = frame2im(F);
    igif = igif+1;
end

%%
gifname = [savedir 'paros.gif'];
for idx = 1:igif-1
    [A,~] = imread(['~/Desktop/temp/Frame' num2str(idx,'%02.0f') '.png']);
    [X,map] = rgb2ind(A,256);
    if idx == 1
        imwrite(X,map,gifname,'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(X,map,gifname,'gif','WriteMode','append','DelayTime',0);
    end
end

%%
% ldir = '/Users/juliafiedler/Documents/Repositories/scarp/data/las/drone/';
% lname = '20191214_00582_TorreyRunup_H3.las';
%
figure

load('/Users/juliafiedler/Documents/Repositories/scarp/mat/lidar/drone/20191209_00582_TorreyRunup_H3_05m.mat')
plot(Processed.x,nanmean(Processed.Zinterp2)-0.05)
hold on
plot(Processed.x,nanmin(Processed.Zinterp2)-0.05)
plot(Processed.x,nanmax(Processed.Zinterp2)-0.05)
plot(Processed.x,mode(Processed.Zinterp2)-0.05)

% 
% c = lasdata([ldir lname]);
% tSecFromSunday = get_gps_time(c);
% 
% [tUTC,tLST] = GPStoUTC(tSecFromSunday,lname);
% [T, x, y, z, amp, sortedInd] = sortLASobject(tUTC,c);
% % rotate data
% [XRl, YRl] = xyRotate(x,y,THETA,XO,YO);
% %
% ind = 1:1000:length(z);
% hold on
% h = scatter(XRl(ind),z(ind));

%
minTime = datenum(dateshift(datetime(nanmin(Processed.t),'ConvertFrom','datenum'),'start','minute'));
maxTime = datenum(dateshift(datetime(nanmax(Processed.t),'ConvertFrom','datenum'),'end','minute'));

indTide = find(t_verified>=minTime & t_verified<=maxTime);

tideMean = mean(verified(indTide));

ind = find(dateVec == dateshift(datetime(nanmean(Processed.t),'ConvertFrom','datenum'),'start','minute'));

 for pnum = 2:10
        scatter(XR(pnum),Pmat(pnum,ind)+paros(pnum).z)
        hold on
        text(XR(pnum),Pmat(pnum,ind)+paros(pnum).z,P.Properties.RowNames(pnum))
        plot([-130 20],[1 1]*tideMean,':k')
    end
    ylim([0.5 2.5])
    title(datestr(dateVec(ind)))
    xlabel('X-shore location (m)')




% print(hFig, '-dpdf', [savedir savename]);




