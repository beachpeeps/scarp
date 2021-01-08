%plot_HsDroneTruck
% script to compare the drone with the truck lidar
% will spit out 2 figures: pcolor of coherence between truck and drone, and
% spectra @ a paros location
%choose hover number and paros number and load data
%choose hover number  and load data
clear all
hovern = 3;
datadir = '/Volumes/FiedlerBot8000/scarp/';

hoverdate = '20191214';
% hoverdate = '20200224';

% switch the limits for analysis, depending on hoverdate, to stop at x=0
if str2double(hoverdate) == 20191214
    xlimit = 1000;
    hovers = [1 3 4];
elseif str2double(hoverdate) == 20200224
    xlimit = 500;
    hovers = 2:5;
end

   
tstack  = load([datadir '/mat/timestacks/' hoverdate '_' num2str(hovern) '.mat']);
load('../mat/paros.mat') 


drone = tstack.TXdrone2;
truck = tstack.TXtruck2;
t = datenum(tstack.tvecHover);
x = tstack.Xgrid(1:xlimit,1);


%%
pp = 12;
dt = 0.1;
[PSI,~] = sleptap(length(t),pp);

for i=1:xlimit
    [f,Sdd(i,:),Stt(i,:),Sdt(i,:)]=mspec(dt,detrend(drone(i,:))',detrend(truck(i,:))',PSI);
end

fSS = f>0.04 & f<=0.25;
fIG = f>0.004 & f<=0.04;
df = f(2)-f(1);

Hss_truck = 4*sqrt(sum(Stt(:,fSS),2).*df);
Hig_truck = 4*sqrt(sum(Stt(:,fIG),2).*df);
Hss_drone = 4*sqrt(sum(Sdd(:,fSS),2).*df);
Hig_drone = 4*sqrt(sum(Sdd(:,fIG),2).*df);
Htot_drone = 4*sqrt(sum(Sdd,2).*df);

%%
for parosn = 3:5
% find which hour to pull from
tstarthour = dateshift(tstack.tvecHover(1),'start','hour');
tstophour = dateshift(tstack.tvecHover(end),'start','hour');
tvecHover2Hz = tstack.tvecHover(1):seconds(1/2):tstack.tvecHover(end);

tind = find(dateshift(paros(parosn).t+minutes(1),'start','hour') == tstarthour);

if tstophour ~= tstarthour
    tind = [tind tind+1];
end


% make a time vector for the paros hourly chunk
for i=1:length(tind)
    tvec(:,i) = paros(parosn).t(tind(i)):seconds(0.5):paros(parosn).t(tind(i))+7166*seconds(0.5);
end
tvec = tvec(:);

if ~isempty(paros(parosn).offset)
    pvec = paros(parosn).etaCorrected(:,tind(:))+paros(parosn).z-paros(parosn).offset;
else
    pvec = paros(parosn).etaCorrected(:,tind(:))+paros(parosn).z;
end

TSparos = interp1(tvec,pvec(:),tvecHover2Hz);

dt = 0.5;
[PSI,~] = sleptap(length(tvecHover2Hz),pp);
[f,Spp]=mspec(dt,detrend(TSparos)',PSI);
fSS = f>0.04 & f<=0.25;
fIG = f>0.004 & f<=0.04;
df = f(2)-f(1);
Hss_paros(parosn) = 4*sqrt(sum(Spp(fSS)).*df);
Hig_paros(parosn) = 4*sqrt(sum(Spp(fIG)).*df);
Htot_paros(parosn) =  4*sqrt(sum(Spp).*df);
xp(parosn) = paros(parosn).crossshore;
clear tvec
end

xp(1:2) = [];


%%
close all
figwidth = 8;
figheight = 4;
nrow = 1;
ncol = 1;
options.margin = 0.03;
options.units = 'centimeters';
options.axesfontsize = 9;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);

colors = lines(3);
axes(ax(1))
plot(x,Hss_drone,'color',colors (1,:))
hold on
plot(x,Hig_drone,'--','color',colors (1,:))
plot(x,Hss_truck,'color',colors (2,:))
plot(x,Hig_truck,'--','color',colors (2,:))
scatter(xp,Hss_paros(3:5),50,colors(3,:),'o','filled','markeredgecolor','k')
scatter(xp,Hig_paros(3:5),100,colors(3,:),'p','filled','markeredgecolor','k')
% hl = legend('Hss drone','Hig drone', 'Hss truck', 'Hig truck');
% hl.NumColumns = 2;
ylim([0 1.5])
ylabel('Hs (m)')
xlabel('cross-shore (m)')
ax(1).Position(2) = 0.15;
ax(1).Position(1) = 0.15;
ax(1).Position(3) = 0.8;

print(gcf, '-dpdf', ['../viz/paperFigures/Hsig.pdf'],'-r300');
