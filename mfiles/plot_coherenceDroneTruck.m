%plot_coherenceDroneTruck
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
clear f Axy COH2xy PHIxy CI_coh CI_phi
pp = 12;
dt = 0.1;
[PSI,~] = sleptap(length(t),pp);

for i=1:xlimit
    [f,Sdd(i,:),Stt(i,:),Sdt(i,:)]=mspec(dt,detrend(drone(i,:))',detrend(truck(i,:))',PSI);
end


%% Calculate Amplitude, Phase, and Coherence
Axy = abs(Sdt); % cross amplitude spectrum
Lxy = real(Sdt); % co spectrum, associated with cosine transform
Qxy = imag(Sdt); % quadrature spectrum, associated with sine transform

PHIxy = atan2(-Qxy,Lxy); % phase spectrum
PHIxy = rad2deg(PHIxy);

COHxy = Axy./sqrt(Sdd.*Stt); % coherence
COH2xy = Axy.^2./(Sdd.*Stt); % squared coherence
%% CONFIDENCE INTERVALS and Standard Devs
df = f(2)-f(1); %frequency resolution
alpha = 0.05; % confidence level
nens = pp-1;
dof = 2.4*nens; % 2.4 because of the hanning window (maybe?) FIX THIS
stdA = (1/sqrt(dof))*Axy.*(1+(COHxy.^-1)).^0.5; % 1 std of A

% Confidence Level (95% Significance) for Coherence
CI_coh = 1-alpha^(2/(dof-2));
fSS = f>0.04 & f<=0.25;
fIG = f>0.004 & f<=0.04;
Hss_truck = 4*sqrt(sum(Stt(:,fSS),2).*df);
Hig_truck = 4*sqrt(sum(Stt(:,fIG),2).*df);
Hss_drone = 4*sqrt(sum(Sdd(:,fSS),2).*df);
Hig_drone = 4*sqrt(sum(Sdd(:,fIG),2).*df);
Htot_drone = 4*sqrt(sum(Sdd,2).*df);


%%
close all
figwidth = 8;
figheight = 16;
nrow = 5;
ncol = 1;
% options.margin = 0.05;
options.units = 'centimeters';
options.axesfontsize = 9;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);
for i=1:5
ax(i).Position(1) = 0.14;
end

axes(ax(1))
ind = find(COH2xy < CI_coh); % find Coherence less than signifance levels
fSig = f;
COHsig = COH2xy; % rename to establish matrix of significant coherence
COHsig(ind) = NaN; %nan out everything that isn't coherent

hold on
hp = pcolor(x,f,COHsig'); shading flat
hc = colorbar(ax(1));
hc.Label.String = {'Coherence squared'};
hc.Location = 'NorthOutside';
hlim = hc.Limits;

% delete(hp)
delete(hc)
% hp2 = pcolor(x,f,nan(size(COHsig'))); shading flat
hc = colorbar(ax(1));
hc.Label.String = {'Coherence squared'};
hc.Location = 'NorthOutside';
hc.Limits = hlim;
% 
ax(1).YScale = 'log';
ylabel('Frequency (Hz)')


% plot(ax(1).XLim,[10^log10(0.04) 10^log10(0.04)],'-.k','linewidth',1)
% plot(ax(1).XLim,[10^log10(0.004) 10^log10(0.004)],'-.k','linewidth',1)
% plot(ax(1).XLim,[0.25 0.25],'-.k', 'linewidth',1)
% 
% text(2,10^-2,'IG','fontsize',8,'fontweight','bold')
% text(2,10^-1.1,'SS','fontsize',8,'fontweight','bold')
ax(1).Color = [0.85 0.85 0.85];

% box on 

  %
% for i=2:5
%     hs(i) = scatter(paros(i).crossshore,5,'vk','filled');
%     hs(i).SizeData = 50;
%     hs(i).MarkerFaceColor = [1 1 1];
%     hs(i).MarkerEdgeColor = [0 0 0];
%         hst(i) = text(paros(i).crossshore,6.5,['P' num2str(i)],...
%             'fontsize',8,'verticalalignment','bottom',...
%             'horizontalalignment','center');
% 
% end
ax(1).YLim(1) = 1e-3;













% check individual locations
p = 3;
for parosn = 3:5
    axes(ax(p))


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
pvec = pvec(:);

xind = knnsearch(x(:),paros(parosn).crossshore);
plot(tstack.tvecHover,drone(xind,:),'linewidth',0.5)
hold on
plot(tstack.tvecHover,truck(xind,:),'linewidth',0.5)
plot(tvec,pvec,'linewidth',0.5)
xlim([tstack.tvecHover(1) tstack.tvecHover(1)+minutes(4)])

TSparos = interp1(tvec,pvec(:),tvecHover2Hz);

dt = 0.5;
[PSI,~] = sleptap(length(tvecHover2Hz),pp);
[fp,Spp]=mspec(dt,detrend(TSparos)',PSI);
df = fp(2)-fp(1);
Hss_paros(parosn) = 4*sqrt(sum(Spp(fSS)).*df);
Hig_paros(parosn) = 4*sqrt(sum(Spp(fIG)).*df);
Htot_paros(parosn) =  4*sqrt(sum(Spp).*df);
xp(parosn) = paros(parosn).crossshore;




p = p+1;
clear tvec pvec
end

axes(ax(2))

cmap = lines(3);
dronecolor = cmap(1,:);
truckcolor = cmap(2,:);
paroscolor = cmap(3,:);

hd1 = plot(x,Hss_drone,'color',dronecolor);
hold on
plot(x,Hig_drone,'--','color',dronecolor)
ht1 = plot(x,Hss_truck,'color',truckcolor);
plot(x,Hig_truck,'--','color',truckcolor)
hp1 = scatter(xp(3:5),Hss_paros(3:5),50,paroscolor,'o','filled','markeredgecolor','k');
scatter(xp(3:5),Hig_paros(3:5),100,paroscolor,'p','filled','markeredgecolor','k')

%paros
ax(2).YLim = [0 1.5];
  %
for i=3:5
    hs(i) = scatter(paros(i).crossshore,1.5,'vk','filled');
    hs(i).SizeData = 50;
    hs(i).MarkerFaceColor = [1 1 1];
    hs(i).MarkerEdgeColor = [0 0 0];
end
xlabel('cross-shore (m)')


hl = legend([hd1 ht1 hp1],'drone','truck','paros');
hl.NumColumns = 2;
% hl.Position = [0.45    0.51    0.2225    0.0097]
hl.Box = 'off';


for i=3:5
    if i<5
    ax(i).XLabel = [];
    ax(i).XTickLabel = [];
    end
    ax(i).Position(4) = 0.09;
end

for i=1:5
    ax(i).Position(3) = 0.8;
end

ax(4).YLabel.String = 'z (m NAVD88)';
ax(2).YLabel.String = 'Hs (m)';

ax(5).Position(2) = 0.08;
%move stuff
    ax(4).Position(2) = ax(5).Position(2) + ax(5).Position(4)+0.01;
    ax(3).Position(2) = ax(4).Position(2) + ax(4).Position(4)+0.01;
    ax(2).Position(2) = ax(3).Position(2) + ax(3).Position(4)+0.07;

% ax(1).Position(3) = 0.8;
%
% ss = {' a ',' b ',' c ',' d ',' e '};
% ii = 1
% for i=[1:5]
% text(0.02,0.83,ss{i},'parent',ax(i),'units','normalized','fontsize',10,'EdgeColor','k','backgroundColor','white','margin',0.5)
% ax(i).FontSize = 8;
% 
% end

% ax(2).Visible = 'off';
ax(1).Position(4) = 0.18;
ax(2).Position(4) = 0.15;

ax(1).Position(2) = ax(2).Position(2)+ ax(2).Position(4)+0.03;
ax(1).Position(4) = 0.23;
ax(1).XTickLabel = [];

hc.Position(2) = 0.9;
% hc.Position
% print(gcf, '-dpdf', ['../viz/paperFigures/coherence3.pdf'],'-r300', '-painters');
print(gcf, '-djpeg', ['../viz/paperFigures/coherence.jpeg'],'-r300');
