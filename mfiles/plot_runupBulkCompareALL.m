clear

F = load('../mat/runupBulkS06_Feb.mat');
D = load('../mat/runupBulkS06_Dec.mat');

Sig_d = [D.Sig_d F.Sig_d];
Sinc_d = [D.Sinc_d F.Sinc_d];
eta_d = [D.eta_d F.eta_d];

err_drone = [D.err_drone; F.err_drone];

Sig_t = [D.Sig_t F.Sig_t];
Sinc_t = [D.Sinc_t F.Sinc_t];
eta_t = [D.eta_t F.eta_t];
err_truck = [D.err_truck; F.err_truck];

Sinc_S06_d = [D.Sinc_S06_d; F.Sinc_S06_d]';
Sig_S06_d = [D.Sig_S06_d; F.Sig_S06_d]';
eta_S06_d = [D.eta_S06_d; F.eta_S06_d]';

Sinc_S06_t = [D.Sinc_S06_t; F.Sinc_S06_t]';
Sig_S06_t = [D.Sig_S06_t; F.Sig_S06_t]';
eta_S06_t = [D.eta_S06_t; F.eta_S06_t]';

tideMean = [D.tideMean F.tideMean];

% get rid of Dec H2 and Fan H1
hovers = [1 3 4 6:9];
Sig_d = Sig_d(hovers);
Sinc_d = Sinc_d(hovers);
eta_d = eta_d(hovers);
err_drone = err_drone(hovers,:);
Sig_t = Sig_t(hovers);
Sinc_t = Sinc_t(hovers);
eta_t = eta_t(hovers);
err_truck = err_truck(hovers,:);

Sinc_S06_d = Sinc_S06_d(hovers);
Sig_S06_d = Sig_S06_d(hovers);
eta_S06_d = eta_S06_d(hovers);

Sinc_S06_t = Sinc_S06_t(hovers);
Sig_S06_t = Sig_S06_t(hovers);
eta_S06_t = eta_S06_t(hovers);

tideMean = tideMean(hovers);


hovers = 1:length(Sinc_S06_t);


%%
close all
figwidth = 8;
figheight = 12;
nrow = 3;
ncol = 1;
options.bmargin = 0.15;
options.visible = 'off';
options.units = 'centimeters';
options.axesfontsize = 10;

[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,options);
cmap = lines(4);

% truckcolor = [255 125 0]./255;
% dronecolor = [0.4660    0.6740    0.1880];
truckcolor = 'k';
dronecolor = 'k';

axes(ax(1))
hd = scatter(hovers,Sig_d,'MarkerEdgeColor',dronecolor,'LineWidth',1,'Marker','x');
hold on
eb(1) = errorbar(hovers,Sig_d,err_drone(:,1), 'vertical', 'LineStyle', 'none','color',dronecolor);

hold on
ht = scatter(0.1+hovers,Sig_t,'MarkerEdgeColor',truckcolor,'LineWidth',1,'Marker','s');
eb(2) = errorbar(0.1+hovers,Sig_t,err_truck(:,1), 'vertical', 'LineStyle', 'none','color',truckcolor);

% hst = scatter(0.05+hovers,Sig_S06_t,'MarkerEdgeColor',cmap(3,:),'LineWidth',1);
% hsd = scatter(0.05+hovers,Sig_S06_d,'MarkerEdgeColor',cmap(4,:),'LineWidth',1);

ylim([0.2 2])
ylabel('S_{ig} (m)')

% hl = legend([hd ht hsd hst],{'drone','truck','S06 \beta_d','S06 \beta_t'});
hl = legend([hd ht],{'hovering','stationary'});

hl.NumColumns = 1;
hl.Location = 'Northeast';
hl.Box = 'on';

% hl.FontSize = 8;
% hl.Position(2) = ax(1).Position(2)+ax(1).Position(4)+0.05;


axes(ax(2))

scatter(hovers,Sinc_d,'MarkerEdgeColor',dronecolor,'LineWidth',1,'Marker','x')
hold on
eb(1) = errorbar(hovers,Sinc_d,err_drone(:,2), 'vertical', 'LineStyle', 'none','color',dronecolor);

scatter(0.1+hovers,Sinc_t,'MarkerEdgeColor',truckcolor,'LineWidth',1,'Marker','s')
eb(2) = errorbar(0.1+hovers,Sinc_t,err_truck(:,2), 'vertical', 'LineStyle', 'none','color',truckcolor);

% scatter(0.05+hovers,Sinc_S06_t,'MarkerEdgeColor',cmap(3,:),'LineWidth',1)
% scatter(0.05+hovers,Sinc_S06_d,'MarkerEdgeColor',cmap(4,:),'LineWidth',1)

% ylim([0.2 1.2])
ylabel('S_{inc} (m)')


axes(ax(3))

scatter(hovers,eta_d-tideMean,'MarkerEdgeColor',dronecolor,'LineWidth',1,'Marker','x')
hold on
scatter(0.1+hovers,eta_t-tideMean,'MarkerEdgeColor',truckcolor,'LineWidth',1,'Marker','s')
% scatter(0.05+hovers,eta_S06_t,'MarkerEdgeColor',cmap(3,:),'LineWidth',1)
% scatter(0.05+hovers,eta_S06_d,'MarkerEdgeColor',cmap(4,:),'LineWidth',1)

ylim([0.1 1])
ylabel('setup (m)')
xlabel('Hover')

for i=1:3
    ax(i).XLim = [0.5 7.5];
    ax(i).Box = 'on';
    ax(i).Position(1) = 0.17;
    ax(i).Position(3) = 0.8;
    if i<3
        ax(i).XTickLabel = [];
    end
    ax(i).FontSize = 8;
end

ax(3).XTick = [1:7];
ax(3).XTickLabel = {'D1' 'D3' 'D4' 'F2' 'F3' 'F4' 'F5'};

ax(2).YTick = [0.4000:0.2:1.6000];


print(hFig, '-dpdf', ['../viz/paperFigures/runupBulkCompareALL_noS06.pdf'],'-r300');
