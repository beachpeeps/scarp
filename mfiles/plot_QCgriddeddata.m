%plot_QCgriddeddata
%not gonna do a loop because whatever
hovern = 2;
vizdir = '../viz/';
%because everything has a different name
filename = ['20200224_00582_TorreyRunup_H' num2str(hovern) '_10cm_ParLot.mat'];
% filename = ['20191214_H' num2str(hovern) '_navd88_geoid12b_10cm_ParLot.mat'];
load(['../mat/lidar/drone/' filename]);

%% timestack figure
hFig = figure;
pcolor(Processed.x,Processed.t(1:1:3000),Processed.Zinterp2(1:1:3000,:));
shading flat
xlabel('m from lidar')
datetick('y','MM:SS')
ylabel('Time (MM:SS)')
title(filename, 'Interpreter', 'none');
hc = colorbar;
caxis([0.5 10])
hc.Label.String = 'Z ';
print(hFig, '-djpeg', [vizdir filename '_timestack_ParLot.jpg'],'-r300');


%% wiggle figure
hFig = figure;

stable_xshore = [45;50]; %m
[wiggle, ~] = remove_wiggles(Processed,stable_xshore);


xind = [10 30 40];

for i=1:length(xind)
    plot(Processed.t,Processed.Zinterp2(:,xind(i))-nanmean(Processed.Zinterp2(:,xind(i))))
    hold on

    ss{i} = [num2str(Processed.x(xind(i))) 'm'];
end

plot(Processed.t,wiggle,'k','linewidth',2)
ylim([-0.3 0.3])
legend([ss 'wiggle'] )
ylabel('Drift: Z-mean(Z) (m)')
datetick('x','MM:SS')
xlabel('Time (MM:SS')
title(filename, 'Interpreter', 'none');
% print(hFig, '-djpeg', [vizdir filename '_wiggle_ParLot.jpg'],'-r300');

t = Processed.t;

%% save the wiggles save the world
save(['../mat/20200224_H' num2str(hovern) '_wiggle_35.mat'],'wiggle','stable_xshore','t')
