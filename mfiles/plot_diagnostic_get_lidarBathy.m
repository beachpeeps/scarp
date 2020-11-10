%plot_diagnostic_get_lidarBathy
figwidth = 5;
figheight = 8;
nrow = 4;
ncol = 1;
units = 'inches';
cbar = 0;
[hFigDiagnostic, ax] = makeFig(figwidth,figheight,nrow,ncol,units,cbar,'hspace',0.01);
hFigDiagnostic.Color = 'w';

for i=1:length(ax)
    hold(ax(i),'on')
    ax(i).XLim = [x(1) x(end)];
end


axes(ax(1))
plot(wave_x,1:max_time,'color',[0.8 0.8 0.8 0.4],'parent',ax(1));

hwave = plot(wave_x(nwaveplot,:),1:max_time,'r','parent',ax(1));

legend(hwave,['wave traj.' num2str(nwaveplot)])


axes(ax(2))
hcp = plot(x(1:xlimit),Cp,'.','parent',ax(2),'color',[0.8 0.8 0.8 0.4]);
h1 = plot(x(1:xlimit),Cp(nwaveplot,:),'r','parent',ax(2));

hm = plot(x(1:xlimit),Cpm,'m','parent',ax(2));

legend([hcp(1) h1 hm],{'C_p',['C_p wave ' num2str(nwaveplot)],'median C_p'})

axes(ax(3))
plot(x(1:xlimit),-h_linear,'parent',ax(3));
plot(x(1:xlimit),-h_bore,'parent',ax(3));

legend('linear','with bore correction')


axes(ax(4))
    r = nan(1,xlimit);
    for i=1:xlimit
        if sum(~isnan(Cp(:,i)))>4
            val = find(~isnan(Cp(:,i)));
            
            mdl = corrcoef(Cp(val,i),Cp_H(val,i));
            r(i) = mdl(2);
        end
    end
    
    plot(x(1:xlimit),r,'.','parent',ax(4))
    hold on
    plot(x, zeros(size(x)),':k','parent',ax(4))
    title('Correlation of C_p and H')






ax(1).Title.String = ['wave trajectories and slope analysis region for wave ' num2str(nwaveplot)];
ax(2).YLabel.String = 'phase speed (m/s)';
ax(4).XLabel.String = 'x (m)';

ax(1).YLabel.String = 'timestep (0.1s)';

ax(3).YLabel.String = 'h (m)';
ax(4).YLabel.String = 'Correlation r';

ax(3).YLim = [-3 0];
ax(3).Legend.Location = 'southeast';


ax(2).YLim = [0 10];
ax(4).YLim = [-1 1];
ax(4).Title.Position(2) =  0.5;

for i=1:length(ax)-1
    ax(i).XTickLabel = [];
end

ss = ['../viz/' hoverdate '_H' num2str(hovern) '_lidarBathyTruck_diagnostic_cresttrack.png'];

print(hFigDiagnostic,'-dpng','-r300',ss)