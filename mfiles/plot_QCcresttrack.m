cmap = lines(2);
clf
v = VideoWriter('../viz/cresttrack.avi');
v.FrameRate = 4;
open(v)
for i=2:45
    clf

    wavex = wave_x(i,~isnan(wave_x(i,:)));
    waveh = wave_h(i,~isnan(wave_x(i,:)));
    
    wavet = wave_t(i,~isnan(wave_x(i,:)));
    [a,b] = unique(wavex);
    %
    wt = dt.*interp1(a,wavet(b),x(1:xlimit));

%model is x = a*t^2 + b*t + c
%so that Cp is dx/dt = 2*a*t + b

[p,stats] = polyfit(wavet(b)*dt,a,2);
polyFit = polyval(p,wt);
dp = p.*[2 1 0];
cpFit = dp(1)*wt+dp(2);

[p2,stats] = polyfit(wavet(b)*dt,a,3);
polyFit2 = polyval(p2,wt);
dp = p2.*[3 2 1 0];
cpFit2 = dp(1)*wt.^2+dp(2)*wt+dp(3);

% if dp(1)>0
%     cpFit = nan(size(cpFit));
% end



% temporary plot checkpoint
subplot(4,1,1:2)
tchunk = round(min(wt)):0.1:round(min(wt))+30;
pcolor(x(1:xlimit),tchunk,eta(1:xlimit,tchunk*10)'); shading flat
hold on
h1 = plot(x(1:xlimit),wt,'color',[0 0 0 0.5]);
hold on
h2 = plot(polyFit,wt,'color',cmap(1,:));
h3 = plot(polyFit2,wt,'color',cmap(2,:));
ylabel('time (sec)')
title('Crest track')
legend('lidar','raw track','2nd degree','3rd degree','location','best')

xlim([-100 -10])
subplot(4,1,3)
plot(x(1:xlimit),cpFit,'color',cmap(1,:))
hold on
plot(x(1:xlimit),cpFit2,'color',cmap(2,:))
title('phase speed (m/s)')
ylabel('C_p (m/s)')
xlim([-100 -10])
subplot(4,1,4)
plot(x(1:xlimit),-cpFit.^2/9.81,'color',cmap(1,:))
hold on
plot(x(1:xlimit),-cpFit2.^2/9.81,'color',cmap(2,:))
xlabel('x (m)')
ylabel('h (m)')
title('Calculated depth - linear')
xlim([-100 -10])
frame = getframe(gcf)
writeVideo(v,frame)
if i<50
pause(0.1)
end
% delete(h1)
% delete(h2)
end



