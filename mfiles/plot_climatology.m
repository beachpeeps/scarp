%plot figure 2: Cortez Climatology
% function plot_figure2
clear
figureDir = '~/Documents/Repositories/scarp/viz/';
figureName = 'climatology2';
% close all
% 

%top: HoLo vs frequency spread, colored by mean peak frequency
%bottom: Ho histogram

%for subplot a (top), need fpeak, sqrtHoLo, fspread, mshift
MOPfile = '~/Documents/Repositories/scarp/data/MOP582.mat';
[sqrtHoLo, mshift, fpeak, fspread, Ho] = prep_plot_climatology(MOPfile);

load(MOPfile)

monthnum = month(MOP.time);
%%
clf
subplot(2,1,1)
hist(MOP.Hs(monthnum == 2),100)
xlabel('Hourly Hs (m) from MOP 582')
ylabel('Counts')
% titlestr = ['February Wave Climate' datestr(MOP.time(1)  datestr(MOP.time(end))]
title('February Wave Climate MOP 582: 2000-2017')

subplot(2,1,2)
[f,x] = ecdf(MOP.Hs(monthnum == 2));
plot(x,f*100)
xlabel('Hourly Hs (m) from MOP 582')
ylabel('cdf')
per90 = x(find(f>0.9,1));
hold on
plot([per90 per90],[0 90],':k')
text(per90,50,['90% of all February waves are ','less than Hs = ' num2str(per90,'%2.2f') ' m'])

% 
% 
% set(0,'defaultaxesfontsize',10)
% 
% Z.label = 'Mean $f_p$ (Hz)';
% Z.labelInterpreter = 'latex';
% Z.array = fpeak;
% Z.binMethod = 'mean';
% Z.cmap = cmocean('haline');
% Z.scale = 'linear';
% 
% Y.label = '$\sqrt{H_0L_0}$ (m)';
% Y.labelInterpreter = 'latex';
% Y.array = sqrtHoLo;
% Y.bin = [0:1:50];
% Y.scale = 'linear';
% 
% X.label = 'Frequency Spread (Hz)';
% X.labelInterpreter = 'latex';
% X.array = fspread;
% X.bin = [0.01:0.005:0.12]; % I made these bounds up
% X.scale = 'linear';
% 
% HC.xtick = [0 1];
% HC.xticklabel= {'summer','winter'};
% HC.label =  'Avg. Season';
% HC.array = mshift;
% HC.binMethod = 'fraction';
% % 
% %%
% figwidth = 14;
% figheight = 23;
% nrow = 3;
% ncol = 2;
% units = 'centimeters';
% [hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,'units',units,'axesfontsize',10);
% hFig.Visible = 'off';
% 
% ax(1).Position = [0.1 0.80 0.6 0.13];
% ax(3).Position = [0.1 0.45 0.6 0.34];
% ax(4).Position = [0.71 0.45 0.24 0.34];
% ax(2).Visible = 'off';
% % 
% ax(13) = axes('Position',ax(3).Position);
% ax(13).Visible = 'off';
% 
% ax(5).Position([3 4]) = ax(3).Position([ 3 4]);
% ax(6).Position([1 3 4]) = ax(4).Position([1 3 4]);
% 
% % ax(6).Visible = 'off';
% 
% % print('-dpdf',[figureDir figureName '.pdf'],'-r300')
% 
% plot_jointpdfMOPS(X,Y,Z,HC,hFig,ax);
% 
% % %% Ho plots
% 
% Y.label = '$H_0$ (m)';
% Y.labelInterpreter = 'latex';
% Y.array = Ho;
% Y.bin = [0:0.1:5.5];
% Y.scale = 'linear';
% 
% 
% %BIN DATA
% yind = discretize(Y.array,Y.bin); % bin ydata
% xind = discretize(X.array,X.bin); % bin xdata
% ny = length(Y.bin);
% nx = length(X.bin);
% 
% % preallocate variables for binning
% zdata = nan(nx,ny);
% counts = nan(nx,ny);
% hdata = nan(nx,ny);
% indZ = cell(nx,ny);
% %
% for i=1:length(X.bin)
%     for j = 1:length(Y.bin)
%         zind = intersect(find(xind==i),find(yind==j));
%         counts(i,j) = numel(zind);
%         indZ(i,j) = {zind};
%         if numel(zind)>0
%             zdata(i,j) = nanmean(Z.array(zind));
%             hdata(i,j) = sum(HC.array(zind))./numel(zind); 
%         end
%     end
% end
% %
% 
% cmap = cmocean('-gray');
% axes(ax(5))
% pcolor(X.bin,Y.bin,zdata'); shading flat
% 
% ax(5).XLabel.Interpreter = X.labelInterpreter;
% ax(5).XLabel.String = X.label;
% ax(5).YLabel.Interpreter = Y.labelInterpreter;
% ax(5).YLabel.String = Y.label;
% 
% 
% colormap(ax(5),Z.cmap)
% ax(5).CLim = [0.04 0.12];
% ax(3).CLim = [0.04 0.12];
% 
% 
% axes(ax(6))
% hy = histogram(Y.array,'BinEdges',Y.bin);
% by = bar(Y.bin(1:end-1),hy.BinCounts,1);
% by.FaceColor = [0.5 0.5 0.5];
% 
% ax(6).View = [90 -90];
% ax(6).Box = 'off'; ax(6).Color = 'none';
% ax(6).YLim = ax(4).YLim;
% 
% 
% 
% 
% 
% 
% %% labeling
% text(0.9,0.93,'a','units','normalized',...
%         'fontsize',18,'parent',ax(1),'fontweight','bold')
% text(0.9,0.93,'b','units','normalized',...
%         'fontsize',18,'parent',ax(3),'fontweight','bold')    
% text(0.1,0.93,'c','units','normalized',...
%         'fontsize',18,'parent',ax(4),'fontweight','bold')     
% text(0.9,0.93,'d','units','normalized',...
%         'fontsize',18,'parent',ax(5),'fontweight','bold')
% text(0.1,0.93,'e','units','normalized',...
%         'fontsize',18,'parent',ax(6),'fontweight','bold')   
%     
%     
% 
% 
% % %% condo event star    
% % % add in condo event
% % MOPfile2 = '~/Documents/Cortez/data/processed/MOPdata.mat';
% % [sqrtHoLoC, mshift, fpeakC, fspreadC, HoC] = prep_plot_climatology(MOPfile2);
% % 
% % 
% % axes(ax(3))
% % hold on
% % hs = scatter(fspreadC(1),sqrtHoLoC(1),200,'pr');
% % hs.MarkerFaceColor = 'y';
% % 
% % axes(ax(5))
% % hold on
% % scatter(fspreadC(1),HoC(1),200,'pr','MarkerFaceColor','y');
% % 
% % %% define large waves on the plots (Hs>2m, HoLo>25m)
% % axes(ax(3))
% % hStorm(1) = plot(ax(3).XLim,[25 25],'r');
% % 
% % axes(ax(5))
% % hStorm(2) = plot(ax(3).XLim,[2 2],'r');
% % 
% % axes(ax(4))
% % hold on
% % hStorm(3) = plot([25 25],ax(4).YLim,'r');
% % 
% % axes(ax(6))
% % hold on
% % hStorm(4) = plot([2 2],ax(6).YLim,'r');
% 
% %%
% ax(5).XLim = ax(3).XLim;
% ax(6).YLabel.String = 'counts ($\times 10^4$)';
% ax(6).YLabel.Interpreter = 'latex';
% ax(1).YLabel.String = 'counts ($\times 10^4$)';
% ax(1).YLabel.Interpreter = 'latex';
% 
% ax(4).YTickLabel = [];
% ax(3).XTickLabel = [];
% ax(3).XLabel.String = {};
% % ax(1).YTickLabel{1} = ' ';
% % ax(3).YTickLabel{1} = ' ';
% ax(6).YRuler.Exponent = 0;
% ax(1).YRuler.Exponent = 0;
% ax(6).YAxis.TickLabels = {' ';'1';'2'};
% ax(1).YAxis.TickLabels = {' ';'1';'2'};
% 
% 
% ax(6).XLim = ax(5).YLim;
% ax(6).XScale = ax(5).YScale;
% ax(6).XTickLabel = [];
%     
% %% print    
% print('-djpeg',[figureDir figureName '.jpeg'],'-r300')
% print('-dpdf',[figureDir figureName '.pdf'],'-r300')
% % close all
% 
% function [zdata, hdata, ZBIN] = plot_jointpdfMOPS(X,Y,Z,HC,hFig,ax);
% % joinpdfMOPS(xdatalabel,xdata,ydatalabel,ydata,Z,hcolor,HC)
% % plots the joint pdf of Z as binned into X and Y, with histograms of values
% % used to form the bins on the x and y axes. Histogram bars are colored
% % by HC.
% %
% % INPUTS: 
% % fileName: location where output figure is to be saved
% % printYN: 0 if not saved, 1 if saved
% %
% % ***Input for the main plots (X,Y,Z)***
% % Z, structure with the following components
% %   .label = string of variable name, example: '$\sqrt (HoLo)$'
% %   .labelInterpreter = 'tex' or 'latex'
% %   .array = variable to put on the z axis, size(1, number of obs)
% %   .binMethod = 'max', 'min', or 'mean'
% %   .scale = 'log' or 'linear'
% % 
% % X/Y, structure with the following components:
% %   .label = string of variable name, example: 'Lowest f peak (Hz)'
% %   .labelInterpreter = 'tex' or 'latex'
% %   .array = variable to put on the x or y axis, size(1, number of obs)
% %   .bin = bins for the x/y axis, example: [245:3:295]
% %   .scale = 'log' or 'linear'
% % 
% % *** Input for the histogram colors (HC) ****
% % HC, structure with the following components
% %   .xtick = array to set bounds of colorbar, example: [0 1]
% %   .xticklabel = label for the colorbar ticks, ex: {'summer','winter'}
% %                 **note these must match the # of ticks!!** 
% % 	.label = string of variable name, example: 'Avg. Season'
% % 	.array = variable to color the histograms, size(1, number of obs)
% %   .binMethod = 'mean' or 'fraction'  
% %   .yLimit = maximum count number, for showing extremes
% %
% % OUTPUTS:
% % zdata: binned Z.array
% % hdata: binned HC.array
%  
% % Julia Fiedler, 2019 jfiedler@ucsd.edu
% 
% 
% %% BIN DATA
% yind = discretize(Y.array,Y.bin); % bin ydata
% xind = discretize(X.array,X.bin); % bin xdata
% nx = length(X.bin);
% ny = length(Y.bin);
% 
% % preallocate variables for binning
% zdata = nan(nx,ny);
% counts = nan(nx,ny);
% hdata = nan(nx,ny);
% indZ = cell(nx,ny);
% %
% for i=1:length(X.bin)
%     for j = 1:length(Y.bin)
%         zind = intersect(find(xind==i),find(yind==j));
%         counts(i,j) = numel(zind);
%         indZ(i,j) = {zind};
%         if numel(zind)>0
%             switch Z.binMethod
%                 case 'mean'
%                     zdata(i,j) = nanmean(Z.array(zind));
%                 case 'max'
%                     zdata(i,j) = nanmax(Z.array(zind));
%                 case 'min'
%                     zdata(i,j) = nanmin(Z.array(zind));
%                 case 'fraction'
%                     zdata(i,j) = sum(Z.array(zind))./numel(zind);
%             end
% 
%             switch HC.binMethod
%                 case 'mean'
%                     hdata(i,j) = nanmean(HC.array(zind));
%                 case 'fraction'
%                     hdata(i,j) = sum(HC.array(zind))./numel(zind);
%             end
%         end
%     end
% end
% %%
% 
% 
% % disp(['Making plot of ' Z.label ' as binned by ' X.label ' and ' Y.label]);
% 
% 
% 
% cmap = cmocean('-gray');
% axes(ax(3))
% pcolor(X.bin,Y.bin,zdata'); shading flat
% 
% ax(3).XLabel.Interpreter = X.labelInterpreter;
% ax(3).XLabel.String = X.label;
% ax(3).YLabel.Interpreter = Y.labelInterpreter;
% ax(3).YLabel.String = Y.label;
% 
% 
% colormap(ax(3),Z.cmap)
% 
% % axes(ax(13))
% % contour(X.bin,Y.bin,counts','linewidth',0.5)
% % ax(13).Color = 'none';
% % colormap(ax(13),gray);
% 
% axes(ax(1))
% hx = histogram(X.array,'BinEdges',X.bin);
% bx = bar(X.bin(1:end-1),hx.BinCounts,1);
% bx.FaceColor = [0.5 0.5 0.5];
% 
% 
% axes(ax(4))
% hy = histogram(Y.array,'BinEdges',Y.bin);
% by = bar(Y.bin(1:end-1),hy.BinCounts,1);
% by.FaceColor = [0.5 0.5 0.5];
% 
% % switch HC.binMethod
% %     case 'mean'
% %         colorindX = nanmean(hdata,2)/max(hdata(:));
% %         colorindY = nanmean(hdata)/max(hdata(:));
% %     case 'fraction'
% %         colorindX = nanmean(hdata,2);
% %         colorindY = nanmean(hdata);
% % end
% % 
% % colorindX(isnan(colorindX)) = 0;
% % colorindY(isnan(colorindY)) = 0;
% % 
% % colorindXY = [colorindX(:) ; colorindY(:)];
% % colorindXYind = round(colorindXY./max(colorindXY)*256);
% % colorindXYind(colorindXYind==0) = 1;
% % 
% % colorX = colorindXYind(1:length(colorindX)-1);
% % colorY = colorindXYind(length(colorindX)+1:end-1);
% % 
% % bx.CData = cmap(colorX,:);
% % by.CData = cmap(colorY,:);
% 
% ax(4).View = [90 -90];
% % if isfield(HC,'yLimit')
% %     ax(4).YLim(2) = HC.yLimit;
% % end
% 
% % ax(16) = axes('Position',[0.71 0.85 0.2 0.03]);
% % ax(16).Visible = 'off'; 
% % plot(hdata)
% % ax(16).CLim = [0 max(colorindXY)];
% % colormap(ax(16),cmap)
% % c = colorbar(ax(16),'location','northoutside');
% % c.Position = ax(16).Position;
% % ax(16).Position([3 4]) = [0 0];
% % ax(16).XTickLabel = []; ax(6).YTickLabel = [];
% % 
% % c.Ticks = ax(16).CLim;
% % c.TickLabels = HC.xticklabel;
% % % c.TickLabels = {num2str(c.Limits(1)), num2str(round(c.Limits(2),1),'%2.1f')};
% % c.Label.String = HC.label;
% % 
% c2 = colorbar(ax(3),'location','northoutside');
% c2.Position = [0.4    0.67    0.2700    0.0300];
% c2.Label.Interpreter = Z.labelInterpreter;
% c2.Label.String = Z.label;
% c2.TickDirection = 'out';
% 
% for i = [1 4]
% ax(i).Box = 'off'; ax(i).Color = 'none';
% end
% 
% ax(1).XTickLabel = [];
% for i= [1 3]
%     ax(i).XScale = X.scale;
% end
% 
% linkaxes([ax(1),ax(3)],'x')
% 
% 
% % in lieu of linkaxes([ax(1),ax(3),ax(13)],'y'):
% ax(3).YScale = Y.scale;
% % ax(13).YScale = Y.scale;
% 
% ax(4).XLim = ax(3).YLim;
% ax(4).XScale = ax(3).YScale;
% ax(4).XTickLabel = [];
% ax(4).YLabel.String = [];
% ax(1).YLabel.String = 'counts (x10^4)';
% 
% ax(3).ColorScale = Z.scale;
% end
% 
% 
