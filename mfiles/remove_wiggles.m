function [wiggle, Processed] = remove_wiggles(Processed,stable_xshore)
% remove wiggles from GPS funk
% correction = remove_wiggles(filename,stable_xshore)
% INPUT:
%   filename: mat file of gridded drone hover
%   stable_xshore: xshore locations (in m) of stable substrate to check for
%   wiggles. should be 2 element vector of min and max stable area.
% OUTPUT:
%   correction: wiggle (in m) with same time vector as gridded drone hover
%   Processed: will automatically apply the correction to the Processed
%   structure, so include carefully!!!
%
% example input:
% filename = '/Volumes/FiedlerBot8000/scarp/mat/lidar/drone/20200224_00582_TorreyRunup_H3_10cm_ParLot.mat';
% stable_xshore = [ 40 55 ]; m, endpoints of stable area to look at
% load(filename,'Processed')

xlocs = knnsearch(Processed.x',stable_xshore); %find indices of stable_xshore
xlocs = xlocs(1):xlocs(end); % just in case stable_xshore>2 pts

Zneat_movingmean = movmedian(Processed.Zinterp(:,xlocs),50,'omitnan'); % 5 second moving mean, TODO: change for non 10Hz collection

Zneat = Zneat_movingmean-nanmean(Zneat_movingmean);

wiggle = nanmedian(Zneat,2);

%THIS WILL CHANGE THE Processed INPUT!!! TAKE NOTE!!!
Processed.Zinterp = Processed.Zinterp-wiggle;
Processed.Zinterp2 = Processed.Zinterp2-wiggle;

