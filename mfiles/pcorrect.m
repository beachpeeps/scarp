function [P,H] = pcorrect(pressure,Hz,fcutoff,burial)
% code takes in pressure at depth, measured with pressure sensor, and
% returns linearly corrected surface pressure time series.

if nargin<2
    fcutoff = 0.33;
end

% get mean wave height
H = nanmean(pressure);

% give error for psensors not under water table
if H<0
    error('Warning: negative H detected - sensor is likely over watertable: skipping hour');
end

%% Apply depth correction to frequencies below fcutoff
% % Take into frequency space

if mod(length(pressure),2)==0
    pressure = pressure(1:end-1);
end

% infill nans, although nans should not exist...
pressure(isnan(pressure)) = nanmean(pressure);

%****** prewhiten pressure******
dp = diff(pressure);
pre = nanmean(pressure(1:5))*ones(600,1)+nanmean(dp(1:10)).*(1:600)'-nanmean(dp(1:10)).*600;
post = nanmean(pressure(end-5:end))*ones(600,1)-nanmean(dp(1:10)).*(1:600)';
ppad = [pre; pressure; post];
%*******************************

%****create frequency vector****
nn=length(ppad);
ny = floor(nn/2);

f = Hz/2*linspace(-1,1,nn);

fY = fft(ppad);
fY = fftshift(fY);
%********************************

% locate frequencies subject to correction
cutoffind = find(abs(f) <= fcutoff  & abs(f)>0);

% do correction in frequency/wavenumber space
k = ones(size(fY));
if H>0 % can only solve for k if water depth is greater than 0   
    ff = cutoffind;
    k(ff) = getk(f(ff),H);
    fY(ff) = fY(ff).*exp(abs(k(ff).*burial)).*cosh(k(ff)*H);  
end
%******************

%***bring back to time domain and remove padding
fN = ifftshift(fY);
Psspad = abs(ifft(fN));

P = Psspad(601:end-600);
% need to solve the ringing problem at jumps/bores!!
%% check timeseries if needed
% clf
% plot(P)
% hold on
% plot(pressure)

%% check spectra if needed
% clf
% nfft = 800; fs = Hz;
% [Amp,~] = calculate_fft2(P,nfft,Hz);
% 
% df = fs/(nfft-1);   % This is the frequency resolution
% nnyq = nfft/2 +1;
% 
% fm = [0:nnyq-1]*df;
% Spp = mean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution
% Spp = Spp(1:nnyq);
% 
% loglog(fm,Spp)
% 
% [Amp,~] = calculate_fft2(pressure,nfft,Hz);
% 
% df = fs/(nfft-1);   % This is the frequency resolution
% nnyq = nfft/2 +1;
% 
% fm = [0:nnyq-1]*df;
% Spp = mean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution
% Spp = Spp(1:nnyq);
% 
% hold on
% loglog(fm,Spp)