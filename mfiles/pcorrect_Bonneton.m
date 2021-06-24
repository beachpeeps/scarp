function [eta_SNL, eta_SNL_discrete] = pcorrect_Bonneton(eta,f_c,Hz,burial)
% Bonneton et al. 2018 "A nonlinear weakly dispersive method for recovering
% the elevation of irrotational surface waves from pressure measurements."

%***** variable definitions *****%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z = eta(x,t) elevation of the surface above the still water level
% z = 0, still water level
% z = -h_0, constant bottom elevation
% P_m(x,z,t) pressure measured at a distance delta_m above the bottom
% delta_m, distance off the bottom

% eta_H, hydrostatic reconstruction
% eta_SL, linear shallow water reconstruction
% eta_SNL, nonlinear shallow water reconstruction

% equations 11,12,13:
% eta_H = (P_m-P_atm)./(rho g) - h_0 + delta_m                      (11)
% eta_SL = eta_H - h_0/(2g)*(1-(delta_m/h_0)^2)*(eta_H)_tt          (12)
% eta_SNL = eta_SL - (1/g)( [eta_SL*(eta_SL)_t]_t -
%           (delta_m/h_0)^2*((eta_SL)_t)^2 )                        (13)  

% Note that since we bury our pressure sensors, delta_m = 0, so that
% changes our equations to after accounting for burial:
% eta_H = (P_m-P_atm)./(rho g) - h_0                                (11)
% eta_SL = eta_H - h_0/(2g)*(eta_H)_tt                              (12)
% eta_SNL = eta_SL - (1/g)( [eta_SL*(eta_SL)_t]_t )                 (13)  

%***********************************************************%


eta = eta(:)'; %make sure eta is oriented 1 x ###
[~,n] = size(eta); norig = n;
dt = 1/Hz; t = (0:dt:(n-1)/2); g = 9.81; %m/s^2


% do a lowpass filter, with cutoff frequency f_c
eta = [fliplr(eta) eta fliplr(eta)];
eta = lowpass(eta,f_c,Hz,'Steepness',0.6); % the steepness factor is chosen for slower rolloff
eta = eta(n+1:n*2);


% make even time series
if mod(length(eta),2) ==1
    eta= eta(1:end-1);
    t1 = t(1:end-1);
elseif mod(length(eta),2) == 0
    t1 = t;
end

% new length
nn = length(eta);

% mean water depth
h_0 = mean(eta);


%% ***** BURIAL CORRECTION *******
% This follows Raubenheimer equation 3 for fully saturated sands and soils
% in the limit of beta >> G (compressibility >> shear modulus), where 
% P_z/P_0 = exp(-kz)
% To do this, we will correct in frequency space

%****create frequency vector****

lpad = length(eta);
pad = [fliplr(eta) eta fliplr(eta)];
fY = fft(pad);
fY = fftshift(fY);
n = length(pad);
f = 2*pi/(n/Hz)* ( -n/2 : n/2-1 );

%********************************

% locate frequencies subject to correction
cutoffind = find(abs(f) <= 0.5  & abs(f)>0);

% do correction in frequency/wavenumber space
k = ones(size(fY));
if h_0>0 % can only solve for k if water depth is greater than 0   
    ff = cutoffind;
    k(ff) = getk(f(ff),h_0);
%     fY(ff) = fY(ff).*exp(abs(k(ff)).*burial);  
    fY(ff) = fY(ff).*exp(k(ff).*burial);  

end
%*

%***bring back to time domain 
fN = ifftshift(fY);
eta_H = abs(ifft(fN));
eta_H = eta_H(lpad+1:lpad*2); % remove padding

% **************END BURIAL CORRECTION************
% eta_H = eta;
%%
h_0 = mean(eta_H)-burial;
% h_0 =mean(eta_H);
[~,eta_H_tt_fourier] = padFourier(eta_H,Hz);


%% do discrete gradients as in Appendix A. Local time discretization
eta_H_t_discrete= zeros(size(eta_H));

eta_H_tt_discrete= zeros(size(eta_H));
for n=2:nn-1
    eta_H_t_discrete(n) = (eta_H(n+1)-eta_H(n-1))./(2/Hz);
    eta_H_tt_discrete(n) = (eta_H(n+1)-2*eta_H(n)+eta_H(n-1))./(1/Hz)^2;
end


eta_SL_discrete = eta_H - h_0/(2*g).*(eta_H_tt_discrete);
eta_SL_fourier = eta_H - h_0/(2*g).*(eta_H_tt_fourier);


eta_SL_t_fourier = padFourier(eta_SL_fourier,Hz);

eta_SL_t_discrete= zeros(size(eta_H));
eta_SL_eta_SL_t_t = zeros(size(eta_H));

for n=3:nn-2
    eta_SL_t_discrete(n) = (eta_SL_discrete(n+1)-eta_SL_discrete(n-1))./(2/Hz);
    eta_SL_eta_SL_t_t(n) = (eta_SL_discrete(n+1).*(eta_SL_discrete(n+2)-eta_SL_discrete(n))-eta_SL_discrete(n-1).*(eta_SL_discrete(n)-eta_SL_discrete(n-2)))./(4*(1/Hz).^2);
end

xout = padFourier((eta_SL_fourier.*eta_SL_t_fourier),Hz);

eta_SNL = eta_SL_fourier - (1/g).*xout;
eta_SNL_discrete = eta_SL_discrete - (1/g).*eta_SL_eta_SL_t_t; 

% replace the removed data point with uncorrected eta for truncated tseries 
if norig>nn
    eta_SNL(end+1) = eta(end);
    eta_SNL_discrete(end+1) = eta(end);
end



