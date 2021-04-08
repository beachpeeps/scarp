function ratio = burialcorrection(G,v,gamma,n,kc,S,beta,TSparos,burial,f,k)
% load ../mat/Dec_H2_P2_check.mat

% function to make burial correction from original Yamamoto equation
%Constants
rho = 1028; %density kg m-2
g = 9.81; %gravity m s-2
% P = (TSparos+burial).*g.*rho;
P = (TSparos+burial).*g.*rho+101600;
disp(mean(P))

% h = (P(tind)-patm(tind))./rho/g - burial;

% 
% Hz_paros = 2; %2 Hz
% nfftP = 5*60*Hz_paros; % 5 minute windows
% [f, SppPP, ~, ~, nens, dof] = get_spectrum(P(tind), nfftP, Hz_paros, 0.05);
% f=f(2:end);
% h = (P(tind)-patm(tind))./rho/g - burial;
% k = getk(f,mean(h));

% SOURCES:
% Yamamoto et al 1978
% Raubenheimer et al 1998
% Guest and Hay 2016
% Michallet et al 2009
% Zhang et al. 2009 https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-18801.pdf

% omega = radian frequency
% k = wavenumber
% that satisfy omega^2 = g k tanh (kh), where h is water depth
% 
% other parameters
% 
% omegapp = omega''
% kp = k'
% kpp = k''
% 
% n = sediment porosity
% kc = hydraulic conductivity (aka coefficient of permeability)
% G = shear modulus of the porous matrix
% v = Poisson's ratio
% betaprime = effective compressability of the pore fluid
% 
% m = n*G*betaprime./(1-2v);
% z = positive downward!!



% % Parameter values
% % Raubenheimer et al 1998 values
% G = 2e8; %N/m^2, shear modulus Mei and Foda 1981
% v = 0.3; %Poisson's ratio
% gamma = 10^4;%kg m-2 s-2, unit weight of the pore fluid
% % n = 0.3; %sediment porosity
% % kc = 1e-4; %m/s, hydraulic conductivity
% 
% %Edited values
% n = 0.3; %sediment porosity
% kc = 4e-4; %m/s, hydraulic conductivity
% S = 0.985; % degree of saturation, Raubenheimer assumes full saturation
% beta = 2.34e-9; %Pa-1, or m2 N-1, from Michallet et al 2019, Table 2
% % beta = 1e-7; %Pa-1, or m2 N-1, from Michallet et al 2019, Table 2
% 

% Ratio of pressure at depth vs bed surface, Eq. 1 Guest and Hay (Yamamoto
% formula)
betaprime = S*beta+(1-S)./mean(P); % Pa-1, or m2 N-1, or m s2 kg-1, effective compressibility, 
% P is P(z), absolute pore-water pressure:
% where P(z) = rho*g*(h+z) + P_atm, where h is mean water depth and z is
% positive downward burial depth, here mean pressure at depth

m = n*G*betaprime./(1-2*v);

kappa = (1-v)./(1-2*v);

B = n*betaprime+ 1./(2*kappa*G);

alpha = gamma./ kc.* B;


omega = 2*pi*f;

omegap = alpha*omega;

omegapp = kappa*(omegap./k.^2);

kp = k.*(1 + 1i*gamma*omega./(kc .* k.^2) .* B).^0.5;

kpp = (kp-k)./k;

A = 1i*m*omegapp./(-kpp + 1i*(1+m)*omegapp);


kprime = k .* (1+(1i*gamma*omega) .* B ./ (kc*k.^2)).^0.5;

% P(z)/P_0 ratio, where P(z) is pressure at sediment depth, and P_0 is the
% bed surface pressure
ratio = (1-A).*exp(-k.*burial)+A .*exp(-kprime.*burial);
ratio = ratio.*conj(ratio);
ratio = real(ratio)';

% 
% clf
% loglog(fmp,SppD(1:301)./SppParos)
% 
% hold on
% hl = plot(f,1./ratio);
% h2 = plot(f,1./ratio2.^2);

