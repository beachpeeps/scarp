function [Spec,Info,Bulk,Tseries] = get_runupStatsLidar(modDir,varargin)
% [Spec,Info,Bulk,Tseries,Sur] = get_runupStats(modDir,inputDate,varargin)
% This function processes the gridded Lidar data from Torrey Pines to get
% the runup line, as defined by Stockdon et al, 2006.
%
% Input
% moddir:       file location where gridded lidar data are located
% inputDate:    not sure why this is needed except for axis labels
% varargin:     additional options for computing the output
%
% Output
% Spec:     frequency (f), spectra (S), spectra at the offshore
%           boundary (Sboundary), confidence intervals (Slo,Sup),
%           degrees of freedom (dof);
% Info:     frequency sampling (Hz), runup threshold (threshold),
%           thinLayer***, timestamp in datetime format (datahour), data
%           directory (grdDir, surDir)
% Bulk:     swashparams and confidence intervals, arranged as
%           {Sig,Sinc,eta}, foreshore beach slope (beta),
%           offshore significant wave height (Ho)
% Tseries:  timeseries of the runupline (T, Zrunup, Xrunup, idxrunup)
% Sur:      xshore distance (x), mean water level (etabar)
%
% TODO: Is the algorithm finding the correct runupline, especially in
% coarse resolution? Could probably use a little help here.
% TODO: Holdovers here from model code, this is tricky, need to clean up!!!

% default options
options.threshold = 0.1;
options.windowlength = 5; % in minutes, for spectra windows
options.dx = 0.1; % for gridding the runup onto something finer
options.g = 9.81;
% options.Tm = 10;
% options.wlev = 0;

options = parseOptions( options , varargin );
%% Load the processed lidar


% load(modDir,'filename','Processed','filedir');
load(modDir,'Processed');


%%
nt = length(Processed.t);
xx = Processed.x;
M = movmin(Processed.Zinterp2,1000,1,'omitnan'); % get moving minimum for foreshore
%%
RunupImage = nan(1,nt);
idxrunup = ones(1,nt);
Zrunup = nan(1,nt);
flag = zeros(1,nt);

%%

L2 = [-40 30; options.threshold options.threshold];
tol = 1e-6;

for ii=1:nt 
    %
    indM = ii-500:ii+500; % take minumum over a window of 1000 pts in time
    indM = indM(indM>1 & indM<nt ); % cut off window length at ends of tseries
    %     wlevtemp =Processed.Zinterp(i,:)-medfilt1(min(M(indM,:)),3);
    %     wlevtemp = medfilt1(wlevtemp,3);
    
    
    if ii>2 && ii<nt
        watline = median(Processed.Zinterp(ii-1:ii+1,:));
        wlevtemp2 = watline-medfilt1(min(M(indM,:)),3);
    else
            wlevtemp2 = Processed.Zinterp(ii,:)-medfilt1(min(M(indM,:)),3);
    end

    L1 = [xx;wlevtemp2];
    %     L1c = [xx;wlevtemp2];
    if ii>4 && ~isnan(RunupImage(ii-1))
        prev3 = nanmean(RunupImage(ii-4:ii-1));
        L2 = [prev3-5 prev3+5; options.threshold options.threshold];
    end
    runupline = InterX(L1,L2);
    
    if ~isempty(runupline)
        %         %make a moving window over which to look for runupline
        RunupImage(ii) = runupline(1);
        %         idxrunup(i) = find(abs(xx-round(runupline(1)*10)/10)<tol);
    else
        %
        %         L2 = [-10 45; options.threshold options.threshold];
        %         %expand the search area if i>5
        if ii>5 && ~isnan(nanmean(RunupImage(ii-4:ii-1)))
            prev3 = nanmean(RunupImage(ii-4:ii-1));
            L2 = [prev3-10 prev3+10; options.threshold+0.01 options.threshold+0.01];
            %         L2 = [-10 45; options.threshold options.threshold];
            
        elseif ii>5 && isnan(nanmean(RunupImage(ii-4:ii-1)))
            L2 = [-40 45; options.threshold options.threshold];
            
        end
        
        runupline = InterX(L1,L2);
        if ~isempty(runupline)
            RunupImage(ii) = runupline(1);
        else
            RunupImage(ii) = NaN;
        end
    end
    %
%         clf
%         plot(L1(1,:),L1(2,:))
%         hold on
%         plot(L2(1,:),L2(2,:))
% %             plot(L1c(1,:),L1c(2,:),'o')
%     
%     % %     pause
%     % %
%         xlim([-40 50])
%         ylim([0 1.5])
%     % %         pause
% %         plot(xx,Processed.Zinterp(i,:))
%         title(ii)
%         ii=ii+1;
%      %
%     pause(0.01)%
% %     pause
    
    
    
    clear wlevtemp L1 runupline
end

%if things do not pass the eye-test, we will hand draw the line!!!





%%
R = medfilt1(RunupImage,5);
RR = inpaint_nans(R);
idxR = nan(size(RR));
for i=1:nt
    Rint = round(RR(i).*10/10);
    if Rint<max(xx)
        idxR(i) = find(xx==Rint);
    else
        idxR(i) = find(xx == max(xx));
    end
end

Xrunup = RunupImage;
for i=1:nt
    Zrunup(i) = Processed.Zinterp(i,idxR(i));
end


Xrunup(RunupImage==1) = nan;
Zrunup(RunupImage==1) = nan;

%%

dt = nanmean(diff(Processed.t*24*60*60));
dt = floor(dt*10)./10; % want accuracy only to 0.1


ZZ = inpaint_nans(Zrunup);


% Get spectra
nfft = options.windowlength*60/dt; % "windowlength" minute chunks
[f, S, Slo, Sup,~,dof] = get_spectrum(detrend(ZZ), nfft, 1/dt, 0.05);

% find freqs in inc and ig
nIG = find(f>0.004 & f<=0.04);
nINC = find(f>0.04 & f<=0.25); % <--carefull, do we want to include 0.04 with SS or IG?? This should be globally declared somewhere...
df = f(2)-f(1);

%
eta = nanmean(Zrunup(:));
%% get beta
stdEta = nanstd(Zrunup);

maxEta = eta+2*stdEta;
minEta = eta-2*stdEta;

foreshore = nanmin(Processed.Zinterp2(:,nanvar(M)<0.01));
foreshorex = xx(nanvar(M)<0.01);

%super kluge rescue
rmInd = find(foreshorex<-30);
foreshore(rmInd) = [];
foreshorex(rmInd) = [];

% all x should be connected:
rmInd = find(diff(foreshorex)>1);

foreshore(rmInd) = [];
foreshorex(rmInd) = [];

% all x should be connected:
% rmInd = find(diff(foreshorex)>1);

% foreshore(rmInd) = [];
% foreshorex(rmInd) = [];

botRange = find(foreshore>=minEta & foreshore<=maxEta);

fitvars = polyfit(foreshorex(botRange), foreshore(botRange), 1);
beta = fitvars(1);


%%
% Sig = 4*sqrt(sum(SLidar(nIG)*df));
% Sinc = 4*sqrt(sum(SLidar(nINC)*df));
% S = sqrt(Sinc.^2 + Sig.^2);

[Sig, Sig_lo, Sig_up, ~ ]= getSWHebounds( S(nIG), dof, 0.95, df );
[Sinc, Sinc_lo, Sinc_up, ~ ]= getSWHebounds( S(nINC), dof, 0.95, df );


swashparams = [Sig Sinc eta];
swashparamsLO = [Sig_lo Sinc_lo eta];
swashparamsUP = [Sig_up Sinc_up eta];

Spec.f = f;
Spec.S = S;
Spec.Slo = Slo;
Spec.Sup = Sup;
Spec.dof = dof;

Info.Hz = 1/dt;
Info.threshold = options.threshold;
Info.datahour = datetime(Processed.t(3),'ConvertFrom','datenum');
Info.duration = seconds((Processed.t(end)-Processed.t(3))*24*60*60);
% Info.rawFilename = filename;
% Info.rawFiledir = filedir;
Info.processedFilename = modDir;

Bulk.swashparams = swashparams;
Bulk.swashparamsLO = swashparamsLO;
Bulk.swashparamsUP = swashparamsUP;
Bulk.swashParamsNames = {'Sig','Sinc','eta'};
Bulk.beta = beta;
Bulk.foreshore = foreshore;
Bulk.foreshoreX = foreshorex;

Tseries.T = Processed.t; % in seconds
Tseries.Zrunup = Zrunup;
Tseries.Xrunup = Xrunup;
Tseries.idxrunup = idxR;





%% HELPER FUNCTIONS
%%
    function options = parseOptions( defaults , cellString )
        
        %INPUT PARSING
        p = inputParser;
        p.KeepUnmatched = true;
        
        names = fieldnames( defaults );
        for ii = 1 : length(names)
            %
            addOptional( p , names{ii}, defaults.(names{ii}));
            %
        end
        %
        parse( p , cellString{:} );
        options = p.Results;
        %
    end
%
    function P = InterX(L1,varargin)
        %INTERX Intersection of curves
        %   P = INTERX(L1,L2) returns the intersection points of two curves L1
        %   and L2. The curves L1,L2 can be either closed or open and are described
        %   by two-row-matrices, where each row contains its x- and y- coordinates.
        %   The intersection of groups of curves (e.g. contour lines, multiply
        %   connected regions etc) can also be computed by separating them with a
        %   column of NaNs as for example
        %
        %         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
        %               y11 y12 y13 ... NaN y21 y22 y23 ...]
        %
        %   P has the same structure as L1 and L2, and its rows correspond to the
        %   x- and y- coordinates of the intersection points of L1 and L2. If no
        %   intersections are found, the returned P is empty.
        %
        %   P = INTERX(L1) returns the self-intersection points of L1. To keep
        %   the code simple, the points at which the curve is tangent to itself are
        %   not included. P = INTERX(L1,L1) returns all the points of the curve
        %   together with any self-intersection points.
        %
        %   Example:
        %       t = linspace(0,2*pi);
        %       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
        %       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
        %       P = InterX([x1;y1],[x2;y2]);
        %       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')
        
        %   Author : NS
        %   Version: 3.0, 21 Sept. 2010
        
        %   Two words about the algorithm: Most of the code is self-explanatory.
        %   The only trick lies in the calculation of C1 and C2. To be brief, this
        %   is essentially the two-dimensional analog of the condition that needs
        %   to be satisfied by a function F(x) that has a zero in the interval
        %   [a,b], namely
        %           F(a)*F(b) <= 0
        %   C1 and C2 exactly do this for each segment of curves 1 and 2
        %   respectively. If this condition is satisfied simultaneously for two
        %   segments then we know that they will cross at some point.
        %   Each factor of the 'C' arrays is essentially a matrix containing
        %   the numerators of the signed distances between points of one curve
        %   and line segments of the other.
        
        %...Argument checks and assignment of L2
        narginchk(1,2);
        if nargin == 1
            L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
        else
            L2 = varargin{1}; hF = @le;
        end
        
        %...Preliminary stuff
        x1  = L1(1,:)';  x2 = L2(1,:);
        y1  = L1(2,:)';  y2 = L2(2,:);
        dx1 = diff(x1); dy1 = diff(y1);
        dx2 = diff(x2); dy2 = diff(y2);
        
        %...Determine 'signed distances'
        S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
        S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
        
        C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
        C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
        
        %...Obtain the segments where an intersection is expected
        [i,j] = find(C1 & C2);
        if isempty(i),P = zeros(2,0);return; end;
        
        %...Transpose and prepare for output
        i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
        L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
        i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
        
        %...Solve system of eqs to get the common points
        P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
            dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
        
        function u = D(x,y)
            u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
        end
    end

%calculate spectrum
    function [fm, Spp, Spplo, Sppup, nens, dof] = get_spectrum(P, nfft, fs,alpha)
        % [fm, Spp, Spplo, Sppup, nens, dof] = get_spectrum(P, nfft, fs,alpha)
        % Makes a spectrum for a single point only
        %
        % Input:
        %   P: timeseries to be analyzed (should be size (nt,1))
        %   nfft: window size for spectral averaging
        %   fs: sampling rate (Hz)
        %   alpha: for confidence intervals (alpha = 0.05 is standard 95% Conf Int)
        %
        % Output:
        %   fm: frequency
        %   Spp: Spectra (inputUnits^2/Hz)
        %   Spplo/Sppup: lower and upper bounds of confidence intervals
        %   nens: number of ensembles used
        %   dof: degrees of freedom
        
        % copyright 2019 J Fiedler jfiedler@ucsd.edu
        
        
        [m,n] = size(P);
        if m<n
            P = P';
        end
        
        
        [Amp,nens] = calculate_fft2(P(:,1),nfft,fs);
        nens2 = 0;
        
        % TODO: FIX THIS, FOR NOW only use one-dim P
        % if n==2
        % [Amp2,nens2] = calculate_fft2(P(:,2),nfft,fs);
        % Amp = [Amp; Amp2];
        % end
        
        df = fs/(nfft-1);   % This is the frequency resolution
        nnyq = nfft/2 +1;
        
        fm = [0:nnyq-1]*df;
        Spp = mean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution
        Spp = Spp(1:nnyq);
        
        Spp=Spp(:);
        
        nens = nens+nens2;
        % Confidence Intervals (CI)
        dof = 8/3*nens; % for Hanning window
        
        a = dof.*sum(Spp).^2./(sum(Spp.^2)); % effective dof
        adof = a/(1+2/dof); % unbiased edof
        
        chi2 = [chi2inv(alpha/2,dof) chi2inv(1-alpha/2,dof)];
        CI_chi2 = [(1/chi2(1)).*(dof*Spp) (1/chi2(2)).*(dof*Spp)];
        Spplo = CI_chi2(:,1);
        Sppup = CI_chi2(:,2);
        
        
    end

% get the error bounds using effective degrees of freedom
    function [ Hs, lb, ub, edof ] = getSWHebounds( e, dof, q, df )
        %[ lb, ub ] = getSWHebounds( e, dof )
        %   e - Energy(time,freq), units: m^2
        %   dof - degree of freedom per freq-band
        %   q - percentile to determine ub,lb (e.g. 0.90 or 0.68)
        %   lb - lower Hs bound
        %   ub - upper Hs
        %   Hs - Hsig
        %   edof - Effective dof used
        % SEE ELGAR 1987, for details of un-biased EDOF estimate
        
        % percent on either side
        a = (1-q)/2;
        
        % Get Hsig
        Hs = 4*sqrt(sum(e.*df));
        
        % Determine edof
        edof = dof*sum(e).^2./sum(e.^2);    %estimate effective DOF
        edof = edof/(1+2/dof);                    %UNBIAS THE ESTIMATE
        
        % Generate normalized ub,lb
        
        lb = edof/chi2inv(1-a,edof);
        ub = edof/chi2inv(a,edof);
        
        % Multiply by Hsig
        lb = sqrt(lb).*Hs;
        ub = sqrt(ub).*Hs;
        
        % % Multiply by Hsig
        % lb = lb.*Hs;
        % ub = ub.*Hs;
        
    end
end
