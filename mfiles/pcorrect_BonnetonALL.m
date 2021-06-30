% clear all
% load('../mat/paros.mat')
clearvars -except paros
% load('../mat/paros.mat')
load('../mat/sandThickness.mat')
%%
%establish new variable within paros
paros(10).etaCorrectedB18 = [];

Hz = 2;
fcutoff = 0.33;

for pnum = 3
    %preallocate
    etaCorrectedB18 = nan(size(paros(pnum).eta));
    
    for tnum=1:length(paros(pnum).t)
        
        pRaw = paros(pnum).eta(:,tnum);
        
        burialInd = find(dtimeHourly == dateshift(paros(pnum).t(tnum)+minutes(1),'start','hour'));
        
        burial = sandThickness10cmInterp(pnum,burialInd);
        if isempty(burial) || burial == 0
            burial = sandThickness1mInterp(pnum,burialInd);
%             burial = 0;
        end
        
        try
            [eta_SNL, ~] = pcorrect_Bonneton(pRaw,fcutoff,Hz,burial);
            eta_SNL(eta_SNL<burial) = paros(pnum).eta(eta_SNL<burial,tnum);
            
            etaCorrectedB18(:,tnum) = eta_SNL;
        catch err
            if ~contains(err.message,'negative H detected')
                disp(err.message)
                break
            end
                
        end
        
    end
    
    paros(pnum).etaCorrectedB18 = etaCorrectedB18;
end
%%
clf
pnum = 5;
plot(paros(pnum).etaCorrectedB18(1e6:5e6))
hold on
plot(paros(pnum).eta(1e6:5e6))

%%
save('../mat/paros.mat','paros','-v7.3')