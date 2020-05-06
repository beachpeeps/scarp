clear all
load('../mat/paros.mat')
load('../mat/sandThickness.mat')
%%
%establish new variable within paros
paros(10).etaCorrected = [];

Hz = 2;
fcutoff = 0.25;

for pnum = 1:10
    %preallocate
    etaCorrected = nan(size(paros(pnum).eta));
    
    for tnum=1:length(paros(pnum).t)
        
        pRaw = paros(pnum).eta(:,tnum);
        
        burialInd = find(dtimeHourly == paros(pnum).t(tnum));
        
        burial = sandThickness10cmInterp(pnum,burialInd);
        if isempty(burial)
            burial = 0;
        end
        
        try
            [P,H] = pcorrect(pRaw,Hz,fcutoff,burial);
            etaCorrected(:,tnum) = P;
        catch
        end
        
    end
    
    paros(pnum).etaCorrected = etaCorrected;
end
%%
clf
pnum = 8;
plot(paros(pnum).etaCorrected(100000:250000))
hold on
plot(paros(pnum).eta(100000:250000))

%%
save('../mat/paros.mat','paros','-v7.3')