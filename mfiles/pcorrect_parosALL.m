clear all
load('../mat/paros.mat')
load('../mat/sandThickness.mat')
%%
%establish new variable within paros
paros(10).etaCorrected = [];

Hz = 2;
fcutoff = 0.33;

for pnum = 1:10
    %preallocate
    etaCorrected = nan(size(paros(pnum).eta));
    
    for tnum=1:length(paros(pnum).t)
        
        pRaw = paros(pnum).eta(:,tnum);
        
        burialInd = find(dtimeHourly == dateshift(paros(pnum).t(tnum)+minutes(1),'start','hour'));
        
        burial = sandThickness10cmInterp(pnum,burialInd);
        if isempty(burial)
            burial = 0;
        end
        
        try
            [P,H] = pcorrect(pRaw,Hz,fcutoff,burial);
            etaCorrected(:,tnum) = P;
        catch err
            if ~contains(err.message,'negative H detected')
                disp(err.message)
                break
            end
                
        end
        
    end
    
    paros(pnum).etaCorrected = etaCorrected;
end
%%
clf
pnum = 2;
plot(paros(pnum).etaCorrected(1e6:5e6))
hold on
plot(paros(pnum).eta(1e6:5e6))

%%
save('../mat/paros.mat','paros','-v7.3')