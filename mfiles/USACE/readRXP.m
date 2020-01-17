function [xyzit,scanNum,numBadGPSpts]=readRXP(filenames,rotationMatrix4x4,varargin)
%% Read GPS reads in all the data from the filenames as xyzit and rotates the data using rotationMatrix4x4

% Optionally pass in the path to the RXPconvert executable
if ~isempty(varargin)
    exe = varargin{1};
else
    exe = 'RXPconvert.exe';
end

data=getRXPdata(filenames,'allpoints','x','y','z','r','tgps_datenum', exe)';

%% Parse Into Scan Number
% calculate scan number by looking at change in theta angles
x_socs=data(:,1);
y_socs=data(:,2);
z_socs=data(:,3);

theta=atan2(sqrt((x_socs).^2+(y_socs).^2),z_socs)*180/pi;
scanNum=nan(length(theta),1);
ind=[0; find(diff(theta)<-10); length(theta)];
for i=1:length(ind)-1
    scanNum(ind(i)+1:ind(i+1))=i;
end

%% Coordinate Transform
xyz1=[data(:,1) data(:,2) data(:,3) ones(size(data(:,1)))];
xyzR=rotationMatrix4x4*xyz1';
x=xyzR(1,:);
y=xyzR(2,:);
z=xyzR(3,:);
%% Handle tGPS dropouts using filename and tinternal
tgps=data(:,5);

numBadGPSpts=sum(isnan(tgps));
tgpsGood=tgps;
if numBadGPSpts>0 % at least some are bad
    if sum(~isnan(tgps(:)))==0 % all is bad
        %%
        %calculate a tGps for each using the filename
        fprintf('Using internal time and filename to fill tgps gaps\n\n');
        
        tgpsAllFromFilename=nan(size(tgps)); %*% should preallocate intelligently
        prevInd=1;
        if iscell(filenames) 
            %if filenames is a structure, need to assemble tint separately
            for iFilename=1:numel(filenames)
                tint=getRXPdata(filenames{iFilename},'last','tint', exe);
                fnameDatenum=filename2Datenum(filenames{iFilename});
                tgpsFromFilename=fnameDatenum+tint./(60*60*24);
                
                ind=numel(tgpsFromFilename);
                tgpsAllFromFilename(prevInd:prevInd+ind-1)=tgpsFromFilename;
                prevInd=prevInd+ind;
                
            end
        else
            tint=getRXPdata(filenames,'last','tint', exe);
            fnameDatenum=filename2Datenum(filenames{1});
            tgpsAllFromFilename=fnameDatenum+tint./(60*60*24);
        end
        tgpsGood=tgpsAllFromFilename;
        
    else %some is bad
        %%
        % use existing points with gps and internal time to calculate a
        % mean offset.
        fprintf('Using internal time and sparse gps time to fill tgps gaps\n\n');

        tintAll=nan(size(tgps)); %*% should preallocate intelligently
        prevInd=1;
        if iscell(filenames)
            %if filenames is a structure, need to assemble tint separately
            for iFilename=1:numel(filenames)
                tint=getRXPdata(filenames{iFilename},'last','tint', exe);
                ind=numel(tint);
                tintAll(prevInd:prevInd+ind-1)=tint;
                prevInd=prevInd+ind;
            end
        else
            tintAll=getRXPdata(filenames,'last','tint', exe)';
        end
        %convertTintAll to datenum 
        %keep sig figs by using big offset
        bigOffset=floor(nanmin(tgps));
        
        tgpsDaySeconds=(tgps-bigOffset)*60*60*24;
        
        meanOffsetSeconds=nanmean(tgpsDaySeconds-tintAll);
        
        tgpsFromInternalDayOffset=(tintAll+meanOffsetSeconds)./(60*60*24);
        
        tgpsFromInternal=(tgpsFromInternalDayOffset)+bigOffset;
        
        tgpsGood=tgpsFromInternal;
    end
end
%meanAccuracy=nanmean(tgpsGood-tgps)*60*60*24 %for accuracy debugging

%% Assemble xyzit variable
xyzit=[x(:) y(:) z(:) data(:,4) tgpsGood]';
end