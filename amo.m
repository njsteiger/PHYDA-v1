function [amo_indx] = amo(SST,amask,lat,smth)
%AMO Compute the AMO index given input data
%   [amo_indx] = amo(SST,amask,acell,nrthlat,smth)
%       Input:  SST = annual SST field (lat,lon,time)
%               amask = Atlantic ocean mask (NaNs and 1s)
%               lat = (code currently assumes regular lat/lon grid)
%               smth = smoothing to apply (0=none, 10=decadal)
%       Output: amo_indx = AMO index
%
%   Nathan Steiger, June 2016

[rw,cl,tm]=size(SST);
SST=reshape(SST,rw*cl,tm);

if length(amask(:))~=rw*cl
    error('Size of input data not compatible...')
end

% MASK OUT THE SST DATA
SST=bsxfun(@times,SST,amask(:));

% DETREND? NOT NECESSARY FOR LAST MILLENNIUM?
%gmt=gmt(:)-mean(gmt(:)); % need to modify input to allow 'gmt'
%SST=bsxfun(@minus,SST,gmt');% use definition of Trenberth et al.
%amo_indx=detrend(amo_indx); % simply remove linear trend

% Cosine-lat weighting 
A=cosd(repmat(lat,[1 cl]));
%disp('AMO computation assumes regular lat/lon grid...')
% COMPUTE ATLANTIC AREA AVERAGE
SST=reshape(SST,rw,cl,tm);
amo_indx=wmean_a(SST,A);

% Apply smothing if desired...
if smth==10
    amo_indx=smooth(amo_indx,10);
end

end














