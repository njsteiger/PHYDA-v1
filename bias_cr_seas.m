function [Mbc] = bias_cr_seas(Od,Md)
%BIAS_CR_SEAS Bias corrects the model input according to the given
% observational data. Bias correction happens for the seasonal cycle
%   [Mbc] = bias_cr_seas(Od,Md)
%   Where 
%       Od[space,time] = observational data
%       Md[space,time] = model data (to be corrected)
%
%       Mbc[space,time] = model bias-corrected data
%
% Nathan Steiger, July 2017, LDEO

if ndims(Od)>2 || ndims(Md)>2
    error('Only use one spatial dimension...')
end

% Detrend observational data for getting the climatology
%Od=detrend(Od')';

[ms,mm]=size(Md);
[os,om]=size(Od);
Md_m12=reshape(Md,ms,12,mm/12);
Od_m12=reshape(Od,os,12,om/12);
Md_ns=bsxfun(@minus,Md_m12,mean(Md_m12,3));% remove climatology
Od_clm=mean(Od_m12,3);% find obs climatology
Md_bc=bsxfun(@plus,Md_ns,Od_clm); % replace climatology
Mbc=reshape(Md_bc,ms,mm);	


end

