function [proxy] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,obs_dtp1,mdl_dtp)
%PSM_UL  Univariate linear PSM

addpath('./regem-imputation')

% Create proxy data structure
proxy=struct([]);

% Throw out proxy data that doesn't overlap instrumental period
% with at least 20 points
%clb_prd=[1872 2010];
c_yrs=clb_prd(1):clb_prd(2);
[~,ip] = intersect(year,c_yrs); % overlap years of proxy data
% Find indices to remove
rmv=find(sum(~isnan(prxs(ip,:)))<20);
if ~isempty(rmv);disp(['Less than 20 regression values for ' num2str(length(rmv)) ' proxies...']);end
% INSTEAD OF REMOVING, JUST REPLACE WITH NANS: this will make the
% regression give NaNs for error term, thus pulling it out of contention to
% be picked and skipped over in the DA routine if that particular proxy is
% not possible to regress against

% prxs(:,rmv)=[];
% plat(rmv)=[];
% plon(rmv)=[];
% pNum=length(plat);
% archv_tp(rmv)=[];
prxs(:,rmv)=NaN*prxs(:,rmv);


%-------------------------
% LINEAR STATISTICAL PSM
%-------------------------



% LOAD OBSERVATIONAL DATA
%obs_dtp1='bearth'; % particular data type specified by input
[O_t,xlat,xlon,x_yrs] = load_obs_S(obs_dtp1,mnavg_p_o,mnavg_p_f);


% Get locations of proxies in obs
if ismember(obs_dtp1,{'bearth','hadcrut4'})
    [pl_o] = nearest_latlon(xlat,xlon,plat,plon);
else % mask out land or ocean
    mask=sum(~isnan(O_t),2);
    mask(mask==0)=NaN;
    mask(~isnan(mask))=1;
    [pl_o] = nearest_latlon(xlat,xlon,plat,plon,mask);
end


% LOAD MODEL DATA
if ismember(mdl_dtp,['s','y']) % load different lat/lon for ocean
    [M_t,xmeta]=load_S(mdl_prior,mdl_dtp,mnavg_p_o,mnavg_p_f,p_yrs);
    olat=xmeta{1}.olat;
    olon=xmeta{1}.olon;
else
    [M_t,xmeta]=load_S(mdl_prior,mdl_dtp,mnavg_p_o,mnavg_p_f,p_yrs);
    mlat=xmeta{1}.lat;
    mlon=xmeta{1}.lon;
end

%else % THIS IS JUST FOR THE PAGES2K RECONSTRUCTIONS
%   if ismember(mdl_dtp,['s','y']) % load different lat/lon for ocean
%       [M_t,~,~,~,olat,olon]=load_S_nobiascorr(mdl_prior,mdl_dtp,mnavg_p_o,mnavg_p_f,p_yrs);
%   else
%       [M_t,mlat,mlon]=load_S_nobiascorr(mdl_prior,mdl_dtp,mnavg_p_o,mnavg_p_f,p_yrs);
%   end
%end

% Locations of proxies in the model
if strcmp(mdl_dtp,'t')
    [pl_m] = nearest_latlon(mlat,mlon,plat,plon);
elseif ismember(mdl_dtp,['s','y']) % for ocean variables
    mask=sum(~isnan(M_t),2);
    mask(mask==0)=NaN;
    mask(~isnan(mask))=1;
    [pl_m] = nearest_latlon(olat,olon,plat,plon,mask);
else % mask out land or ocean
    mask=sum(~isnan(M_t),2);
    mask(mask==0)=NaN;
    mask(~isnan(mask))=1;
    [pl_m] = nearest_latlon(mlat,mlon,plat,plon,mask);
end

% Years of overlap
[~,io] = intersect(x_yrs,c_yrs);

%---------------------------------------------------
% Bias-correct model output at proxy locations only
%---------------------------------------------------

%O_tf=O_t(pl_o,:); % use all data for bias correction
%O_tf=O_t(pl_o,io); % use only calibration data for bias correction


% UNCOMMENT BELOW IF YOU WANT TO BIAS CORRECT BOTH THE MEAN AND VARIANCE
% c2_yrs=1950:1999;[~,io2] = intersect(x_yrs,c2_yrs);disp('---Using 1950:1999 for bias correction---')
% O_tf=O_t(pl_o,io2); % use only calibration data for bias correction
%
% if sum(isnan(O_tf(:))>0)
%     % infill random missing values in obs data with regem
%     options.maxit=10; options.regress='iridge';%options.regress='ttls';
%     % quite a few missing values in PDSI obs and 'ttls' isn't stable
%     if strcmp(obs_dtp1,'dai_pdsi');options.maxit=10;end %options.regress='iridge';end 
%     O_tf=regem(O_tf,options);
% end
% Ensures that the model data has same mean and variance as obs, so that the linear fit is meaningful
%  (Issues surrounding whether some obs data is anomalies vs real values) 
%[M_tbc] = bias_correct(O_tf,M_t(pl_m,:),'n');

%disp('Mean correction, no variance correction directly in the UL PSM...')
M_tbc1=bsxfun(@minus,M_t(pl_m,:),nanmean(M_t(pl_m,:),2));
M_tbc=bsxfun(@plus,M_tbc1,nanmean(O_t(pl_o,io),2));

% Uni-variate linear fit with 'observed' data
% transfer those fits to the model data
b=zeros(2,pNum);
r=zeros(pNum,length(io));
pcr=zeros(pNum,1);
for k=1:pNum
    o_t=O_t(pl_o(k),io);
    pd=prxs(ip,k);
    [b(:,k),~,r(k,:)]=regress(pd(:),[ones(size(o_t(:))) o_t(:)]);
    % Compute proxy--data correlation
    pcr(k)=nancorr(pd(:),o_t(:));
end

%-------------------------------------------
% APPLY LINEAR FIT TO MODEL DATA AND OUTPUT
%-------------------------------------------

for jj=1:pNum
    
    % Obs error estimate (obs error variance)
    proxy{jj}.R=nanmean(r(jj,:).^2).*ones(1,itr); % MSE = evar + bias^2
    proxy{jj}.Rvec=r(jj,:); % save all error terms
    
    proxy{jj}.HXb=b(1,jj) + b(2,jj)*M_tbc(jj,:);
    proxy{jj}.ptype=archv_tp{jj}; % save type
    proxy{jj}.pcr=pcr(jj); % save proxy--data correlation
    
    % Make each age ensemble the same (will be modified later?)
    age_yng=repmat(year',[1 itr]);
    age_old=repmat(year',[1 itr]);
    % no age ensemble
    %age_yng=year';
    %age_old=year';
    
    % Put proxies into structure
    proxy{jj}.age_yng=age_yng;
    proxy{jj}.age_old=age_old;
    proxy{jj}.data=repmat(prxs(:,jj)',[itr,1]); % row vector ensemble
    proxy{jj}.location=[plat(jj) plon(jj)];
    
    %ii=ii+1;
end

%timestp=datestr(now);timestp=strrep(timestp,' ','_');
%save(['debug_psm_ul_' timestp '.mat'])

end

