function [proxy] = psm_bl(prxs,plat,plon,pNum,archv_tp,itr,...
    mdl_prior,p_yrs,clb_prd,year,mon_avg_o,mon_avg_f,obs_dtp1,obs_dtp2)
%PSM_BL  Bi-variate linear PSM for coral d18O proxies

addpath('./regem-imputation')

% Create proxy data structure
proxy=struct([]);

% Throw out proxy data that doesn't overlap instrumental period
% with at least 20 points
%clb_prd=[1911 1995];
c_yrs=clb_prd(1):clb_prd(2);
[~,ip] = intersect(year,c_yrs); % overlap years of proxy data
% Find indices to remove because of too few data points for accurate
% regressions
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
%obs_dtp1='soda_sst';% sea surface temperature
%[O_sst,xlat,xlon,x_yrs] = load_obs_S(obs_dtp1,mon_avg_o,mon_avg_f);
% LOAD *MONTHLY* OBSERVATIONAL DATA
[O_sst,xlat,xlon,x_yrs] = load_obs_Smon(obs_dtp1);

%obs_dtp2='soda_sss';% sea surface salinity
%[O_sss] = load_obs_S(obs_dtp2,mon_avg_o,mon_avg_f);
% LOAD *MONTHLY* OBSERVATIONAL DATA
[O_sss] = load_obs_Smon(obs_dtp2);

% Just use contemporary period for bias-correction 
disp('---Using 1950:1999 for PSM bias correction---')
c2_yrs=1950:1999;[~,io2] = intersect(x_yrs,c2_yrs);
% Using indices of years above, deduce months to extract
iom=bsxfun(@plus,(io2'-1)*12+1,repmat(0:11,[length(io2),1])');
iom=iom(:);

% Get locations of proxies in obs
% mask out land lat/lon values
mask=sum(~isnan(O_sst),2);
mask(mask==0)=NaN;
mask(~isnan(mask))=1;
[pl_o] = nearest_latlon(xlat,xlon,plat,plon,mask);

% Pull out obs data at proxy locations
O_sstpmon=O_sst(pl_o,iom);
O_ssspmon=O_sss(pl_o,iom);
clear O_sst O_sss

% LOAD MODEL DATA
%[M_sst,~,~,~,olat,olon]=load_S(mdl_prior,'s',mon_avg_o,mon_avg_f,p_yrs);
%[M_sss]=load_S(mdl_prior,'y',mon_avg_o,mon_avg_f,p_yrs);

% LOAD *MONTHLY* MODEL DATA
[M_sst,~,~,~,olat,olon] = load_Smon(mdl_prior,'s',mon_avg_o,mon_avg_f,p_yrs);
[M_sss] = load_Smon(mdl_prior,'y',mon_avg_o,mon_avg_f,p_yrs);


% Locations of proxies in the model
% mask out land lat/lon values
mask=sum(~isnan(M_sst),2);
mask(mask==0)=NaN;
mask(~isnan(mask))=1;

[pl_m] = nearest_latlon(olat,olon,plat,plon,mask);

% Pull out model data at proxy locations
M_sstpmon=M_sst(pl_m,:);
M_ssspmon=M_sss(pl_m,:);
clear M_sst M_sss

%--------------------------------------------------
% Bias-correct model output at proxy locations only
%--------------------------------------------------

% BIAS-CORRECT THE SEASONAL CYCLE (use this approach because SODA doesn't provide it's own climatology file)
[M_sstbc_mon] = bias_cr_seas(O_sstpmon,M_sstpmon);
[M_sssbc_mon] = bias_cr_seas(O_ssspmon,M_ssspmon);

% Take the annual mean
if max(diff(p_yrs)>1);error('Annualizing SSS and SST data only applies for continuous, non-random-in-time prior');end
M_sstbc0=mon2ann(M_sstbc_mon,mon_avg_o,mon_avg_f);
M_sssbc0=mon2ann(M_sssbc_mon,mon_avg_o,mon_avg_f);

% LOAD ANNUAL OBSERVATIONAL DATA
[O_sst] = load_obs_S(obs_dtp1,mon_avg_o,mon_avg_f);
[O_sss] = load_obs_S(obs_dtp2,mon_avg_o,mon_avg_f);

% BIAS-CORRECT THE ANNUAL VARIABILITY
%[M_sstbc] = bias_correct(O_sst(pl_o,:),M_sst(pl_m,:),'n');
%[M_sssbc] = bias_correct(O_sss(pl_o,:),M_sss(pl_m,:),'n');

% Years of overlap
[~,io] = intersect(x_yrs,c_yrs);

% disp('Mean correction, no annual variance correction directly in the BL PSM...')
M_sstbc1=bsxfun(@minus,M_sstbc0,nanmean(M_sstbc0,2));
M_sstbc=bsxfun(@plus,M_sstbc1,nanmean(O_sst(pl_o,io),2));
M_sssbc1=bsxfun(@minus,M_sssbc0,nanmean(M_sssbc0,2));
M_sssbc=bsxfun(@plus,M_sssbc1,nanmean(O_sss(pl_o,io),2));



% Bi-variate linear fit with 'observed' SST and SSS and
% transfer those fits to the model data
b=zeros(3,pNum);
r=zeros(pNum,length(io));
for k=1:pNum
    o_sst=O_sst(pl_o(k),io);
    o_sss=O_sss(pl_o(k),io);
    pd=prxs(ip,k);
    [b(:,k),~,r(k,:)]=regress(pd(:),[ones(size(o_sst(:))) o_sst(:) o_sss(:)]);
end

%save('pages2kv2_corals_b.mat','b','r','plat','plon')


%-------------------------------------------
% APPLY LINEAR FIT TO MODEL DATA AND OUTPUT
%-------------------------------------------

for jj=1:pNum
    
    % Obs error estimate (obs error variance)
    proxy{jj}.R=nanmean(r(jj,:).^2).*ones(1,itr); % MSE = evar + bias^2
    proxy{jj}.Rvec=r(jj,:); % save all error terms    

    proxy{jj}.HXb=b(1,jj) + b(2,jj)*M_sstbc(jj,:) + b(3,jj)*M_sssbc(jj,:);
    proxy{jj}.ptype=archv_tp{jj}; % save type
    proxy{jj}.b=b(:,jj); % save type
    
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





end

