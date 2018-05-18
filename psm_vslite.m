function [proxy] = psm_vslite(prxs,plat,plon,pNum,archv_tp,itr,...
    mdl_prior,p_yrs,clb_prd,year,prxy_dtst)
%PSM_VSLITE  Use VS-Lite as the PSM

addpath('./regem-imputation')

% Create proxy data structure
proxy=struct([]);

% CALIBRATE VS-LITE: observational data available from 1901 through 2014
%clb_prd=[1911 1995]; % for PAGES2k temperature recons
addpath('./vs-lite')

% Calibration parameters of the VS-lite model based on observations
%[T1,T2,M1,M2,~,~,MSE] = tree_vslite_calib(clb_prd,prxs,plat,plon,year,pNum);
%[T1,T2,M1,M2,~,~,MSE,~,vsl_cr] = tree_vslite_calib(clb_prd,prxs,plat,plon,year,pNum);

if prxy_dtst==2
    
    %save(['vs_lite_historical_params_calib_' num2str(clb_prd(1)) '_'...
    %    num2str(clb_prd(2)) '_v2_pages2kv13_1_v4.mat'],'T1','T2','M1','M2','MSE','vsl_cr')
    %load('vs_lite_historical_params_calib_1911_1995_v2_pages2kv13_1_v4.mat','T1','T2','M1','M2','MSE','vsl_cr')
    load(['vs_lite_historical_params_calib_' num2str(clb_prd(1)) '_'...
        num2str(clb_prd(2)) '_v2_pages2kv13_1_v4.mat'],'T1','T2','M1','M2','MSE','vsl_cr')
    
elseif prxy_dtst==4
    
    %save('vs_lite_historical_params_calib_1911_1995_v2_pages2kv13_1_fulldata_v4.mat','T1','T2','M1','M2','MSE','vsl_cr')
    load('vs_lite_historical_params_calib_1911_1995_v2_pages2kv13_1_fulldata_v4.mat','T1','T2','M1','M2','MSE','vsl_cr')
    
else
    error('Unsupported proxy database...')
end



% Check that number of parameters match the number of proxies
if length(T1)~=pNum || length(T1)~=size(prxs,2)
    error('Tree proxies do not match parameter values')
end



% LOAD PRIOR DATA FOR RUNNING VS-LITE
p1=bsxfun(@plus,repmat(p_yrs*12,[12 1]),(0:11)'); % get indices for months
p_yrsmon=p1(:);
[M_t,mlat,mlon]=load_S(mdl_prior,'t',[],[],p_yrsmon); % monthly temperature
[M_p]=load_S(mdl_prior,'p',[],[],p_yrsmon); % monthly precip

% Get the proxy locations in model space
[pl_m] = nearest_latlon(mlat,mlon,plat,plon);

%---------------------------
% BIAS CORRECT MODEL OUTPUT
%---------------------------

% Load monthly observational data
obs_dtp1='cru_ts3_tmp';% temperature
[Od_t,xlat,xlon,x_yrs] = load_obs_S(obs_dtp1,[],[]);

obs_dtp1='cru_ts3_pre';% precipitation
[Od_p] = load_obs_S(obs_dtp1,[],[]);

% Get locations of proxies in obs
mask=ones(size(Od_t(:,end))); % mask out ocean lat/lon values
mask(isnan(Od_t(:,end)))=NaN;
[pl_o] = nearest_latlon(xlat,xlon,plat,plon,mask);

% Bias-correct temperature
% infill random missing values in obs data with regem
disp('---Using 1950:1999 for bias correction---')
c2_yrs=1950:1999;[~,io2] = intersect(x_yrs,c2_yrs);

options.maxit=10; options.regress='iridge'; %options.regress='ttls';
%Od_tf=regem(Od_t(pl_o,:),options); % bias correct based on all data
Od_tf=regem(Od_t(pl_o,io2),options); % bias correct on recent data
[M_tbc] = bias_correct(Od_tf,M_t(pl_m,:),'n');

% Bias-correct precipitation
%odp=Od_p(pl_o,:); % all precip data
odp=Od_p(pl_o,io2); % just modern data
mdp=double(M_p(pl_m,:));
odp(odp==0)=0.01; % replace zero exactly with value less than error in obs
mdp(mdp<1e-13)=1e-13; % replace a dozen or so really low values
% infill missing values with regem
options.maxit=10; options.regress='iridge'; %options.regress='ttls';
odp_f=regem(odp,options);
odp_f(odp_f<0)=0.01; % make sure there aren't any zeros
[M_pbc] = bias_correct(odp_f,mdp,'g');

%-------------
% RUN VS-LITE
%-------------

% Reshape data
k=1;l=12;
t_cal=zeros(pNum,12,length(p_yrs));
p_cal=zeros(pNum,12,length(p_yrs));
for n=1:length(p_yrs)
    t_cal(:,:,n)=M_tbc(:,k:l);
    p_cal(:,:,n)=M_pbc(:,k:l);
    k=k+12; l=l+12;
end

M_trw=zeros(pNum,length(p_yrs));
%gT=zeros(pNum,12,length(p_yrs)); % growth from T and M
%gM=zeros(pNum,12,length(p_yrs));%,gT(q,:,:),gM(q,:,:)
for q=1:pNum
    [M_trw(q,:)] = VSLite_v2_3(p_yrs(1),p_yrs(end),plat(q),T1(q),T2(q),...
        M1(q),M2(q),squeeze(t_cal(q,:,:)),squeeze(p_cal(q,:,:)),...
        'lbparams',[0.76,0.01,0.093,4.886,5.80,1000,0.2,1]);
end

% Throw out locations that are too cold to model
tr_nanids=find(isnan(M_trw(:,2))); % 2nd year of all

if ~isempty(tr_nanids)
    disp([num2str(length(tr_nanids)) ' tree rings'...
        ' could not initiate growth...'])
       
    prxs(:,tr_nanids)=NaN*prxs(:,tr_nanids);
    MSE(tr_nanids)=NaN*MSE(tr_nanids);
    
    %disp('Saving removed tree locations...')
    %save('M_trw_nanids.mat','tr_nanids')
end

%-----------------------
% ADD PROXIES TO OUTPUT
%-----------------------

for jj=1:pNum
    
    % Obs error estimate (obs error variance)
    %proxy{ii}.R=eVar(jj).*ones(1,itr); % MSE = evar + bias^2, if bias = 0 => MSE = evar
    proxy{jj}.R=MSE(jj).*ones(1,itr);
    
    proxy{jj}.HXb=M_trw(jj,:); % VS-lite TRW estimates
    proxy{jj}.ptype=archv_tp{jj}; % save type
    proxy{jj}.pcr=vsl_cr(jj); % save proxy--data correlation
    
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

%save('vsl_output.mat','-v7.3')

end

