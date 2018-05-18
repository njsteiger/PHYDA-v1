
% PARAMETERS FOR DATA ASSIMILATION BASED RECONSTRUCTIONS

%rng('shuffle') % Shuffle random number seed

%=================================
% Choose the model for the prior
%=================================

% CESM option: 002 through 010
mdl_prior='cesm_lme010';

disp(['Prior model = ' mdl_prior])

%----------------------------------------
% Base time scale for the reconstruction
%----------------------------------------
% Jan=1, Feb=2, Mar=3, Apr=4, May=5, Jun=6
% Jul=7, Aug=8, Sep=9, Oct=10, Nov=11, Dec=12

% Annual average (respecting the climate system, not the calendar)
%mon_avg_o=4;mon_avg_f=3;
%mon_o='Apr';mon_f='Mar';

% DJF seasonal average
%mon_avg_o=12;mon_avg_f=2; 
%mon_o='Dec';mon_f='Feb'; 

% JJA seasonal average
mon_avg_o=6;mon_avg_f=8; 
mon_o='Jun';mon_f='Aug'; 

disp(['Seasonality = ' mon_o '--' mon_f])

% Enter state vector data types (no spaces)
%  surface (t)emperature, 500hPa (g)eopotential,
%  (p)recip, surface p(r)essure, (e) P-E
%  (a)moc ind., am(o), ocean (h)eat OR N.(h)eat trans, (u)wind200mb
%  (s)ea level pressure or (s)st, (w)ind, soi(l) moisture, (z)onal precip
%  OR itc(z), stream (f)unction, 
%  PDSI (r) [full field], PDSI (q) [land only]
%  SPEI12 (i) [full field], SPEI12 (j) [land only] 
%  PDSI NAS(w) index, PDSI (e)ast Africa index
%  (n)ino 3.4, sea surface salinit(y), gmt (2)m
state_tp = 't2qjozn';
%state_tp = 'tj';
%state_tp = '2ozn';
disp(['State type = ' state_tp ])

% Get model information
[gcm_syr,gcm_eyr]=modelyrs(mdl_prior);

% Number of reconstruction iterations
itr=1;
%disp(['MC iterations = ' num2str(itr)])

% Switchblade proxies? y/n
%swtchbld=30;
%swtchbld=50;
%swtchbld=75;
swtchbld=100; % use all proxies but random ordering
disp(['Switchblade ' num2str(swtchbld) '% of proxies.'])

% Percentiles from posterior
xa_prct=[5 50 95];
disp(['Will save ' num2str(xa_prct) ' percentiles' ])

% Save analysis ensemble?
%ens_save= 'a'; % Save all posterior ensembles: ONLY USE THIS FOR INDEX RECONS!!!
%ens_save= 'y';
ens_save= 'n';
%ens_save='s'; sub_ens=100; % Save subsample of posterior
disp(['Save analysis ensemble? ' ens_save])

% Save output?  y/n
outsave = 'y';
%outsave = 'n';
disp(['Save output? ' outsave])
% Set location of output
outpath='./output-recon/';

if outsave=='y'; disp(['Output will be saved to ' outpath]); end

if strcmp(mdl_prior(1:8),'cesm_lme') 
    
    %rng('default');disp('Using consistent random prior...')
    %rng('shuffle');disp('Using different random prior...')
    %p_yrs=randsample((gcm_syr:gcm_eyr-1)-gcm_syr+1,900); % rand sample
    
    % almost all 1000 years (can't use the beg/end years b/c of year definition and var quirks)
    p_yrs=[851:1848]-gcm_syr+1; 
    
    % Years of the reconstruction and calibration period
    r_o = 1;r_f = 2000; clb_prd=[1920 2000]; % do a verif of 1871:1919 (50 yrs)
    disp(['PSM calib period = ' num2str(clb_prd)]) 
    
    reconYrs=length(r_o:r_f);
else
    error('Must add more options for prior type')
end

disp(['Prior size = ' num2str(length(p_yrs))])

disp(['Reconstruction years = ' num2str(r_o) '-' num2str(r_f) ])

disp('Loading data...')

% LOAD THE MODEL DATA
[Xb_o,xmeta]=load_S(mdl_prior,state_tp,mon_avg_o,mon_avg_f,p_yrs);

% State vector length
svl=size(Xb_o,1);

if ens_save=='a' && svl > 100
   error('State vector is too large to save all posterior ensembles!')
elseif ens_save=='s' && svl*sub_ens*reconYrs*8 > 2e10
   error('State vector will be over 20 GB! Cannot save posterior ensembles.')
end 

% Inflation/deflation factor (if = 1, prior is not inflated)
infl = 1; if infl~=1; disp(['Inflation = ' num2str(infl)]); end

%-----------------------
% PROXY NETWORK OPTIONS
%-----------------------

disp('Loading proxies...')

% Proxy dataset
prxy_dtst=5; % 5 = LMR full proxy database
disp(['Proxy dataset = ' num2str(prxy_dtst)])

% Standardize proxy data? standard (y), (g)aussianize, (n)o
stnd_p='y';
disp(['Standardize proxy data = ' stnd_p])

% Proxy type
%  ice cores (i), corals (c), speleothems (s)
%  tree rings (t), (a)ll
prx_tp='tca'; 
%prx_tp='t'; 
%prx_tp='c';
%prx_tp='a'; 
disp(['Proxy types = ' prx_tp])

% PSM modeling options for trees
% 1 = PSM with just univariate linear T
% 2 = PSM with just univariate linear PDSI
% 3 = pick the PSM with the best correlation among T & PDSI
psm_ot=3;

% PSM modeling options for corals
% 1 = Univariate SST for all
% 2 = Bivariate SST & SSS for all
% 3 = SST & SSS for d18O and SST for everything else
psm_oc=3;

%---------------------
% Load proxy database
%---------------------

% Fully updated to inlude more than pages proxy data
[proxy]=load_prxs_v5(prxy_dtst,stnd_p,prx_tp,itr,...
   mon_avg_o,mon_avg_f,mdl_prior,p_yrs,clb_prd,psm_ot,psm_oc);

disp(['Proxies to assimilate N = ' num2str(length(proxy))])

% file name corresponding to current time stamp
timestp=datestr(now);
timestp=strrep(timestp,' ','_'); % remove space

fln=[mdl_prior '_r' num2str(r_o) num2str(r_f) '_p' num2str(length(p_yrs))...
  '_state_' state_tp '_avg_' mon_o mon_f '_prxydtst_' ...
  num2str(prxy_dtst) '_prxtp_' prx_tp '_' num2str(length(proxy)) ...
  '_swtchbld' num2str(swtchbld) '_' timestp '.mat'];

disp(['outfile = ' fln])












