function [proxy]=load_prxs_v5(prxy_dtst,stnd_p,prx_tp,itr,...
    mon_avg_o,mon_avg_f,mdl_prior,p_yrs,clb_prd,psm_ot,psm_oc)
% LOAD_PRXS Loads proxies for reconstructions as well
% as computes the proxy forward models and errors.

addpath('./regem-imputation')

% Create proxy data structure
proxy=cell(0);

year=[]; % matlab has a 'year' function, so define variable here

% LOAD PROXY DATA
if prxy_dtst==1
    
    % PAGES2K TEMPERATURE RECONSTRUCTION DATABASE --- screened!!!!
    
    load(['./proxy-data/PAGES2k_v2/v2_0_0/proxy_ama_2.0.0_calib-selection_1881_1916_1995_0.67_infilled_DINEOF_PAGES-crit-regional+FDR.mat'])
    disp('Using SCREENED PAGES2k temperature proxy data v 2.0.0')
    
    % Convert proxy lons into 0 to 360
    nlns=find(p_lon<0);
    p_lon(nlns)=p_lon(nlns)+360;
    
    if mon_avg_o==1 && mon_avg_f==12
        error('Unsupported proxy averages...')
        %prx_i=proxy_ann;
    elseif mon_avg_o==4 && mon_avg_f==3
        prx_i=proxy_ama; % Use Apr-Mar averaged proxies
    else
        error('Unsupported proxy averages...')
    end

elseif prxy_dtst==2
    
    % PAGES2K TEMPERATURE RECONSTRUCTION DATABASE
    
    % proxy_ama_2.0.0_highres.mat % ALL HI-RES PROXIES
    load(['./proxy-data/PAGES2k_v2/v2_0_0/proxy_ama_2.0.0_highres.mat'])
    disp('Using full-selected PAGES2k temperature proxy data v 2.0.0')
    
    % Convert proxy lons into 0 to 360
    nlns=find(p_lon<0);
    p_lon(nlns)=p_lon(nlns)+360;
    
    if mon_avg_o==1 && mon_avg_f==12
        error('Unsupported proxy averages...')
        %prx_i=proxy_ann;
    elseif mon_avg_o==4 && mon_avg_f==3
        prx_i=proxy_ama; % Use Apr-Mar averaged proxies
    else
        error('Unsupported proxy averages...')
    end
    
elseif prxy_dtst==4
    
    % PAGES2K FULL RECONSTRUCTION DATABASE
    
    load(['./proxy-data/PAGES2k_v2/v1_13_1/'...
        'pages2kTSv1_13_1_unpack_AMA_allproxies_njs_25-Jan-2017_14:12:28.mat'])
    disp('Using PAGES2k temperature proxy data v 1.13.1')
    
    % Convert proxy lons into 0 to 360
    nlns=find(p_lon<0);
    p_lon(nlns)=p_lon(nlns)+360;
    
    if mon_avg_o==1 && mon_avg_f==12
        error('Unsupported proxy averages...')
        %prx_i=proxy_ann;
    elseif mon_avg_o==4 && mon_avg_f==3
        prx_i=proxy_ama; % Use Apr-Mar averaged proxies
    else
        error('Unsupported proxy averages...')
    end
    
    
elseif prxy_dtst==5
    
    % LMR FULL RECONSTRUCTION DATABASE
    
    if mon_avg_o==1 && mon_avg_f==12
        
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_v0.2.0.mat'])
        disp('Using full LMR proxy data v 0.2.0 using the calendar year')
        
    elseif mon_avg_o==4 && mon_avg_f==3
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_aprmar_v0.2.0.mat'])
        disp('Using full LMR proxy data v 0.2.0 using the "tropical" year')
        
    else
        % This option kicks in for seasonal reconstructions: just use Apr - Mar coral proxy averaging
        disp('No annual averaging chosen for proxies, so using the default full LMR proxy data v 0.2.0 with the "tropical" year')
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_aprmar_v0.2.0.mat'])
    end
    
    
    % dims = [time,proxy sites]
    prx_i=lmr2k_data;
    
    year=year(:)'; % needs to be row vector
elseif prxy_dtst==3
    
    % Select out PAGES2k proxies from LMR database: this will allow me to
    % get all the metadata from the LMR database that's missing for PAGES2k
    
    if mon_avg_o==1 && mon_avg_f==12
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_v0.2.0.mat'])
    elseif mon_avg_o==4 && mon_avg_f==3
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_aprmar_v0.2.0.mat'])
    else 
        disp('No annual averaging chosen for proxies, so using the default "tropical" year')
        load(['./proxy-data/LMR_NCDC_database_v.0.2.0/'...
            'lmr2k_proxydata_aprmar_v0.2.0.mat'])
    end
    
    disp('Using full PAGES2k data via LMR database')
    aa=regexp(lmr2k_names,'^PAGES2kv2\w*');
    p2k_ids=find(~cellfun(@isempty,aa));
    
    
    % dims = [time,proxy sites]
    prx_i=lmr2k_data(:,p2k_ids);
    p_lat=p_lat(p2k_ids);
    p_lon=p_lon(p2k_ids);
    archive=archive(p2k_ids);
    
    year=year(:)'; % needs to be row vector
else
    error('Unsupported proxy dataset type...')
end



% Standardize or normalize proxy data?
if stnd_p=='y'
    % standardize proxy data
    prx_im=bsxfun(@minus,prx_i,nanmean(prx_i,1));
    prx_d=bsxfun(@rdivide,prx_im,nanstd(prx_im,0,1));
elseif stnd_p=='g'
    % Impose normality on proxy data
    prx_d=NaN(size(prx_i));
    for i=1:size(prx_i,2)
        prx_d(:,i)=gaussianize(prx_i(:,i));% WORKS FOR VECTORS BUT NOT MATRICES...
    end
else
    prx_d=prx_i;
end



%----------------------------------------------
% Load and estimate each proxy type at a time
%----------------------------------------------



%ii=1;
for j=1:length(prx_tp)
    switch prx_tp(j)
        case 't' % trees
            % SELECT PROXY DATA
            p_idx = find(ismember(archive,{'tree','Tree Rings'})); % 'ice', 'coral', 'documents', 'glacier ice'
            
            prxs=prx_d(:,p_idx);
            plat=p_lat(p_idx);
            plon=p_lon(p_idx);
            pNum=length(p_idx);
            archv_tp=archive(p_idx);
            
            if psm_ot==4 % = pick the PSM with the best correlation among T, PDSI, SPEI, VSL
                error('VS-lite option needs updating...')
                disp('%%%%%% TREE PSMs USING T, PDSI, SPEI, VS-LITE %%%%%%')
                
                % VS-Lite PSM
                disp('Calculating VS-Lite PSM...')
                [proxy_vsl] = psm_vslite(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,prxy_dtst);
                
                % PDSI PSM
                disp('Calculating PDSI PSM...')
                obs_dtp1='dai_pdsi';mdl_dtp='r';
                %mnavg_p_o=6;mnavg_p_f=8; % NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_pdsi] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                % SPEI PSM
                disp('Calculating SPEI PSM...')
                obs_dtp1='my_spei';mdl_dtp='i';
                %mnavg_p_o=6;mnavg_p_f=8;% NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_spei] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                % Temperature PSM
                disp('Calculating T PSM...')
                obs_dtp1='bearth';mdl_dtp='t';
                %mnavg_p_o=6;mnavg_p_f=8;% NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_t] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                mxcr=zeros(4,pNum);
                for kk=1:pNum
                    
                    %  Find maximum correlation
                    mxcr(1,kk)=proxy_vsl{kk}.pcr;
                    mxcr(2,kk)=proxy_pdsi{kk}.pcr;
                    mxcr(3,kk)=proxy_spei{kk}.pcr;
                    mxcr(4,kk)=proxy_t{kk}.pcr;
                    
                    [~,mi]=max(abs(mxcr(:,kk)));
                    
                    switch mi
                        case 1
                            proxy_i{kk}=proxy_vsl{kk};
                            proxy_i{kk}.psm_mthd='v';
                        case 2
                            proxy_i{kk}=proxy_pdsi{kk};
                            proxy_i{kk}.psm_mthd='p';
                        case 3
                            proxy_i{kk}=proxy_spei{kk};
                            proxy_i{kk}.psm_mthd='s';
                        case 4
                            proxy_i{kk}=proxy_t{kk};
                            proxy_i{kk}.psm_mthd='t';
                        otherwise
                            error('Wrong minimum value...')
                    end
                    
                end
                
                % Concatenate proxy types together
                proxy=[proxy,proxy_i];
                
                
            elseif psm_ot==3 % = pick the PSM with the best correlation among T, PDSI
                
                disp('%%%%%% TREE PSMs USING ONLY PDSI AND T %%%%%%')
                
                % PDSI PSM
                disp('Calculating PDSI PSM...')
                obs_dtp1='dai_pdsi';mdl_dtp='r';
                %mnavg_p_o=6;mnavg_p_f=8; % NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_pdsi] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                % Temperature PSM
                disp('Calculating T PSM...')
                obs_dtp1='bearth';mdl_dtp='t';
                %mnavg_p_o=6;mnavg_p_f=8;% NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_t] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                mxcr=zeros(2,pNum);
                for kk=1:pNum
                    
                    mxcr(1,kk)=proxy_pdsi{kk}.pcr;
                    mxcr(2,kk)=proxy_t{kk}.pcr;
                    
                    [~,mi]=max(abs(mxcr(:,kk)));
                    
                    switch mi
                        case 1
                            proxy_i{kk}=proxy_pdsi{kk};
                            proxy_i{kk}.psm_mthd='p';
                        case 2
                            proxy_i{kk}=proxy_t{kk};
                            proxy_i{kk}.psm_mthd='t';
                        otherwise
                            error('Wrong minimum value...')
                    end
                    
                end
                
                % Concatenate proxy types together
                proxy=[proxy,proxy_i];
                
            elseif psm_ot==1 % PSM with just linear T
                
                disp('%%%%%% TREE PSMs USING ONLY T %%%%%%')
                
                % Temperature PSM
                %disp('Calculating T PSM...')
                obs_dtp1='bearth';mdl_dtp='t';disp('Using BEarth')
                %obs_dtp1='hadcrut4';mdl_dtp='t';disp('Using Hadcrut4')
                %mnavg_p_o=6;mnavg_p_f=8;% NH growing season
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_t] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                % Concatenate proxy types together
                proxy=[proxy,proxy_t];
                
            elseif psm_ot==2 % PSM with just linear PDSI
                
                disp('%%%%%% TREE PSMs USING ONLY PDSI %%%%%%')
                
                % PSM
                obs_dtp1='dai_pdsi';mdl_dtp='r';
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
                
                [proxy_t] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                    obs_dtp1,mdl_dtp);
                
                % Concatenate proxy types together
                proxy=[proxy,proxy_t];
            else
                error('Wrong PSM type for tree rings!')
            end
            
        case 'c' % corals
            
            
            % SELECT PROXY DATA
            p_idx = find(ismember(archive,{'coral','sclerosponge','Corals and Sclerosponges','Bivalve'}));
            
            prxs=prx_d(:,p_idx);
            plat=p_lat(p_idx);
            plon=p_lon(p_idx);
            pNum=length(p_idx);
            archv_tp=archive(p_idx);
            
            %-------------------------------------------------------------
            % Do bivariate fit for d18O, everything else a univariate fit
            %-------------------------------------------------------------
            
            % option for both bivariate and univariate
            if psm_oc==3
                
                disp('%%%%%% CORALS: BIVARIATE AND UNIVARITE %%%%%%')
		% Find corals that are d18O
                aa=regexp(msrmt(p_idx),'^d18O\w*');
                d18o_ids=find(~cellfun(@isempty,aa));
                
                % Bivariate linear coral PSM
                obs_dtp1='soda_sst';% sea surface temperature
                obs_dtp2='soda_sss';% sea surface salinity
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f; % use same seasonal averaging as prior
                [proxy_b] = psm_bl(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,obs_dtp1,obs_dtp2);
                
                
                % Univariate coral temperature PSM
                %obs_dtp1='soda_sst';% sea surface temperature
                %mdl_dtp='s'; % SST in model
                obs_dtp1='bearth';mdl_dtp='t';disp('Using BEarth')
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f; % use same seasonal averaging as prior
                [proxy_u] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,obs_dtp1,mdl_dtp);
                
                
                % Initialize with temperature only, replace d18O proxies
                proxy_i=proxy_u;
                for oo=1:length(d18o_ids);
                    proxy_i{d18o_ids(oo)}=proxy_b{d18o_ids(oo)};
                    proxy_i{d18o_ids(oo)}.psm_mthd='sst_sss';
                end
                
            elseif psm_oc==1
                
                disp('%%%%%% Assuming univariate SST fit for all corals %%%%%%')
                
                % Univariate coral temperature PSM
                %obs_dtp1='soda_sst';mdl_dtp='s'; % SST in model
                %obs_dtp1='hadcrut4';mdl_dtp='t';disp('Using Hadcrut4')
                obs_dtp1='bearth';mdl_dtp='t';disp('Using BEarth')
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f; % use same seasonal averaging as prior
                [proxy_i] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,obs_dtp1,mdl_dtp);
                
                
            elseif psm_oc==2
                
                disp('%%%%%% Assuming bivariate SST & SSS fit for all corals %%%%%%')
                
		% Bivariate linear coral PSM
                obs_dtp1='soda_sst';% sea surface temperature
                obs_dtp2='soda_sss';% sea surface salinity
                mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f; % use same seasonal averaging as prior
                [proxy_i] = psm_bl(prxs,plat,plon,pNum,archv_tp,itr,...
                    mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,obs_dtp1,obs_dtp2);
                
                
            else
                error('Wrong coral proxy PSM choice...')
            end
            
            % Concatenate proxy types together
            proxy=[proxy,proxy_i];
            
            %save('coral_output2.mat')
            
            
        case 'a' % all others
            
            %=======================
            % NOT YET ACCOUNTED FOR
            % case 'i' % ice cores
            % {'glacier ice'}
            %=======================
            
            disp('%%%%%% ALL OTHER PROXY TYPES: UNIVARITE WITH TEMPERATURE %%%%%%')
            
            % Do proxies that AREN'T the following:
            p_idx = find(~ismember(archive,{'coral','tree','sclerosponge',...
                'marine sediment','Corals and Sclerosponges','Tree Rings','Bivalve'}));
            %disp('Throwing out marine sediment proxies...')
            %p_idx = find(ismember(archive,{'marine sediment'}));
            
            
            % Testing univariate linear fit for other proxy types
            %othr_idx = find(~ismember(archive,{''}));disp('All proxies are UL fit')
            %othr_idx = find(ismember(archive,{'tree'}));
            %othr_idx = find(ismember(archive,{'coral','sclerosponge'}));
            
            prxs=prx_d(:,p_idx);
            plat=p_lat(p_idx);
            plon=p_lon(p_idx);
            pNum=length(p_idx);
            archv_tp=archive(p_idx);
            
            
            % Throw out extreme outlier proxies
            maxp=max(abs(prxs));
            outl=find(maxp > 7.5); % throw out if greater than 7.5 std units
            prxs(:,outl)=[];
            plat(outl)=[];
            plon(outl)=[];
            pNum=length(plon);
            archv_tp(outl)=[];
            
            % Temperature PSM
            obs_dtp1='bearth';mdl_dtp='t'; disp('Using BEarth')
            %obs_dtp1='hadcrut4';mdl_dtp='t';disp('Using Hadcrut4')
            %mnavg_p_o=6;mnavg_p_f=8;% NH growing season
            mnavg_p_o=mon_avg_o;mnavg_p_f=mon_avg_f;
            
            [proxy_i] = psm_ul(prxs,plat,plon,pNum,archv_tp,itr,...
                mdl_prior,p_yrs,clb_prd,year,mnavg_p_o,mnavg_p_f,...
                obs_dtp1,mdl_dtp);
            
            % Concatenate proxy types together
            proxy=[proxy,proxy_i];
            
    end
end

% Initial proxies
disp(['Initially ' num2str(length(proxy)) ' proxies'])

ii=1;
% Pull out proxies that are just NaNs (from PSM fits)
for i=1:length(proxy)
    if isnan(proxy{i}.R(1)) || sum(~isnan(proxy{i}.data(1,:)))==0 || sum(proxy{i}.HXb)==0
        rmidx(ii)=i;
        ii=ii+1;
    elseif proxy{i}.R(1) < 1e-3 % remove if R is unrealistically small (bad regression)
        rmidx(ii)=i;
        ii=ii+1;
    elseif sum(isnan(proxy{i}.HXb(1,:)))>0 % sometimes sporadic NaNs in HXb
        rmidx(ii)=i;
        ii=ii+1;
    end
end

if ii>1
    proxy(rmidx)=[];
    disp([num2str(ii-1) ' proxies removed because of PSM problems' ]);
end

end


