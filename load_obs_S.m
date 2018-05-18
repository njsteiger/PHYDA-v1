function [X_obs,xlat,xlon,x_yrs] = load_obs_S(obs_dtp,mon_avg_o,mon_avg_f)
%LOAD_S  Load observational data for statistical PSMs
%   [X,xlat,xlon,id_X] = load_obs_S(model,state_tp,mon_avg_o,mon_avg_f,s_yrs)
%     Given a model code 'obs_dtp', state vector variable type 'state_tp',
%     start and end months for averaging, and the years from the model
%     'm_yrs', this function will load the data saved to the computer
%     'humid'. The outputs are the state vector 'X', the latitudes and
%     longitudes of the gridded data, and the locations of the state
%     variables 'id_X' which are in the same order as given in 'state_tp'.
%     The state variables are annually averaged according to the start
%     month and end moth inputs, unless these are empty in which case the
%     original monthly data is returned.
%     The possible data source codes are:
%     'cru_ts3_tmp'
%     'hadcrut4'
%     'bearth'
%     'mlost'
%     'cru_ts3_pre'
%     'soda_sst'
%     'soda_sss'
%     'my_spei'
%     
%
%   Nathan Steiger, Sept 2016


switch obs_dtp
    
    case 'cru_ts3_tmp'
        
        %------------------------
        % TEMPERATURE, LAND ONLY
        %------------------------
        
        X=ncread('./input-recon/cru_ts3.23.1901.2014.tmp.dat.nc','tmp');
        X=permute(X,[2 1 3]);
        xlon = double(ncread('./input-recon/cru_ts3.23.1901.2014.tmp.dat.nc','lon'));
        xlat = double(ncread('./input-recon/cru_ts3.23.1901.2014.tmp.dat.nc','lat'));
        [X,xlon]=regrid2gcm(X,xlon);
        xlat=xlat(:);xlon=xlon(:);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        x_yrs=1901:(1901+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        % Remove the time mean to be consistent with prior
        %X_obs=X_obs-repmat(nanmean(X_obs,2),[1,size(X_obs,2)]);
        
    case 'hadcrut4'
        
        
        % contains both temperature and lats & lons
        load('./input-recon/HadCRUT4.3_GraphEM_SP80_18502014.mat')
        
        
        tas_inf=mon2ann(tas,mon_avg_o,mon_avg_f);
        tas_inf(tas_inf==0)=NaN; % make fill values nans
        
        x_yrs=1850:(1850+size(tas_inf,3)-1);
        
        
        X_obs=reshape(tas_inf,size(tas_inf,1)*size(tas_inf,2),size(tas_inf,3));
        
        xlat=lat(:);xlon=lon(:);
        
    case 'bearth'
        
        
        tas=ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','temperature');
        tas=tas(:,:,1:1992); % just take from 1850 to 2015
        tas=permute(tas,[2 1 3]);
        lon = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
        lat = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
        [tas,lon]=regrid2gcm(tas,lon);
        xlat=lat(:);xlon=lon(:);
        
        tas_inf=mon2ann(tas,mon_avg_o,mon_avg_f);
        tas_inf(tas_inf==0)=NaN; % make fill values nans
        
        %x_yrs=1880:(1880+size(tas_inf,3)-1); % for infilled data
        x_yrs=1850:(1850+size(tas_inf,3)-1);
        
        X_obs=reshape(tas_inf,size(tas_inf,1)*size(tas_inf,2),size(tas_inf,3));
        
        
    case 'cru_ts3_pre'
        
        X=ncread('./input-recon/cru_ts3.23.1901.2014.pre.dat.nc','pre');
        X=permute(X,[2 1 3]);
        xlon = double(ncread('./input-recon/cru_ts3.23.1901.2014.pre.dat.nc','lon'));
        xlat = double(ncread('./input-recon/cru_ts3.23.1901.2014.pre.dat.nc','lat'));
        [X,xlon]=regrid2gcm(X,xlon);
        xlat=xlat(:);xlon=xlon(:);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        x_yrs=1901:(1901+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        
    case 'soda_sst'
        
        X=squeeze(ncread('./input-recon/SODAsi_mon_1871_2011.nc','TEMP_MN'));
        X=permute(X,[2 1 3]);
        xlon = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LON'));
        xlat = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LAT'));
        xlat=xlat(:);xlon=xlon(:);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        x_yrs=1871:(1871+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        
    case 'soda_sss'
        
        X=squeeze(ncread('./input-recon/SODAsi_mon_1871_2011.nc','SALT_MN'));
        X=permute(X,[2 1 3]);
        xlon = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LON'));
        xlat = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LAT'));
        xlat=xlat(:);xlon=xlon(:);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        x_yrs=1871:(1871+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        
    case 'dai_pdsi'
        
        X=ncread('./input-recon/pdsi_mon_mean_selfcalibrated_dai.nc','pdsi');
        X(X<-100)=NaN; % replace fill values
        X=permute(X,[2 1 3]);
        xlon = double(ncread('./input-recon/pdsi_mon_mean_selfcalibrated_dai.nc','lon'));
        xlat = double(ncread('./input-recon/pdsi_mon_mean_selfcalibrated_dai.nc','lat'));
        [X,xlon]=regrid2gcm(X,xlon);
        xlat=xlat(:);xlon=xlon(:);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        x_yrs=1850:(1850+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        
    case 'my_pdsi'
        
        error('Not available yet...')
        
    case 'my_spei'
        
        load('./input-recon/spei_cruts3_scl_12_krnl_e_02-Sep-2016_16:54:56.mat');
        X=spei_f;
        xlon=lon;
        xlat=lat;
        [X,xlon]=regrid2gcm(X,xlon);
        
        Xann=mon2ann(X,mon_avg_o,mon_avg_f);
        Xann(Xann==0)=NaN; % add back fill nans
        
        %x_yrs=syr:eyr;
        x_yrs=syr:(syr+size(Xann,3)-1);
        
        X_obs=reshape(Xann,size(Xann,1)*size(Xann,2),size(Xann,3));
        
        
end

