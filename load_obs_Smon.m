function [X_obs,xlat,xlon,x_yrs] = load_obs_Smon(obs_dtp)
%LOAD_OBS_SMON  Load observational data for statistical PSMs, *MONTHLY*
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
%   Nathan Steiger, LDEO July 2017; Dec 2017


switch obs_dtp
    
        
    case 'bearth'
        
        
        tas=ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','temperature');
        tas=tas(:,:,1:1992); % just take from 1850 to 2015
        tas=permute(tas,[2 1 3]);
        lon = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
        lat = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
        [Xv_mon,lon]=regrid2gcm(tas,lon);
        xlat=lat(:);xlon=lon(:);
        
        x_yrs=1850:(1850+size(Xv_mon,3)/12-1);
        
        X_obs=reshape(Xv_mon,size(Xv_mon,1)*size(Xv_mon,2),size(Xv_mon,3));
        
        
    case 'soda_sst'
        
        Xv_mon=squeeze(ncread('./input-recon/SODAsi_mon_1871_2011.nc','TEMP_MN'));
       	Xv_mon=permute(Xv_mon,[2 1 3]);
        x_yrs=1871:(1871+size(Xv_mon,3)/12-1);
        
        xlon = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LON'));
        xlat = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LAT'));
        xlat=xlat(:);xlon=xlon(:);
        
        X_obs=reshape(Xv_mon,size(Xv_mon,1)*size(Xv_mon,2),size(Xv_mon,3));
        
    case 'soda_sss'
        
        Xv_mon=squeeze(ncread('./input-recon/SODAsi_mon_1871_2011.nc','SALT_MN'));
       	Xv_mon=permute(Xv_mon,[2 1 3]);
        x_yrs=1871:(1871+size(Xv_mon,3)/12-1);
        
        xlon = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LON'));
        xlat = double(ncread('./input-recon/SODAsi_mon_1871_2011.nc','LAT'));
        xlat=xlat(:);xlon=xlon(:);
        
        X_obs=reshape(Xv_mon,size(Xv_mon,1)*size(Xv_mon,2),size(Xv_mon,3));
        
        
        
end

