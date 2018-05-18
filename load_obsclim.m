function [X_clim] = load_obsclim(obs_dtp,mlat,mlon)
%LOAD_OBS_BC  Load observational data climatology for mean state bias correction 
% 	and interpolate the observations to the resolution of the model lat/lon
%	[X_clim] = load_obsclim(obs_dtp,mlat,mlon)
%
%	Inputs: observation data type, climate model lat/lon (for interpolation)
%	Outputs: X obs climatology [mlat,mlon,12 months], using calendar year
%	
%	The possible data source codes are:
% 	'gpcp' (global precipitation)
%	'bearth' (global temperature)
%
%
%   Nathan Steiger, July 2017


switch obs_dtp
    
   case 'gpcp'

      % Load climate data
      xclim=ncread('./input-recon/gpcp_v2.3_climatology.nc','precip');
      xclim=permute(xclim,[2 1 3]);
      xlon = double(ncread('./input-recon/gpcp_v2.3_climatology.nc','lon'));
      xlat = double(ncread('./input-recon/gpcp_v2.3_climatology.nc','lat'));

      % Interpolate to model resolution
      [X_clim] = interpclim(xlat,xlon,xclim,mlat,mlon);


   case 'bearth'
     
      % Load climate data and regrid to climate model format
      xclim=ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      xclim=permute(xclim,[2 1 3]);
      xlon1 = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('./input-recon/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [xclim,xlon]=regrid2gcm(xclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Interpolate to model resolution
      [X_clim] = interpclim(xlat,xlon,xclim,mlat,mlon);


end

