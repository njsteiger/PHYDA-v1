function [itcz_mx,itcz_cm,imeta] = itcz_lon(pr,lat,lon,lonext)
%ITCZ Compute the latitude of the ITCZ over certain longitudinal range
%   [itcz_mx,itcz_cm] = itcz_lon(pr,lat,lon,lonext)
%       Input:  pr = annual pr field (lat,lon,time)
%               lat/lon = (code currently assumes regular lat/lon grid)
%               lonext = longitudinal extent (two element vector), left to
%               right bounds (e.g., lonext=[343 52] for Africa)
%       Output: itcz_mx = ITCZ expectation maximum latitude
%		itcz_cm = ITCZ center of mass latitude
%               imeta = metadata of ITCZ variable and also the two above variables
%
% ITCZ_max = latitude of expected latitudes using a weighting function of 
% an integer power N=10 of the area-weighted precipitation P integrated over
% the tropics. See Adam et al 2016: https://doi.org/10.1175/JCLI-D-15-0512.1
%
% ITCZ_cent = precipitation centroid "defined as the median of the zonal 
% average precipitation from 20S to 20N. The precipitation is 
% interpolated to a 0.1deg grid over the tropics to allow the precipitation 
% centroid to vary at increments smaller than the grid spacing."
%
%   Nathan Steiger, November 2016
%
%   Significantly modified by N. Steiger June 2017 to account for an 
%   improved definition of the ITCZ location; mods Nov 2017
%

% Create variable metadata
imeta.lon1=lonext(1);imeta.lon2=lonext(2);
%imeta.indxnm='';

% Find nearest lat/lon equal to or outside bounds
ln_a=find((abs(lon-lonext(1))-min(abs(lon-lonext(1))))<1e-5);
ln_b=find((abs(lon-lonext(2))-min(abs(lon-lonext(2))))<1e-5);
% if two equal mins, pick one outside of bounds
l1=min(ln_a);l2=max(ln_b);

% Get indices for extracting precip
if l1>l2
    lon_idx=[l1:length(lon) 1:l2];
elseif l2>l1
    lon_idx=l1:l2;
end

% get tropical precip (using +-30 to account for "monsoon" land areas)
trop_idx=find(lat<= 30 & lat >= -30);
trop_pr=pr(trop_idx,lon_idx,:);

% interpolate to high resolution
tlat=lat(trop_idx);
% spacing = 0.1,
%sp=0.1;x1=min(tlat);x2=max(tlat);
%npts=(sp-x1+x2)/sp;
hhlat=min(tlat):0.1:max(tlat);
%hhlat=linspace(-30,30,601); % 0.1 degree grid in latitude

[~,~,t]=size(trop_pr);

[X,Y] = meshgrid(1:length(lon(lon_idx)),hhlat);
[X1,Y1] = meshgrid(1:length(lon(lon_idx)),tlat);

cdata=zeros(length(hhlat),length(lon(lon_idx)),t);
for j=1:t
   cdata(:,:,j) = interp2(X1,Y1,trop_pr(:,:,j),X,Y);
end

% Take zonal mean
P=squeeze(mean(cdata,2));

% Compute expected latitudes following equation 1 in
% Adam et al 2016 
phi=hhlat(:);

itcz_mx=zeros(length(t),1);
for i=1:t
   itcz_mx(i)=sum(phi.*power(cosd(phi).*P(:,i),10))/sum(power(cosd(phi).*P(:,i),10));
end

itcz_cm=zeros(length(t),1);
for i=1:t
   itcz_cm(i)=sum(phi.*cosd(phi).*P(:,i))/sum(cosd(phi).*P(:,i));
end

% Put variable into the structure
imeta.itcz_mx=itcz_mx;
imeta.itcz_cm=itcz_cm;


end

