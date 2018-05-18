function [X,lon] = regrid2gcm(Xo,lon_o)
%REGRID2GCM Regrids lat-lon data sets into lons that resemble GCM output
%   Simply switch the bounding box of the data such that lons go from 0 to
%   360 instead of -180 to 180. Xo is in units of (lat,lon,time), and this
%   also assumes that the grid is regular lat-lon.

X=zeros(size(Xo));
lon_mid=size(Xo,2)/2;

% Put last half into first half, then first half into last half
X(:,1:lon_mid,:)=Xo(:,(lon_mid+1):size(Xo,2),:);
X(:,(lon_mid+1):size(Xo,2),:)=Xo(:,1:lon_mid,:);

% Convert lons as well into 0 to 360
lon=zeros(size(lon_o));
nlns=find(lon_o<0);
lon_o(nlns)=lon_o(nlns)+360;

lon(1:lon_mid)=lon_o((lon_mid+1):end);
lon((lon_mid+1):end)=lon_o(1:lon_mid);

end

