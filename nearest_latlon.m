function [nlm,nlat,nlon] = nearest_latlon(mlat,mlon,plat,plon,varargin)
%NEAREST_LATLON Find the nearest lat and lon in model space given vectors
%   or grids of lat/lon data. 'mlat/mlon' are from the model fields,
%   'plat/plon' are the points that we want to query. The function finds
%   the closest spherical distance, then selects those points that are
%   nearest the query points. 'nlm' is the nearest location in the model
%   state vector (the vector defined by lat*lon, not [lat,lon]), and
%   'nlat/nlon' are the nearest lat/lon values in the model space. Also
%   includes optional search over non-masked locations when input includes
%   'mask' matrix/vector that contains NaNs over masked points.
%
%   Nathan Steiger, LDEO, Sept 2016


if nargin > 4
    
    if isvector(mlat) && isvector(mlon)
        mlat=mlat(:);mlon=mlon(:);
        qlat=repmat(mlat,length(mlon),1);
        qlon=reshape(repmat(mlon,1,length(mlat))',length(mlon)*length(mlat),1);
    else
        qlat=mlat(:);qlon=mlon(:);
    end
    
    %----------------------
    % Account for NaN mask
    %----------------------
    mask=varargin{1}(:); % vectorize mask
    % Find the non-nan-masked locations to do distance calculation over
    d_ids=find(~isnan(mask));
    qnlat=qlat(d_ids);
    qnlon=qlon(d_ids);
    
        
    nlm=zeros(length(plat),1);
    nlat=zeros(length(plat),1);
    nlon=zeros(length(plat),1);
    for i=1:length(plat)
        % Great circle distances (Haversine formula)
        A=sind((qnlat-plat(i))/2).^2+cosd(plat(i)).*cosd(qnlat).*...
            sind((qnlon-plon(i))/2).^2;
        Di=6371.0*2*atan2(sqrt(A),real(sqrt(1-A)));
        
        D=NaN(size(mask));
        D(d_ids)=Di; % put un-masked values into full variable
        
        % Find minimum distance (omits nans)
        m_id=find(D==min(D));
        
        nlm(i)=m_id(1); % pick first if there are equivalent values
        % pull out lat/lon from original, unmasked values
        nlat(i)=qlat(nlm(i)); 
        nlon(i)=qlon(nlm(i));
    end
    
    
    
    
else
    % Assume all grid points are available
    if isvector(mlat) && isvector(mlon)
        mlat=mlat(:);mlon=mlon(:);
        qlat=repmat(mlat,length(mlon),1);
        qlon=reshape(repmat(mlon,1,length(mlat))',length(mlon)*length(mlat),1);
    else
        qlat=mlat(:);qlon=mlon(:);
    end
    nlm=zeros(length(plat),1);
    nlat=zeros(length(plat),1);
    nlon=zeros(length(plat),1);
    for i=1:length(plat)
        % Great circle distances (Haversine formula)
        A=sind((qlat-plat(i))/2).^2+cosd(plat(i)).*cosd(qlat).*...
            sind((qlon-plon(i))/2).^2;
        D=6371.0*2*atan2(sqrt(A),real(sqrt(1-A)));
        
        % Find minimum distance
        m_id=find(D==min(D));
        
        nlm(i)=m_id(1); % pick first if there are equivalent values
        nlat(i)=qlat(nlm(i));
        nlon(i)=qlon(nlm(i));
    end
    
end

end

