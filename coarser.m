function [cdata] = coarser(olat,olon,odata,clat,clon)
% COARSER Interpolation of model data onto coarser resolution
%   function [cdata] = coarser(olat,olon,odata,clat,clon)
%       where olat/olon are the finer original resolution lats/lons,
%       odata is the original data [lat,lon,time], and clat/clon are
%       coarser lats/lons.
%
%   Interpolates model data to desired coarser resolution while taking
%   account of possible mismatches in the end longitude values. For
%   example, if we take olon and clon values to be
%   olon = [0.5, 1, 1.5, ... 358.5, 359,359.5]
%   clon = [0, 2, ... 356, 358]
%   there will be a missing olon value at 0 deg for interpolation purposes.

[~,~,t]=size(odata);

%olat1=olat(:); olon1=olon(:);
%clat1=clat(:); clon1=clon(:);
olat=olat(:); olon=olon(:);
clat=clat(:); clon=clon(:);

% Make finer instead of coarser if lat/lons require
%if length(clat)>length(olat) || length(clon)>length(olon)
%   clat=olat1; clon=olon1;
%   olat=clat1; olon=clon1;
%   disp('Making data finer, not coarser...')
%else % Just coarsen like normal
%   clat=olat1; clon=olon1;
%   olat=clat1; olon=clon1;
%end

if min(clon) > 10 && min(olon) > 10
    disp('COARSER: Assuming that the data are not full lat/lon grids...')
end

if min(clon) < min(olon) && min(clon) < 10 && min(olon) < 10 % Kludge for grids that are partial, not full grids
    
    % add value on end
    clon(end+1)=clon(end)+clon(2)-clon(1);
    
    % repeat end value in higher resolution data
    odata(:,end+1,:)=odata(:,1,:);
    
    [X,Y] = meshgrid(clon,clat);
    [X1,Y1] = meshgrid(olon,olat);
    
    cdata=zeros(length(clat),length(clon),t);
    for j=1:t
        cdata(:,:,j) = interp2(X1,Y1,odata(:,:,j),X,Y);
    end
    
    % replace end value at beginning, then get rid of extra end
    cdata(:,1,:)=cdata(:,end,:);
    cdata(:,end,:)=[];
    
    % interpolation seems to flip latitudes...
    % cdata=flipdim(cdata,1);
    
else
    
    [X,Y] = meshgrid(clon,clat);
    [X1,Y1] = meshgrid(olon,olat);
    
    cdata=zeros(length(clat),length(clon),t);
    for j=1:t
        cdata(:,:,j) = interp2(X1,Y1,odata(:,:,j),X,Y);
    end
    
    % interpolation seems to flip latitudes...
    % cdata=flipdim(cdata,1);
    
end

end

