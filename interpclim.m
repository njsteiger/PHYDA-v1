function [cdata] = interpclim(olat,olon,odata,clat,clon)
% INTERPCLIM Interpolation of model data onto different resolution
%   function [cdata] = interpclim(olat,olon,odata,clat,clon)
%       where olat/olon are the original resolution lats/lons,
%       odata is the original data [lat,lon,time], and clat/clon are
%       destination lats/lons. Note that this function doesn't account
% 	for grid mismatches at the end points. It just assumes that
%	both grids extend to the same lat/lon extents.
%	
%   Nathan Steiger, July 2017

[~,~,t]=size(odata);

ot=linspace(0,1,length(olat));
on=linspace(0,1,length(olon));
ct=linspace(0,1,length(clat));
cn=linspace(0,1,length(clon));
    
[X,Y] = meshgrid(cn,ct);
[X1,Y1] = meshgrid(on,ot);

cdata=zeros(length(clat),length(clon),t);
for j=1:t
  cdata(:,:,j) = interp2(X1,Y1,odata(:,:,j),X,Y);
end
    

end

