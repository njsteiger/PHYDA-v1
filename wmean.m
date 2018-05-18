function [M] = wmean(X,lat)
%WMEAN Weighted mean value from global data
%   Given spatial, global or hemispheric data, this function computes the 
%   weighted global mean value defined by: 
%       mean = sum(weights.*data)/sum(weights)
%
%   wmean(X,lat) takes the weighted mean of 'X' according to the given
%   latitude vector 'lat', which is in +/- degrees. 'X' can have any of 
%   the following dimensions: X(lat,lon,time), X(lon,lat,time), 
%   X(lat,lon), X(lon,lat)

if nargin < 2
    error('Not enough input arguments.');
end

lat=lat(:);
[rows,cols,time]=size(X);

% Align X consistently and compute the weighted mean
if time > 1
    if length(lat) ~= rows
        X=permute(X,[2,1,3]);
        cols=rows;
    end
    % Cosine-latitude weighting
    W=cosd(repmat(lat,[1 cols]));
    % Weighted mean
    M=zeros(time,1);
    for i=1:time
        Xi=X(:,:,i);
        M(i)=nansum(W(:).*Xi(:))/nansum(W(:));
    end
else
    if length(lat) ~= rows
        X=permute(X,[2,1]);
        cols=rows;
    end
    % Cosine-latitude weighting
    W=cosd(repmat(lat,[1 cols]));
    % Weighted mean
    M=nansum(W(:).*X(:))/nansum(W(:));
end

end

