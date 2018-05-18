function [M] = wmean_a(X,A)
%WMEAN_A Weighted mean value from global data
%   Given spatial data, this function computes the weighted global mean 
%   value defined by:
%       mean = sum(weights.*data)/sum(weights)
%
%   wmean_a(X,A) takes the weighted mean of X according to the given
%   area matrix A, which is the area of the grid cells of X.

if nargin < 2
  error('Not enough input arguments.');
end

[rows1,cols1]=size(A); [rows,cols,time]=size(X);
if rows1~=rows || cols1~=cols
  error('X is not compatible with A.')
end

% Compute the weighted mean
M=zeros(time,1);
for i=1:time
  Xi=X(:,:,i);
  % Make weights have nans where X has them
  W=A; W(isnan(Xi))=nan;
  M(i)=nansum(W(:).*Xi(:))/nansum(W(:));
end

end

