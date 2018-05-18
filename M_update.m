
function [Xa,Xa_m]=M_update(Xb,HXb,y,R,infl)
% function [Xa,Xa_m]=M_update(Xb,HXb,y,R,infl)
% M_UPDATE  Matrix update: Data assimilation of all observations for a given time step.
%	Note that this function does not include covariance localization. The update 
% 	equations are taken from Whitaker and Hamill 2002: Eq. 2, 4, 10, & Sec 3. 
%
%   INPUT:
%       Xb = background (prior) [statevector ens]
%       y = observation [vector]
%       HXb = model estimate of observations H(X_b) [length(y) ens]
%       R = observation error covariance matrix [length(y) length(y)]
%       infl = inflation factor, almost always = 1 [scalar]
%
%   OUTPUT:
%       Xa = analysis (posterior) [statevector ens]
%       Xa_m = analysis mean [statevector]
%
%   Nathan Steiger, Columbia University, July 2017
%

% Size of ensemble
[~,n_ens]=size(Xb);

% Decompose X_b and HXb (for Eqs. 4 & Sec 3)
Xb_m=mean(Xb,2);
Xb_p=bsxfun(@minus,Xb,Xb_m);

HXb_m=mean(HXb,2);
HXb_p=bsxfun(@minus,HXb,HXb_m);

% Apply inflation (if chosen to be ~= 1)
if infl~=1; Xb_p=infl*Xb_p; end

% Kalman gain for mean and matrix covariances (Eq. 2)
PbHT=(Xb_p*HXb_p')./(n_ens-1);
HPbHTR=(HXb_p*HXb_p')./(n_ens-1)+R;
K=PbHT*inv(HPbHTR);

% Kalman gain for the perturbations (Eq. 10)
sHPbHTR=sqrtm(HPbHTR);
%sR=sqrtm(R); % Use for non-diagonal R matrix (much slower)
sR=sqrt(R); % Assumes diagonal R
Kt_n=PbHT*(inv(sHPbHTR))';
Kt_d=inv(sHPbHTR+sR);
Kt=Kt_n*Kt_d;

% Update mean and perturbations (Eq. 4 & Sec 3)
Xa_m=Xb_m+K*(y-HXb_m);
Xa_p=Xb_p-Kt*HXb_p;

% Reconstitute the full prior
Xa=bsxfun(@plus,Xa_p,Xa_m);



