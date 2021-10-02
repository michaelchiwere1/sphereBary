%HEALBARYINTERP
% Implements a bivariate barycentric formula for interpolating on the sphere
% for interpolation data given at HEALPx points on a sphere.
% Usage:
%    S = HEALBARYINTERP(lb, th, lbk, thk, fjk, WK).
%       lb, th : column vectors of evaluation points.
%       lbk, thk : matrice of HEALPix points.
%       fjk : function values at the grid values.
%       Wk : barycentric weight for interpolating in the latitude direction

function S = healBaryInterp(lb, th, lbk, thk, fjk, Wk)
M = numel(th); S=zeros(M,1);
[~,J] = size(fjk);
% Condition requiring points in the even direction be even
if mod(J,2) ~= 0
    error('The number of grid points in longitude must be even');
end
d = J/2; Nside = J/4; fjkp = cell(J-1,1); fjkm = cell(J-1,1); lbj = cell(J-1,1);

% Split the known function into fjkplus and fjkminus
for j=1:Nside-1
    k = 2*j;
    % north rings of the function
    fjkp{j} = (1/2)*(fjk(j,1:k) + fjk(j,k+1:2*k));
    fjkm{j} = (1/2)*(fjk(j,1:k) - fjk(j,k+1:2*k));
    lbj{j} = lbk(j,k+1:2*k);
    % southern rings of functions
    fjkp{J-j} = (1/2)*(fjk(J-j,1:k) + fjk(J-j,k+1:2*k));
    fjkm{J-j} = (1/2)*(fjk(J-j,1:k) - fjk(J-j,k+1:2*k));
    lbj{J-j} = lbk(J-j,k+1:2*k);
end
% equatorial rings of the function
for j=Nside:3*Nside
    fjkp{j} = (1/2)*(fjk(j,1:d) + fjk(j,d+1:end));
    fjkm{j} = (1/2)*(fjk(j,1:d) - fjk(j,d+1:end));
    lbj{j} = lbk(j,d+1:end);
end
% interpolation step (begin with variable lambda)
gkplus = zeros(J-1,M); gkminus = gkplus;

for j=1:J-1
    wght = (-1).^(1:length(fjkp{j}))';
    denom = wght .* cot(lb'-lbj{j}');
    numplus = denom .* fjkp{j}';
    numminus = wght .* csc(lb'-lbj{j}') .* fjkm{j}';
    gkplus(j,:) = (sum(numplus,1)./sum(denom,1)).';
    gkminus(j,:) = (sum(numminus,1)./sum(denom,1)).';
end
% interpolating in the theta direction
csth = cos(th); csthk = cos(thk); sth= 1./sin(thk);
for k=1:M
    cdiff = csth(k) - csthk;
    denom = Wk./cdiff;
    numplus = denom .* gkplus(:,k);
    numminus = sin(th(k))* (sth .* denom .* gkminus(:,k));
    S(k) = (sum(numplus,1)-sum(numminus,1))./sum(denom,1);
    
end
end