function S = sphereBaryInterpSEQ(lb,th,lbk, thk, fjk)
% SPHEREBARYINTERPSEQ
% Implementation of a barycentric interpolation of a unit sphere where the
% data is given at mid points of the grids on both variables,
% i.e theta (thk) = ((0:n-1)+1/2)*pi/n, lambda(thk) = -pi+2pi/m*((0:m-1)+1/2).
% Usage:
%       S = SPHEREBARYINTERPSEQ(lb,th,lbk,thk, fjk)
%           S - column vector of interpolated function.
%           lb, th - column vectors of nongrid points to interpolate at.
%           thk, lbk - Matrices of grid values.
%           fjk - Matrix of value of the function on the grid points.

% Author: Michael Chiwere

M = numel(th);
[n,J] = size(fjk);

% A condition requiring number of points in latiude direction be even
if mod(J,2) ~= 0
    error('The number grid points in longitude must be even');
end
d = J/2;
% computing the known function values fjkplus and fjkminus
fjkp = (1/2)*(fjk(:,1:d) + fjk(:,d+1:end));   % fjkplus
fjkm = (1/2)*(fjk(:,1:d) - fjk(:,d+1:end));   % fjkminus


% Initializing container for interpolant in the latitude direction:  g_jplus(th) and g_jminus
gjp = zeros(M,d);   % g_jplus
gjm = zeros(M,d);   % g_jminus

% Precomputing values o be reused
thj = thk(:,1);
wght = (-1).^(0:n-1)';
costh = 1./(cos(th') - cos(thj));
sintj = sin(thj);
sinth = sin(th);
wght2 = (-1).^(1:n)';

% Coefficients for the even interpolant
denomt = wght .* (sintj .* costh);
den = sum(denomt,1);
% Coefficients for the odd interpolant
temp = wght2 .*   costh;
coefdm  = sintj.* temp;
den2 = sum(coefdm,1);

% evaluating interpolant in the longitude direction
for j= 1:d
    % Even interpolant in theta; g_kplus
    sjp = fjkp(:,j);
    num = sjp .* denomt;
    gjp(:,j) = (sum(num,1)./den).';
    
    
    % Odd interpolant in the theta; gkminus
    sjm = fjkm(:,j);
    coefnm = sjm .* temp;
    gjm(:,j) = sinth .* (sum(coefnm,1)./den2).';
end

% Evaluation the whole formula by interpolating in the longitude direction
S = zeros(M,1);
lbj = lbk(1,:)';
wght = (-1).^(1:d)';

% Interpolant when number of grid points in longitude d is even
if(mod(d,2)== 0)
    for i = 1:M
        diff = lb(i) - lbj;
        cdiff = cot(diff);
        denom = wght .* cdiff;
        denum = cdiff .* gjp(i,:)' - csc(diff) .* gjm(i,:)';
        num = wght .* denum;
        S(i) = sum(num)./sum(denom);
    end
    
else
    % Interpolant when number of grid points in longitude d is odd
    for i = 1:M
        diff = lb(i) - lbj;
        cdiff = csc(diff);
        denom = wght .* cdiff;
        denum = cdiff .* gjp(i,:)' - cot(diff) .* gjm(i,:)';
        num = wght .* denum;
        S(i) = sum(num,1)./sum(denom,1);
    end
end
end


