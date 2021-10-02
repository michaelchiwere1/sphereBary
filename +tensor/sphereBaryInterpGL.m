% SPHEREBARYINTERPGL
% Implements a bivariate barrycentric formula for interpolating a
% function on a sphere. Here we assume the data is given on an equally spaced
% grid in lambda (lbk) consisting of an even number points and we prefer the
% points in latitude form a Gauss-Legendre grid.
% Usage:
%   S = SPHEREBARYINTERPGL(lb, th, lbk, thk, fjk, Wk)
%       S - column vector of interpolated function.
%       lb, th - column vectors of nongrid points to interpolate at.
%       thk, lbk - Matrices of grid values.
%       fjk - Value of the function on the grid points.
%       Wk - Weights for interpolating in the theta direction. Weights are
%       computed as Wk = sphereBaryWeights(thk).


function S = sphereBaryInterpGL(lb,th,lbk,thk,fjk,Wk)
M = numel(th);
[~,J] = size(fjk);
% A condition ensuring the number of grid points in longitude is even
if mod(J,2) ~= 0
    error('The number grid points in longitude must be even');
end
d = J/2;
% computing the known function values fjkplus and fjkminus
fjkp = (1/2)*(fjk(:,1:d) + fjk(:,d+1:end));   % fjkplus
fjkm = (1/2)*(fjk(:,1:d) - fjk(:,d+1:end));   % fjkminus

% Initializing containers the  g_kplus(th) and g_kminus for the interpolants in the latitude direction
gjp = zeros(M,d);   % g_kplus
gjm = zeros(M,d);   % g_kminus

% Precomputing values to be reused in interpolation
thj = thk(:,1);
coeff = cos(th') - cos(thj);
coefdm =  Wk./coeff;
denom = sum(coefdm,1);
if(ismember(0,thj)||ismember(pi,thj))
    coefdm2 = coefdm.*sin(thj);
    sinth = 1./sin(th);
else
    coefdm2 = coefdm./sin(thj);
    sinth = sin(th);
end

% Evaluating the interpolants in latitude direction
for j= 1:d
    
    % evaluating the even interpolant in theta (g_kplus)
    sjp = fjkp(:,j);
    coefgn = sjp .* coefdm;
    gjp(:,j) = (sum(coefgn,1)./denom).';
    
    % evaluating the odd interpolant in theta (g_kminus)
    sjm = fjkm(:,j);
    coefnm = sjm .* ( coefdm2 );
    gjm(:,j) = sinth.*(sum(coefnm,1)./denom).';
    
end

% Evaluating the  whole formula by interpolating in the longitude direction
S = zeros(M,1);
lbj = lbk(1,:)';
n = length(lbj);
wght = (-1).^(1:n)';

% Interpolant in longitude when number of grid points d is even
if(mod(n,2)== 0)
    for i = 1:M
        diff = lb(i) - lbj;
        cdiff = cot(diff);
        denom = wght .* cdiff;
        denum = cdiff .* gjp(i,:)' - csc(diff) .* gjm(i,:)';
        num = wght .* denum;
        S(i) = sum(num)./sum(denom);
    end
    
else
% Interpolant in longitude when number of grid points d is odd 
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
