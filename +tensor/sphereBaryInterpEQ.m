function S = sphereBaryInterpEQ(lb,th,lbk, thk, fjk)
% SPHEREBARYINTERPEQ
% Computes a bivariate interpolation of a function defined on a sphere when
% the independent variables are given on an equispaced grid in both latitude and longitude. We require
% that the number of grid points on the latitude is even.
% Usage:
%   S = SPHEREBARYINTERPEQ(lb, th, lbk, thk, fjk)
%       S - column vector of interpolated function.
%       lb, th - column vectors of nongrid points to interpolate at.
%       thk, thk - Matrices of grid values.
%       fjk - Matrix of value of the function on the grid points.



M = numel(th);
[n,J] = size(fjk);

% Condition requiring number of grid points  to be even in longitude
if mod(J,2) ~= 0
    error('The number grid points in longitude must be even');
end
d = J/2;
% computing the known function values fjkplus and fjkminus
fjkp = (1/2)*(fjk(:,1:d) + fjk(:,d+1:end));   % fjkplus
fjkm = (1/2)*(fjk(:,1:d) - fjk(:,d+1:end));   % fjkminus

% Initializing containers for the interpolants in the latitude direction the  g_kplus(th) and g_kminus
gjp = zeros(M,d);   % g_kplus
gjm = zeros(M,d);   % g_kminus

% Precomputing terms to be reused in the interpolation
thj = thk(:,1);
coeff = 1./(cos(th') - cos(thj));
wght1 = (-1).^ (0:n-1)';
wght2 = (-1).^(1:n)';
sinth = sin(th);
sintk = sin(thj);

% coefficients for the even interpolant
coefgd = wght1 .* coeff;
coefgd(1,:) = (1/2) * coefgd(1,:);
coefgd(n,:) = (1/2) * coefgd(n,:);
deno = sum(coefgd,1);

% coefficients for the odd interpolant
temp = wght2 .* (sintk .* coeff);
coefdm  = sintk .* temp;
deno2 = sum(coefdm,1);
%interpolating in the latitude direction

for k= 1:d
    % Evaluating the even interpolant g_kplus
    sjp = fjkp(:,k);
    coefgn = sjp .* coefgd;
    gjp(:,k) = (sum(coefgn,1)./deno).';
    
    % Evaluating the odd interpolant g_kminus
    sjm = fjkm(:,k);
    coefnm = sjm .* temp;
    gjm(:,k) = sinth .* (sum(coefnm,1)./deno2).';
    
    
end
% Evaluation of the whole interpolation by interpolating in the longitude
% direction
S = zeros(M,1);
lbj = lbk(1,:)';
wght = (-1).^(1:d)';

% Interpolant in the longitude direction when number of interpolation
% points in longitude direction is even
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
    % Interpolant in the longitude direction when number of interpolation points
    % is odd
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

