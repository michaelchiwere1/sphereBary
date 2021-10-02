%FLOATERHORMANNWGHT     computes weights necessary evaluate a rational
% trigonometric interpolating polynomial of the Floater-Hommann class.
% Takes in matrix or a column vector of node thk, and blending parameter d.
% the weights computed are prod(1/(cos(x_k)-cos(x_j))).
% Usage:
%   W_k = FLOATERHORMANNWGHT(thk,d)
%       d - a blending parameter, an integer.
%       thk - matrix of grid values in longitude direction.
%       Wk - weights necessary to interpolate in the latitude direction.


function Wk = FloaterHormannWght(thk, d)
%get length of the input
[n,~] = size(thk);
x = cos(thk(:,1));
%enforce necessary conditions
if(n<=d)
    error("d too large for the problem size");
end
% Periodically extend the nodes by d points < 0 and > 2*pi
x = [x(n-d+1:n)-2*pi;x;x(1:d)+2*pi];
% Number of points has increased by 2*d
n = n + 2*d;

Wk = zeros(n,1);
% Evaluating the weights one by one while enforcing necessary conditions
for k=1:n
    imin = max(k-d,1);
    if(k>=n-d);imax=n-d;else;imax=k;end
    if(mod(imin,2));temp=1;else;temp=-1;end
    
    sum = 0;
    for i=imin:imax
        jmax=min(i+d,n);
        coeff = 1.0;
        for j=i:jmax
            if(j~=k)
                coeff = coeff * (x(k)- x(j));
            end
        end
        coeff = temp/coeff;
        temp = -temp;
        sum = sum + coeff;
    end
    Wk(k) = sum;
    
end
% Just take the weights corresponding to the original values over [0 pi)
Wk = Wk(d+1:n-d);
end