clear; clc;
% Test function for the bivariate interpolation formula
% on the Gauss-Legendre (GL) grid.
s = spherefun(@(la,th) cosh(sin(cos(la).*sin(th)+50*(cos(la).*sin(th)).*(sin(la).*sin(th)).*cos(th))));

n = 400; m = 2*n;
% Full lat-lon grid
thk = acos(legpts(n+1)); 
lbj =  -pi+(2*pi/m)*(0:m-1); 
[lbj, thj] = meshgrid(lbj,thk);
% Data at the grid
fjk = s(lbj,thj);
% lat-lon grid where longitude is restricted to [0 pi].
lbk = lbj(:,floor(m/2)+1:end); thk = thj(:,floor(m/2)+1:end);
% Barycentric weights in theta
Wk = tensor.sphereBaryWeights(thk);
% Evaluation points over the sphere
N=m*n;
lb = (1-2*rand(N,1))*pi;
th = rand(N,1)*pi;
tic
S = tensor.sphereBaryInterpGL(lb,th,lbk,thk,fjk,Wk);
toc
% Error
F = s(lb,th);
norm(S(:) - F(:),inf)/norm(F(:),inf)
