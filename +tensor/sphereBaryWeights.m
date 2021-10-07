
function [Wk,x] = sphereBaryWeights(th)
% SPHEREBARYWEIGHTS
% This funcction computes weights necessary to compute a bivariate
% barycentric interpolation.
% It takes the matrix containing the grid points th(theta) as input and
% returns the barycentric weights for unequally spaced grid e.g. Gauss-Legendre grid.

[m,~] = size(th);

x = cos(th(:,1));
w = x - x';
W = speye(m) + w;
Wk = 1./prod(W,2);

end
