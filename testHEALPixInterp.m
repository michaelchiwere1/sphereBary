% Example test function for interpolating on the HEALPiX grid.
clear; clc;
s = spherefun(@(la,th) cosh(sin(cos(la).*sin(th)+5*(cos(la).*sin(th)).*(sin(la).*sin(th)).*cos(th))));

% initializing resolution parameter N_{side}Wk([1 m]) = 1/2 * Wk([1 m]);
err = zeros(11,1); NN = err;
for t = 3:9
Nside = 2^t; n = 4*Nside; m=n-1; fjk = zeros(m,n); lbk=fjk;
% north rings of the function
N = zeros(Nside -1,1); es=zeros(2*Nside+1,1); S = N; 
i=1;
for j = 1:Nside-1
    k = 4*j; v = 2*j;
    th = acos(1 -j^2/(3*Nside^2)); N(j)=th;
    lb =  -pi + pi*((0:k-1)+(1/2))/v;
    fjk(j,1:k) = s(lb,th);
    lbk(j,1:k) = lb;
end
% equator rings of the function
for j = Nside:3*Nside
    k = 4*Nside; v = 2*Nside;
    th = acos(2*(v - j)/(3*Nside));
    es(i) = th; i = i+1;
    lb = -pi + pi*((0:k-1)+mod(j+1,2)/2)/v;
    fjk(j,:) = s(lb,th);
    lbk(j,:)=lb;
end
% southern rings of functions
for j = 1:Nside-1
    k = 4*j; v = 2*j;
    th = acos(-(1 -j^2/(3*Nside^2)));
    lb =  -pi + pi*((0:k-1)+(1/2))/v;
    fjk(n-j,1:k) = s(lb,th);
    S(Nside-j) = th;
    lbk(n-j,1:k) = lb;
end

% latitude grid
thk = [N; es; S];

% Evaluation points over the sphere
N=1000; 
% rng(291901)
lb = (1-2*rand(N,1))*pi;
th = rand(N,1)*pi;

% interpolating using trigonometric weights
%Wk1 = tensor.sphereBaryWeights(thk);
tic
wk = FloaterHormannWght(thk,3);
S = HEALPix.healBaryInterp(lb, th, lbk, thk, fjk,wk);
toc

% Evaluating error of approximations
F = s(lb,th);
%fprintf("Error in using trigonometric weights :%1e\n",norm(S(:) - F(:),inf)./norm(F(:),inf));
err(t) = norm(S(:) - F(:),inf)./norm(F(:),inf);
NN(t) = N;
end

figure(1)
semilogy((1:11)', err,'k+-')
xlabel('t'), ylabel('Error')
set(gca, 'FontSize',14)
set(gcf,'color',[1 1 1]*1)
ylim([1e-15, 1e0 ])
hold off
