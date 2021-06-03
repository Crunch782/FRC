% Main file for solving the heat eqn with FFT

param;
[X,Y] = meshgrid(xc,ym);
X=X';
Y=Y';

% Initialize
u = zeros(nxm, nym);
ue = exactsol(X,Y);
Q = zeros(nxm, nym);

Q = -sourcef(X, Y, 'L');
u = fft(Q);

% Optimize ADI
kl = (2./(dx.^2))*(cos(2*pi*(ic-1)./nxm) - 1);

aa = ones(nxm, nym);
aa = aa*(1./(dy.^2));

bb = ones(1, nym);
bb = (1./(dy.^2))*(-2 + (dy.^2)*kl)'*bb;

cc = ones(nxm, nym);
cc = cc*(1./(dy.^2));

u(1,1) = 0;
aa(1,1) = 0;
bb(1,1) = 1;
cc(1,1) = 0;

[al, cl, alpl, xs2l] = POI_init(aa, bb, cc);
phi = POI_step(al, cl, alpl, xs2l, u);

u = ifft(phi);
u = real(u);

C = ue(1,1) - u(1,1);
u = u + C;

nn = 1:500;
epsn = 1:500;

plotVisuals(u, ue, X, Y, epsnn, nn, 11);

