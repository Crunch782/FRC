% Main file for Solving 2D NSE w/ Periodic BCs for the KH test
clear all;
clc;

% Set the parameters
param;

% Define the variables
u = zeros(nxm, nym);
v = zeros(nxm, nym);
p = zeros(nxm, nym);
c = zeros(nxm, nym);

% Set the IC
%u = IC(xc,ym,U0,Pj,Ly,Rj,Ax,lx);
c = IC(xm,ym,U0,Pj,Ly,Rj,0.0,lx);

[u,v] = ICRand(nxm,nym,dx,dy,lxmin,lxmax,lymin,lymax);
scale = 0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2);
scale = sqrt(e0./scale);

u = scale*u;
v = scale*v;

fprintf('\nInitial KE = %.3f\n',0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2));


% Compute timestep
%dt = compute_dt(u,v,dx,dy,cfl);
dt=0.5/(1/(dx*dx)+1/(dy*dy))*50;

% Optimize the two Momentum ADIs and the Poisson ADI
% Optimize x ADI
bx = (dt./(2*Re))*(1./dx.^2);
aa = ones(1, nxm);
aa = aa*(-bx);
ab = ones(1, nxm);
ab = ab*(1+2*bx);
ac = ones(1, nxm);
ac = ac*(-bx);

[ax, cx, alpx, xs2x] = ADI_init(aa,ab,ac);

% Optimize y ADI
by = (dt./(2*Re))*(1./dy.^2);
aa = ones(1, nym);
aa = aa*(-by);
ab = ones(1, nym);
ab = ab*(1+2*by);
ac = ones(1, nym);
ac = ac*(-by);

[ay, cy, alpy, xs2y] = ADI_init(aa,ab,ac);

% Optimize poisson ADI
kl = (2./(dx.^2))*(cos(2*pi*(ic-1)./nxm) - 1);

aa = ones(nxm, nym);
aa = aa*(1./(dy.^2));

bb = ones(1, nym);
bb = (1./(dy.^2))*(-2 + (dy.^2)*kl)'*bb;

cc = ones(nxm, nym);
cc = cc*(1./(dy.^2));

aa(1,1) = 0;
bb(1,1) = 1;
cc(1,1) = 0;

[al, cl, alpl, xs2l] = POI_init(aa, bb, cc);

% Main Loop to integrate from t=0 to t=T
maxits = ceil(T./dt);

fprintf('\n\n Beginning Integration . . . \n\n');

[u, v, p, c] = NSE_Solver(u,v,p,c,dt,maxits,ax,cx,alpx,xs2x,ay,cy,alpy,xs2y,al,cl,alpl,xs2l);

fprintf('\n\n Integration Complete! \n\n');

%w = vort(u,v,dx,dy,im,jc,ic,jm);
%plotVort(w,xc,yc);



