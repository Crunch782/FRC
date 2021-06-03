% Parameters for the Problem

global Lx Ly
global Vol
global nx ny nxm nym
global dx dy
global xm ym xc yc
global T nt dt Td nd nf
global a b
global ip im jp jm ic jc
global cfl
global Re Pe
global e0
global lxmax lxmin lymax lymin
global s
global ax cx alpx xs2x
global ay cy alpy xs2y
global al cl alpl xs2l
global tag

Lx = 2*pi;
Ly = 2*pi;
Vol = Lx*Ly;

nx = 129;
ny = 129;
nxm = nx-1;
nym = ny-1;

dx = Lx./nxm;
dy = Ly./nym;

ic = 1:nxm;
jc = 1:nym;

xc = (ic-1)*dx;
yc = (jc-1)*dy;

xm = (ic-1./2)*dx;
ym = (jc-1./2)*dy;

T = 1;
dt = 0.002;
nt = T./dt;
Td = 1;
nd = Td./dt;
nf = floor(nd./(Td*24));

a = 2*pi./Lx;
b = 2*pi./Ly;

ip = 1:nxm;
ip = ip + 1;
ip(nxm) = 1;
im = 1:nxm;
im = im - 1;
im(1) = nxm;

jp = 1:nym;
jp = jp + 1;
jp(nym) = 1;
jm = 1:nym;
jm = jm- 1;
jm(1) = nym;

Re = 50;
Pe = 50;
e0 = 3e-2;

cfl = 0.2;

lxmax = 2;
lxmin = 1;
lymax = 2;
lymin = 1;

s = 5;

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

tag = append(num2str(s),'_',num2str(T),'_',num2str(Re));
