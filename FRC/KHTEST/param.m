% Parameters for the Problem

global Lx Ly
global Vol
global nx ny nxm nym
global dx dy
global xm ym xc yc
global T 
global a b
global ip im jp jm ic jc
global cfl
global Re Pe
global e0
global U0 Pj Rj lx Ax
global lxmax lxmin lymax lymin

Lx = 2;
Ly = 1;
Vol = Lx*Ly;

nx = 257;
ny = 257;
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

Re = 500;
Pe = 500;
e0 = 3e-2;

cfl = 0.2;

U0 = 1;
Pj = 20;
Rj = Ly./4;
Ax = 0.25;
lx = 0.5*Lx;
lxmax = Lx;
lxmin = 0.1*Lx;
lymax = Ly;
lymin = 0.1*Ly;
