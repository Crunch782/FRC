% Parameters for the Problem

global Lx Ly
global Vol
global nx ny nxm nym
global dx dy
global xm ym xc yc
global a b
global ip im jp jm ic jc
global cfl
global epsilon

Lx = 1;
Ly = 2;
Vol = Lx*Ly;

nx = 65;
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

cfl = 100;
epsilon = 5e-5;