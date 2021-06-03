% Main file for the Optimization Problem

format long;
clc;
clear all;

tic
param;

% Define the variables
u = zeros(nxm, nym);
v = zeros(nxm, nym);
p = zeros(nxm, nym);
c = zeros(nxm, nym);

% Set the IC and Normalize
[u,v] = ICRand(nxm,nym,xc,ym,xm,yc,dx,dy,lxmin,lxmax,lymin,lymax);
scale = 0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2);
scale = sqrt(e0./scale);
u = scale*u;
v = scale*v;

% Convert Fields
Xg = [u(:) ; v(:)];

% Start Optimization Algorithm
fprintf('\n\n Starting Optimization Algorithm . . . \n\n');
mkdir(['Checkpoints/' tag '/']);
[JDJ, Xopt, neval, conv]  = gboptim('J_eval',Xg,-1,'conjgrad',norm(Xg),2,tag);
save(strcat('IC/X_',tag,'.dat'),'Xopt','-ASCII');
toc

fprintf('\n\n Optimization Algorithm Complete! \n\n');

if conv == 1
    fprintf('\n\n Convergence Reached After %d Loops! \n\n',neval);
else
    fprintf('\n\n Solution not converged after %d  loops . . . \n\n',neval);
end

