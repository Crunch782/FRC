% Main file for solving the heat eqn with implicit method

param;
[X,Y] = meshgrid(xc,ym);
X=X';
Y=Y';

%Initialize
u = zeros(nxm, nym);
du = zeros(nxm,nym);
conv = zeros(nxm, nym);
dt = compute_dt()

% Optimize x ADI
bx = (dt./2)*(1./dx.^2);
a = ones(1, nxm);
a = a*(-bx);
b = ones(1, nxm);
b = b*(1+2*bx);
c = ones(1, nxm);
c = c*(-bx);

[ax, cx, alpx, xs2x] = ADI_init(a,b,c);

% Optimize x ADI
bx = (dt./2)*(1./dx.^2);
a = ones(1, nxm);
a = a*(-bx);
b = ones(1, nxm);
b = b*(1+2*bx);
c = ones(1, nxm);
c = c*(-bx);

[ax, cx, alpx, xs2x] = ADI_init(a,b,c);
% Optimize y ADI
by = (dt./2)*(1./dy.^2);
a = ones(1, nym);
a = a*(-by);
b = ones(1, nym);
b = b*(1+2*by);
c = ones(1, nym);
c = c*(-by);

[ay, cy, alpy, xs2y] = ADI_init(a,b,c);

%Exact Solution
ue = exactsol(X,Y);

%Main Loop
iters = 0;
itersmax = 10000;
t = 0;
epsn = 1;
nn = [1];
epsnn = [L2Norm(dt.*(sourcef(X,Y,'NL')),0)];
fprintf('\n\n Starting Integration \n\n');

while((epsnn(iters+1) > 1e-6)&&(iters <= itersmax))
    
    
    rhs = RHS(u, sourcef(X,Y,'NL'), dt);
    rhs = rhs - (dt./2)*conv;
    conv = nonlinear(u);
    rhs = rhs + (3*dt./2)*conv;
    
    du1 = ADI_step(ax, cx, alpx, xs2x, rhs');
    du = ADI_step(ay, cy, alpy, xs2y, du1');
    
    u = u + du;
    
    espn = L2Norm(du, 0);
    
    iters = iters+1;
    t = t + dt;
    
    nn(iters+1) = iters+1;
    epsnn(iters+1) = espn;
        
end

fprintf('\n\n Integration Complete at time %.2f after %.d iterations \n\n',t,iters);

plotVisuals(u, ue, X, Y, epsnn, nn, 11);