% Main file for solving the heat eqaution with implicit method

param;
[X,Y] = meshgrid(xc,ym);
X=X';
Y=Y';

%Initialize
u = zeros(nxm, nym);
du = zeros(nxm,nym);
dt = compute_dt();

%Exact Solution
ue = exactsol(X,Y);

%Main Loop
iters = 0;
itersmax = 10000;
t = 0;
epsn = 1;

nn = [1];
epsnn = [L2Norm(dt.*(sourcef(X,Y,'L')),0)];
fprintf('\n\n Starting Integration \n\n')

while((epsnn(iters+1) > 5e-5)&&(iters <= itersmax))
   
    du = dt.*(sourcef(X,Y,'L') + NSE_calc_lap(u));
    u = u + du;
    
    espn = L2Norm(du, 0);
    
    iters = iters+1;
    t = t + dt;
    
    nn(iters+1) = iters+1;
    epsnn(iters+1) = espn;
        
end

fprintf('\n\n Integration Complete at time %.2f after %.d iterations \n\n',t,iters);

plotVisuals(u, ue, X, Y, epsnn, nn, 11);



