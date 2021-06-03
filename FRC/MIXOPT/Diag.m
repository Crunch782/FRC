% Computes the diagnostics for the optimal solution

format long;
set(0,'DefaultFigureVisible','off')
param;

fprintf('\n\n Starting Diagnostics for (s,T,Re) = (%.1f, %.d, %.1f) \n\n',s,T,Re);

% Convert to fields
Xopt = load(strcat('IC/X_',tag,'.dat'),'Xopt','-ASCII');
u = reshape(Xopt(1:nxm*nym),nxm,nym);
v = reshape(Xopt(nxm*nym+1:end),nxm,nym);
w = vort(u,v,dx,dy,im,jc,ic,jm);
vl = floor(min(min(w)));
vu = ceil(max(max(w)));

% Set the Passive Scalar
c = ICScalar(xm,ym,dx,dy,Vol);
cd = ICScalar(xm,ym,dx,dy,Vol);

% Now initialize Pressure
p = ICPress(u,v,dx,dy,im,ic,ip,jm,jc,jp,al,cl,alpl,xs2l);

% Check Normalization
fprintf('\nInitial KE = %.2f\n',0.5*(L2Norm(u,1).^2 + L2Norm(v,1).^2));

%Create directory for the diagnostic data
mkdir(['Data/' tag '/']);
mkdir(['Data/' tag '/Figures/']);
mkdir(['Data/' tag '/Movies/']);


% Main Loop to integrate from t=0 to t=T
[u, v, p, c] = NSE_SolverDiag(u,v,p,c,cd,vl,vu);

fprintf('\n\n Diagnostics Complete! \n\n');



