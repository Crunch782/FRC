% Gaussian Noise

function f = noise2(nx,ny,dx,dy,lxmin,lxmax,lymin,lymax)

% Compute vector in phase space
dk = 2*pi*1/(nx*dx);
dh = 2*pi*1/(ny*dy);
kk = dk*(-nx/2:nx/2-1);
hh = dh*(-ny/2:ny/2-1);
[k, h] = meshgrid(kk,hh);

% Random amplitude and phase
ra = 0.5*rand(ny,nx).*exp(1i*2*pi*rand(ny,nx));

% Gaussian shape in 2d... 
kmax = 2*pi/lxmin;
kmin = 2*pi/lxmax;
hmax = 2*pi/lymin;
hmin = 2*pi/lymax;
kmean = (kmin+kmax)/2;
hmean = (hmin+hmax)/2;
sk = (kmax-kmin)/2;
sh = (hmax-hmin)/2;
if sk > sh
    skh = min(sk,sh) + max(sk,sh)*(k./sqrt(k.^2+h.^2+eps)).^2;
else
    skh = min(sk,sh) + max(sk,sh)*(h./sqrt(k.^2+h.^2+eps)).^2;    
end
r2 = ((k/kmean).^2+(h/hmean).^2);
gs = exp(-1/2*((r2-1)./skh).^2);

% Hermitian conjugate for real signal
% TO FIX: THIS GIVES A IMAGINARY PART AS WELL...   %%%%Florence :: did he
% fix this ??
gp = ra.*gs;
gmx = [zeros(1,nx) ; fliplr(gp(2:end,:))];
gmy = [zeros(ny,1) flipud(gp(:,2:end))];
gm = gmx + gmy;
g = gm + conj(gp);


% Normalization for int(g) = int(f) = 1
g = g/sqrt(g(:)'*g(:)*dh*dk/((kk(end)-kk(1))*(hh(end)-hh(1))));

% Inverse fourier transform
[f, ~] = oifft2(g,k,h);

% Normalization
f = f ./ sqrt(nx*ny);

% Real vector  
f = real(f);