% Ordered fft2
% g: fourier transform of f
% k: vector of frequencies along x
% h: vector of frequencies along y

function [g, k, h] = offt2(f,x,y,nx,ny,dx,dy)

    dk = 2*pi*1/(nx*dx);
    kk = dk*(-nx/2:nx/2-1);

    dh = 2*pi*1/(ny*dy);
    hh = dh*(-ny/2:ny/2-1);

    [k, h] = meshgrid(kk,hh);
    k=k';
    h=h';
    g = fft2(f)./(nx*ny);
    g = fftshift(g);

end