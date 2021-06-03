% Ordered fft2
% f: inverse fourier transform of g
% x: vector of space along k
% y: vector of space along h


function [f x y] = oifft2(g,k,h)

    [nh, nk] = size(k);
    dk = diff(k(1,:)); dk=dk(1,1);
    dh = diff(h(:,1)); dh=dh(1,1);

    dx = 2*pi*1/(nk*dk);
    xx = dx*(-nk/2:nk/2-1);

    dy = 2*pi*1/(nh*dh);
    yy = dy*(-nh/2:nh/2-1);

    [x, y] = meshgrid(xx,yy);

    g = fftshift(g);
    f = (nk*nh)*ifft2(g);

end